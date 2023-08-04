# Purpose: make ribbon plot of cell state vs identity
# Author: Salomé Carcy
# Date: May 2023


# **************
# 1. IMPORT ####
# **************

library(Seurat)
library(GGally)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(ggalluvial)
source("./scripts-final/colors_universal.R")


# Import data
seur <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
DimPlot(seur, group.by = "new_clusters", label=T, repel=T, reduction="UMAP_50")+
  scale_color_manual(values=cols_integrated)


# Import GEPs
gep_topgenes <- read.csv("./data/human-thymus/HumanData_17_GEPsOnParkData/genes_per_GEP_df_2023-04-07.csv", row.names=1)
dim(gep_topgenes) # 12 GEPs

# table(gep_topgenes$GEP_1[1:200] %in% gep_topgenes$GEP_5[1:200])




# ********************
# 2. COMPUTE GEPs ####
# ********************

gep_list <- list()
for (i in 1:ncol(gep_topgenes)){
  print(paste0("GEP", i))
  gep_allgenes <- gep_topgenes[1:200,i][!is.na(gep_topgenes[1:200,i])]
  # gep_allgenes <- na.omit(gep_topgenes[,i])
  print(length(gep_allgenes))
  gep_list[[colnames(gep_topgenes)[i]]] <- gep_allgenes
}

seur.geps     <- AddModuleScore(seur, name = "GEP", features=gep_list)
# seur.geps <- readRDS("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/seuratobj_gepscores_allgenes.rds")

# Sanity check
SCpubr::do_FeaturePlot(seur.geps, reduction="UMAP_50", features=paste0("GEP", 1:length(gep_list)), ncol=6,
                       viridis_color_map = "B", order=T)
# ggsave("./scripts-in-progress/human-thymus/HumanData_20_RibbonPlotCellStateToID/plots/gapin_cNMF_geps_allgenes.jpeg", width=40, height=20)




# *************************************
# 3. THRESHOLD CELLS BASED ON GEPs ####
# *************************************

gep_pbmc <- c("GEP1", "GEP4", "GEP5", "GEP6", "GEP8", "GEP11", "GEP12")
df <- as.data.frame(seur.geps@meta.data[seur.geps@meta.data$Tissue=="PBMC",gep_pbmc])
# Check the distribution
ggplot(pivot_longer(df, cols=gep_pbmc, names_to="gep", values_to="score"))+
  geom_histogram(aes(x=score), bins = 100)+
  facet_wrap(~gep, ncol=4)+
  # _______
  # geom_violin(aes(x=gep, y=score), width=1)+
  # geom_boxplot(aes(x=gep, y=score), outlier.shape = NA, width=0.05)+
  # geom_jitter(aes(x=gep, y=score), size=0.1, width = 0.05)+
  # _______
  labs(y="GEP score", title="Raw GEP score")


## 3.5. Method 5: Use AUCell to do a gene set (GEP) enrichment ####
# BiocManager::install("AUCell")
library(AUCell)

# Get gene sets
gep_fulllist <- list()
for (i in 1:ncol(gep_topgenes)){
  print(paste0("GEP", i))
  gep_allgenes <- gep_topgenes[,i][!is.na(gep_topgenes[,i])]
  print(length(gep_allgenes))
  gep_fulllist[[colnames(gep_topgenes)[i]]] <- gep_allgenes
}
# gep_fulllist <- gep_list

# Define expression matrix (features as rows, cells as columns) and gene sets to score
seur.pbmc <- subset(seur, subset=Tissue=="PBMC")
exprMatrix <- seur.pbmc@assays$RNA@counts
geneSets <- gep_fulllist[gsub("^(GEP)(.*)$", "\\1_\\2",  gep_pbmc)]
# names(geneSets)
# Add a "negative control" of 504 random genes
# set.seed(123)
# geneSets$random <- sample(rownames(exprMatrix), 504)
# lapply(geneSets, function(x) length(x))

# cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=TRUE) # there are at least 489 genes detected per cell

# Calculate enrichment scores
cells_AUC <- AUCell_run(exprMatrix, geneSets, aucMaxRank=504) # checking if genes from gene set are present within top 504 most expressed genes per cell

# # Optional: Set the assignment thresholds
par(mfrow=c(3,3))
set.seed(123)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=4, assign=TRUE)

# Optimize the thresholds
cellsUMAP <- seur.pbmc[["UMAP_50"]]@cell.embeddings
selectedThresholds <- getThresholdSelected(cells_assignment)
jpeg("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/maxrank504_thresholds_manual2.jpeg", width=2000, height=5000, res=250)
par(mfrow=c(7,3))
AUCell_plotTSNE(tSNE=cellsUMAP, exprMat=exprMatrix,
                cellsAUC=cells_AUC, thresholds=selectedThresholds)
dev.off()
selectedThresholds[1] <- 0.05 # GEP1 threshold
selectedThresholds[2] <- 0.08 # GEP4 threshold
selectedThresholds[3] <- 0.35 # GEP5 threshold
selectedThresholds[4] <- 0.11 # GEP6 threshold
selectedThresholds[5] <- 0.09 # GEP8 threshold
selectedThresholds[6] <- 0.10 # GEP11 threshold
selectedThresholds[7] <- 0.10 # GEP12 threshold

# Update cell assignment
cells_assignment_new <- AUCell_assignCells(cells_AUC, thresholds=selectedThresholds)
# Check by heatmap cell assignment
assignmentTable <- reshape2::melt(lapply(cells_assignment_new, function(x) x$assignment), value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
# head(assignmentTable)
length(unique(assignmentTable$cell)) # only 37,170 cells out of 41,238 are assigned...
table(duplicated(assignmentTable$cell))
assignmentTable |>
  group_by(cell) |>
  mutate(geneSet=replace(geneSet, n_distinct(geneSet)>1, "undefined")) |>
  ungroup() |>
  group_by(geneSet) |> count()


# Get max AUC per cell (doesn't work because GEP5 is higher than everyone else)
# auc_test <- data.frame("max"=apply(getAUC(cells_AUC), 2, max),
#                        "max2nd"= apply(getAUC(cells_AUC), 2, function(x) sort(x, decreasing=T)[2]))
# auc_test$maxdiff <- auc_test$max-auc_test$max2nd
# hist(auc_test$maxdiff, breaks=100)
# 
# auc_test2 <- t(getAUC(cells_AUC)) |>
#   as.data.frame() |>
#   rownames_to_column("cellid") |>
#   as_tibble() |>
#   rename(GEP1=GEP_1, GEP4=GEP_4, GEP5=GEP_5, GEP6=GEP_6, GEP8=GEP_8, GEP11=GEP_11, GEP12=GEP_12) |>
#   pivot_longer(cols=gep_pbmc, names_to="gep", values_to="auc") |>
#   # get max AUC per cell
#   group_by(cellid) |>
#   mutate(score_max=max(auc),
#          gep_max=gep[which.max(auc)],
#          score_2ndmax=sort(auc, decreasing=T)[2],
#          gep_2ndmax=paste0(gep[which(auc==score_2ndmax)], collapse=";")) |>
#   # create column with GEP assignment
#   mutate(gep_assign=ifelse(score_max-score_2ndmax>0.05, gep_max, "undefined")) |>
#   ungroup()
# 
# # Check distribution of difference between max and 2nd max score
# ggplot(auc_test2 %>% select(cellid, score_max, score_2ndmax) %>% distinct())+
#   geom_histogram(aes(x=score_max-score_2ndmax), bins=1000)+
#   geom_vline(xintercept=0.05)+
#   labs(y="# cells", title="AUCell")
# 
# # Quick look on UMAP
# seur.geps_aucell <- seur.geps
# auc_test2.temp <- auc_test2 |>
#   select(cellid, gep_assign) |>
#   distinct()
# tabl(rownames(seur.geps_aucell@meta.data) %in% auc_test2.temp$cellid)
# tabl(rownames(seur.geps_aucell@meta.data[auc_test2.temp$cellid,]) == auc_test2.temp$cellid)
# seur.geps_aucell@meta.data[auc_test2.temp$cellid,"gep_assign"] <- auc_test2.temp$gep_assign
# DimPlot(seur.geps_aucell, group.by = "gep_assign", repel=T, reduction="UMAP_50") +
#   scale_color_manual(values=c(brewer.pal(8, "Dark2"), "#f0f0f0"))


