# Purpose: make ribbon plot of cell state vs identity
# Author: Salomé Carcy
# Date: August 2023


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

# gep_list <- list()
# for (i in 1:ncol(gep_topgenes)){
#   print(paste0("GEP", i))
#   # gep_allgenes <- gep_topgenes[1:200,i][!is.na(gep_topgenes[1:200,i])]
#   gep_allgenes <- gep_topgenes[,i][!is.na(gep_topgenes[,i])]
#   print(length(gep_allgenes))
#   gep_list[[colnames(gep_topgenes)[i]]] <- gep_allgenes
# }
# 
# seur.geps     <- AddModuleScore(seur, name = "GEP", features=gep_list, seed=123)
# # seur.geps <- readRDS("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/seuratobj_gepscores_allgenes.rds")
# 
# # Sanity check
# gep_pbmc <- c("GEP1", "GEP4", "GEP5", "GEP6", "GEP8", "GEP12")
# SCpubr::do_FeaturePlot(subset(seur.geps, Tissue=="PBMC"), reduction="UMAP_50", features=gep_pbmc, ncol=3,
#                        viridis_color_map = "B", order=T)
# # ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/gapin_cNMF_geps_allgenes_pbmc_orderTRUE.jpeg", width=17, height=15)
# 



# ******************************
# 3. GET GEP USAGE PER CELL ####
# ******************************

# Import GEP usage
gep_usage <- read.table("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/cNMF_output/imputed_cNMF.usages.k_12.dt_0_02.consensus.txt", header=T)
dim(gep_usage)
nrow(gep_usage)==ncol(seur) # rows are cells
head(gep_usage)

# Check that every row sums to 1
hist(rowSums(gep_usage), breaks=100)
min(rowSums(gep_usage))
max(rowSums(gep_usage))

# Add to seurat object
# seur.geps <- seur
colnames(gep_usage) <- paste0("gep", 1:12, "_usage")
table(rownames(gep_usage)==rownames(seur.geps@meta.data))
seur.geps@meta.data <- cbind(seur.geps@meta.data, gep_usage)

# Plot
gepusage_pbmc <- c("gep1_usage", "gep4_usage", "gep5_usage", "gep6_usage", "gep8_usage", "gep11_usage", "gep12_usage")
SCpubr::do_FeaturePlot(subset(seur.geps, Tissue=="PBMC"), reduction="UMAP_50", features=gepusage_pbmc,
                       ncol=3, viridis_color_map = "B", order=F)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/gapin_cNMF_geps_usage_pbmc_orderFALSE.jpeg", width=17, height=15)




# ****************************
# 4. ASSIGN CELLS TO GEPs ####
# ****************************

# Keep only PBMCs and remove gep7 (because it's batch driven)
df <- gep_usage[grep("PBMC", rownames(gep_usage), value=T), colnames(gep_usage)[!colnames(gep_usage) %in% "gep7_usage"] ]
dim(df) # 41,238 cells

# Get the max
df$score_max <- apply(df, 1, max)
df$gep_assign <- gsub("_.*", "", colnames(df)[apply(df, 1, which.max)])
head(df)
tabl(df$score_max<0.5) # only ~8,000 cells have a max score below 0.5
tabl(df$gep_assign)

# Look at distribution of max score
# ggplot(df)+
#   geom_density(aes(x=score_max, color=gep_assign))

# How many have a very little difference between the max score and 2nd max score?
df$score_2ndmax <- apply(df[,1:7], 1, function(x) sort(x, decreasing = T)[2])
tabl(df$score_max-df$score_2ndmax<0.1) # 7% of the cells have less than 0.1 difference in score between top score and 2nd top score
ggplot(df)+
  geom_histogram(aes(x=score_max-score_2ndmax), bins=1000)+
  labs(y="# cells", title="assign cells to GEPs based on cNMF usage")+
  xlim(0,1) + ylim(0,300)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/method1_rawscore_diff_btw_max_2ndmax.jpeg", width=6, height=4)

# First look
table(rownames(seur.geps@meta.data[seur.geps$Tissue=="PBMC",]) == rownames(df))
seur.geps@meta.data[seur.geps$Tissue=="PBMC","gep_assign"] <- df$gep_assign

# Final look
DimPlot(seur.geps, group.by = "gep_assign", repel=T, reduction="UMAP_50") + scale_color_manual(values=c(brewer.pal(7, "Dark2"), "#f0f0f0"))
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/method5_gepusage_umap.jpeg", width=6, height=6)
# VlnPlot(seur.geps[,seur.geps$Tissue=="PBMC"], group.by = "gep_assign", features=gep_pbmc, same.y.lims=T)
# ggsave("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/method_comparison/vlnpreview_gepscomparisononrawscores.jpeg", width=10, height=10)




# *******************
# 5. RIBBON PLOT ####
# *******************

# Format data
counts <- seur.geps@meta.data %>%
  as_tibble() %>%
  filter(Tissue=="PBMC") %>%
  # get nb of cells per gep assignment
  group_by(cell.ident, gep_assign) %>%
  summarise(ncells=n()) %>%
  # get %cells in each gep assignment
  ungroup() %>%
  group_by(cell.ident) %>%
  mutate(totalcells=sum(ncells),
         freq = ncells*100/totalcells) %>%
  ungroup() %>%
  # rename a few variables
  mutate(gep_assign=toupper(gep_assign),
         cell.ident=replace(cell.ident, cell.ident=="NKT", "iNKT"))


# Main figure
counts %>%
  mutate(gep_assign=replace(gep_assign, !gep_assign%in%c("GEP1", "GEP4", "GEP5", "GEP6"), "other"),
         gep_assign = factor(gep_assign, levels=c("GEP1", "GEP4", "GEP5", "GEP6", "other"))) %>%
  filter(gep_assign != "other") %>%
ggplot(aes(axis1=cell.ident, axis2=gep_assign, y=freq)) +
  geom_alluvium(aes(fill=cell.ident))+
  geom_stratum()+
  geom_text(stat="stratum", aes(label=after_stat(stratum)), size=8)+
  scale_fill_manual(values=c("#2ca25f", "#dd1c77", "#045a8d", "#8856a7", "#bdc9e1"), name="cell type")+
  # scale_x_discrete(limits=c("Cell type", "GEP")) + theme_classic()
  theme_void()+
  theme(legend.position="none")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/method5_ribbon.pdf", width=6, height=6)



# Supp figure
counts %>%
  mutate(gep_assign = factor(gep_assign, levels=c("GEP1", "GEP3", "GEP4", "GEP5", "GEP6", "GEP8", "GEP12"))) %>%
ggplot(aes(axis1=cell.ident, axis2=gep_assign, y=freq)) +
  geom_alluvium(aes(fill=cell.ident))+
  geom_stratum()+
  geom_text(stat="stratum", aes(label=after_stat(stratum)))+
  # ggrepel::geom_text_repel(
  #   aes(label = gep_assign),
  #   stat = "stratum", size = 4, direction = "y", nudge_x = .3
  # ) +
  scale_fill_manual(values=c("#2ca25f", "#dd1c77", "#045a8d", "#8856a7", "#bdc9e1"), name="cell type")+
  # scale_x_discrete(limits=c("Cell type", "GEP")) + theme_classic()
  theme_void()+
  labs(title="")

