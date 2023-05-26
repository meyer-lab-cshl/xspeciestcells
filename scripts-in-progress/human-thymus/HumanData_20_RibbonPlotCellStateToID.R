# Purpose: make ribbon plot of cell state vs identity
# Author: Salomé Carcy
# Date: May 2023


# **************
# 1. IMPORT ####
# **************

library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
library(ggalluvial)


# Import data
seur <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
cols_integrated <- c("0" = "#f4c40f", "1" = "#b75347", "2" = "#d8443c", "3" = "#e09351", "4" = "#2b9b81", 
                     "5" = "#421401", "6" = "#92c051", "7" = "#9f5691", "8" = "#17154f", "9" = "#74c8c3", 
                     "10" = "#5a97c1", "11" = "gold", "12" = "#a40000", "13" = "#72bcd5", "14" = "grey50",
                     "15" = "orange", "16" = "blueviolet", "17" = "#0a2e57", "18" = "#bdbdbd")
DimPlot(seur, group.by = "new_clusters", label=T, repel=T, reduction="UMAP_50")+
  scale_color_manual(values=cols_integrated)


# Import GEPs
gep_topgenes <- read.csv("./data/human-thymus/HumanData_17_GEPsOnParkData/genes_per_GEP_df_2023-04-07.csv", row.names=1)
dim(gep_topgenes) # 12 GEPs and 100 genes per GEP



# ********************
# 2. COMPUTE GEPs ####
# ********************

gep_list <- list()
for (i in 1:ncol(gep_topgenes)){
  print(paste0("GEP", i))
  gep_allgenes <- gep_topgenes[1:200,i]
  # gep_allgenes <- na.omit(gep_topgenes[,i])
  print(length(gep_allgenes))
  gep_list[[colnames(gep_topgenes)[i]]] <- gep_allgenes
}

seur.geps     <- AddModuleScore(seur, name = "GEP", features=gep_list)


# Sanity check
SCpubr::do_FeaturePlot(seur.geps, reduction="UMAP_50", features=paste0("GEP", 1:length(gep_list)), ncol=6,
                       viridis_color_map = "B", order=T)
# ggsave("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/gapin_cNMF_geps_top200genes.jpeg", width=40, height=20)




# *************************************
# 3. THRESHOLD CELLS BASED ON GEPs ####
# *************************************

## 3.1. Threshold based on GEP with highest score ####
df <- seur.geps@meta.data[,c("GEP1", "GEP4", "GEP5", "GEP6")]
df$score_max <- apply(df, 1, max)
df$gep_max <- colnames(df)[apply(df, 1, which.max)]
head(df)

# First look
# table(rownames(seur.geps@meta.data) == rownames(df))
# seur.geps@meta.data$gep_max <- df$gep_max
# DimPlot(seur.geps, group.by = "gep_max", label=T, repel=T, reduction="UMAP_50")


## 3.2. Threshold based on score value ####

# Look at distribution of counts
hist(df$GEP1, breaks=1000) # 0.1
hist(df$GEP4, breaks=1000) # 0.15
hist(df$GEP5, breaks=1000) # 0.1?...
hist(df$GEP6, breaks=1000) # 0.1?...

# Remove any score below 0.1
table(df$score_max < 0.05)
df[df$gep_max == "GEP1" & df$score_max < 0.05, "score_max"] <- 0.0
df[df$gep_max == "GEP4" & df$score_max < 0.05, "score_max"] <- 0.0
df[df$gep_max == "GEP5" & df$score_max < 0.05, "score_max"] <- 0.0
df[df$gep_max == "GEP6" & df$score_max < 0.05, "score_max"] <- 0.0
table(df$score_max == 0)
df[df$score_max == 0, "gep_max"]   <- "other"

# Final look
# table(rownames(seur.geps@meta.data) == rownames(df)) # sanity check
seur.geps@meta.data$gep_max <- df$gep_max
DimPlot(seur.geps, group.by = "gep_max", repel=T, reduction="UMAP_50") + scale_color_manual(values=c(brewer.pal(4, "Dark2"), "#f0f0f0"))
VlnPlot(seur.geps, group.by = "gep_max", features=c("GEP1", "GEP4", "GEP5", "GEP6"))
# VlnPlot(seur.geps, group.by = "new_clusters", split.by="gep_max", features=c("GEP1", "GEP4", "GEP5", "GEP6"))

# How many PBMCs are unassigned?
table(seur.geps@meta.data[seur.geps@meta.data$Tissue=="PBMC", "gep_max"]) # 4,621 unassigned vs ...29,661 assigned

# Take a look at unassigned cells on UMAP
DimPlot(subset(seur.geps, subset= Tissue=="PBMC" & gep_max=="other"), group.by = "new_clusters", repel=T, reduction="UMAP_50") +
  scale_color_manual(values=cols_integrated)
FeaturePlot(subset(seur.geps, subset= Tissue=="PBMC" & gep_max=="other"),
            features="GEP1",
            reduction="UMAP_50")

# Look in 2D scatter plot
head(df)
df$tissue <- seur.geps@meta.data$Tissue
df <- df[df$tissue=="PBMC",]
library(GGally)
ggpairs(df, columns=1:4, aes(color=gep_max, alpha=0.9), legend=c(1,1),
        upper=list(continuous=wrap("points", size=0.05)),
        lower=list(continuous=wrap("points", size=0.05)))+
  scale_color_manual(values=c(brewer.pal(4, "Dark2"), "#f0f0f0")) +
  scale_fill_manual(values=c(brewer.pal(4, "Dark2"), "#f0f0f0")) +
  theme_cowplot()+
  theme(legend.position="right")
# ggsave("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/ggpairs_pbmc_allgeps_colorbygepcategory.jpeg", width=10, height=10)

VlnPlot(subset(seur.geps, subset= Tissue=="PBMC"), group.by = "gep_max", features=paste0("GEP", 1:12), same.y.lims=T, ncol=4)
# ggsave("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/vlnplot_pbmc_allgeps_gepcat.jpeg", width=20, height=18)




# *******************
# 4. RIBBON PLOT ####
# *******************

# Format data
counts <- seur.geps@meta.data %>%
  as_tibble() %>%
  filter(Tissue=="PBMC") %>%
  # get nb of cells per gep assignment
  group_by(cell.ident, gep_max) %>%
  summarise(ncells=n()) %>%
  # get %cells in each gep assignment
  ungroup() %>%
  group_by(cell.ident) %>%
  mutate(totalcells=sum(ncells),
         freq = ncells*100/totalcells) %>%
  ungroup() %>%
  # rename
  mutate(gep_max = ifelse(gep_max=="GEP1", "Th17?\n(GEP1)",
                          ifelse(gep_max=="GEP4", "Temra\n(GEP4)",
                                 ifelse(gep_max=="GEP5", "Tnaive\n(GEP5)",
                                        ifelse(gep_max=="GEP6", "Tcm\n(GEP6)", "other")))),
         cell.ident=replace(cell.ident, cell.ident=="NKT", "iNKT"))
  

# Plot
ggplot(data=counts, aes(axis1=cell.ident, axis2=gep_max, y=freq)) +
  geom_alluvium(aes(fill=cell.ident))+
  geom_stratum()+
  geom_text(stat="stratum", aes(label=after_stat(stratum)))+
  scale_fill_manual(values=c("#2ca25f", "#dd1c77", "#045a8d", "#8856a7", "#bdc9e1"), name="cell type")+
  # scale_x_discrete(limits=c("Cell type", "GEP")) + theme_classic()
  theme_void()
# ggsave("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/ribbon_cellstateproportions.jpeg", width=6, height=6)

# Sanity check
# VlnPlot(subset(seur.geps, subset= Tissue=="PBMC"), group.by = "cell.ident", features=c("GEP1", "GEP4", "GEP5", "GEP6"),
#         same.y.lims=T, ncol=4)+
#   scale_fill_manual(values=c("#2ca25f", "#dd1c77", "#045a8d", "#8856a7", "#bdc9e1"))


