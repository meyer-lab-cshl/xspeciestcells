# Purpose: plot CD1d and MR1 expression in Park dataset
# Author: Salomé Carcy
# Date: August 2023


# **************
# 1. IMPORT ####
# **************

library(Seurat)
library(ggplot2)
library(tidyverse)
# library(dplyr)
library(cowplot)
library(RColorBrewer)
library(scico)
source("./scripts-final/colors_universal.R")


# Import data
seur.park <- readRDS("~/Projects/HumanThymusProject/data/human-thymus/HumanData_16_ParkIntegration/park_full_seu_gene_names.rds")
seur_metadata <- read.csv("~/Projects/HumanThymusProject/data/human-thymus/HumanData_16_ParkIntegration/park_fullmetadata.csv")


# Incorporate metadata in seurat object
table(seur_metadata$index == rownames(seur.park@meta.data), useNA="ifany")
seur.park@meta.data <- cbind(seur.park@meta.data, seur_metadata)



# ******************
# 2. PLOT UMAPs ####
# ******************

# Clusters
# colnames(seur.park@meta.data)
table(seur.park@meta.data$Anno_level_3, useNA="ifany")
seur.park@meta.data[seur.park@meta.data$Anno_level_3=="TEC(neuro)", "Anno_level_3"] <- "mTEC"
seur.park@meta.data[seur.park@meta.data$Anno_level_3=="TEC(myo)", "Anno_level_3"] <- "mTEC"
table(seur.park@meta.data$Anno_level_3=="mTEC", useNA="ifany") # should have 6928 mTEC


Idents(seur.park) <- "Anno_level_3"
p1 <- SCpubr::do_DimPlot(seur.park,
                   reduction="UMAP",
                   # group.by="Anno_level_fig1",
                   # group.by="Anno_level_3",
                   # cells.highlight = rownames(seur.park@meta.data[seur.park@meta.data$Anno_level_3 %in% c("cTEC", "DN", "DP"),]),
                   idents.keep=c("cTEC", "mTEC", "DN", "DP", "CD8αα(I)", "Mono", "T_naive", "γδT", "B_naive", "B_memory"),
                   na.value="grey90",
                   legend.position="right", legend.ncol=1, font.size=30, legend.icon.size=10)+
  scale_color_manual(values=grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(11))
  # scale_color_manual(values=RColorBrewer::brewer.pal(8, "Accent"))
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/umap_clusters_level3_highlight.jpeg", width=10, height=8)


# Plot CD1d expression
p2 <- SCpubr::do_FeaturePlot(seur.park, 
                       features = "CD1D", order = T,
                       # plot.title = "CD1D",
                       border.color = "black", border.size = 2,
                       reduction = "UMAP", pt.size = 1.2, legend.position = "right", font.size=30) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0.1, end = 0.85, direction = -1)
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/umap_cd1d_orderT.jpeg", width=10, height=8)


p1 | p2
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/umap_combined.jpeg", width=20, height=8)


VlnPlot(seur.park, features = "CD1D", group.by = "Anno_level_3", raster=F)+
  scale_fill_manual(values=grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(44))+
  labs(y="CD1D (normalized expression)", title="", x="")+
  theme(legend.position="none",
        axis.text.x=element_text(size=15),
        axis.title.y=element_text(size=15))
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/vlnplt_cd1d_level3.jpeg", width=12, height=6)
