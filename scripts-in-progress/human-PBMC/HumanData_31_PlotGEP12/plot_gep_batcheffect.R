# Purpose:
# Author: Salomé Carcy
# Date:




# **************
# 1. IMPORT ####
# **************

# Import librairies
library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)
library(Seurat)
source("./scripts-final/colors_universal.R")

# Import data
seur <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_08_28_23.RDS")
DimPlot(seur, group.by = "new_clusters", label=T, repel=T, reduction="UMAP_50")+
  scale_color_manual(values=cols_integrated)

gep_usage <- read.table("./data/human-PBMC/HumanData_20_RibbonPlotCellStateToID/cNMF_output/non_imputed_cNMF_allcells.usages.k_12.dt_0_02.consensus.txt", header=T)
dim(gep_usage)
colnames(gep_usage) <- paste0("gep", c(2,5,3,1,4,12,6,7,8,10,9,11), "_usage")
head(gep_usage)
table(rownames(gep_usage)==rownames(seur@meta.data))
seur@meta.data <- cbind(seur@meta.data, gep_usage)



# ****************
# 2. ANALYSIS ####
# ****************

#________________________________
## 3.1. Plot GEP12 per batch ####

ggrastr::rasterise(SCpubr::do_FeaturePlot(seur,
                       features="gep12_usage",
                       split.by="Batch",
                       ncol=3,
                       order=F,
                       use_viridis = T,
                       viridis.palette = "B",
                       legend.title="GEP12 (cNMF usage)")+
                     scale_color_viridis_c(option="B"),
                   layers="Point", dpi=300)
ggsave("./scripts-in-progress/human-PBMC/HumanData_31_PlotGEP12/plots/gep12_usage_orderF.jpeg", width=10, height=11)

## /end ####

