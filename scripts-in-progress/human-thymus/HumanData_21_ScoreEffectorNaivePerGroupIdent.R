# Purpose: plot effectorness & naiveness per group.ident
# Author: Salomé Carcy
# Date: June 2023


# **************
# 1. IMPORT ####
# **************

library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)


# Import data
seur <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")




# **************************
# 2. SCORE EFFECTORNESS ####
# **************************

# Define signature genes
effectorness_genes <- c("HOPX", "GZMB", "GZMK", "ZEB2", "NKG7", "GNLY", "TBX21", "EOMES", "PRF1", "KLRB1", "GZMH", "GZMA", "KLRD1", "CST7", "KLF6",
                        "IFNG", "CCR6", "CCR5", "NCR3", "KLRG1", "IFNGR1")

# naive_genes <- c("SATB1", "TCF7", "LEF1", "CCR7", "SELL", "MYC", "EIF3E", "SOX4", "ID3", "BACH2")
naive_genes <- c("SATB1", "TCF7", "CCR7", "SELL", "SOX4", "ID3", "BACH2", "FOXP1", "TOX2", "STAT5A", "BCL11B", "ETS1", "CD28", "LEF1")


# Run gene scoring
seur <- AddModuleScore(seur, features=list("effectorscore"=effectorness_genes, "naivescore"=naive_genes), name=c("effectorscore", "naivescore"))


# Plot
SCpubr::do_FeaturePlot(seur,
                       split.by="group.ident",
                       order=T,
                       features="effectorscore1",
                       viridis_color_map = "B")
ggsave("./data/human-thymus/HumanData_12_AnalysisByLineage/umap_effectorscore_pergroup2.jpeg", width=20, height=18)

SCpubr::do_FeaturePlot(seur,
                       split.by="group.ident",
                       order=T,
                       features="naivescore2",
                       viridis_color_map = "B")
ggsave("./data/human-thymus/HumanData_12_AnalysisByLineage/umap_naivescore_pergroup.jpeg", width=20, height=18)
