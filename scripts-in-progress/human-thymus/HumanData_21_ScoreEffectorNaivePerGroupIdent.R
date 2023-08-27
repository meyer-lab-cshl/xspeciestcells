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
                       min.cutoff=0,
                       max.cutoff=1,
                       features="effectorscore1",
                       viridis_color_map = "B")
# ggsave("./data/human-thymus/HumanData_12_AnalysisByLineage/umap_effectorscore_pergroup2.jpeg", width=20, height=18)
ggsave("./data/human-thymus/HumanData_24_PlotsForFOCISconf/umap_effectorscore_pergroup2.jpeg", width=20, height=18)

SCpubr::do_FeaturePlot(seur,
                       split.by="group.ident",
                       order=T,
                       min.cutoff=0,
                       features="naivescore2",
                       viridis_color_map = "B")
# ggsave("./data/human-thymus/HumanData_12_AnalysisByLineage/umap_naivescore_pergroup.jpeg", width=20, height=18)
ggsave("./data/human-thymus/HumanData_24_PlotsForFOCISconf/umap_naivescore_pergroup.jpeg", width=20, height=18)




# ******************************
# 3. PLOT GENES OF INTEREST ####
# ******************************

SCpubr::do_FeaturePlot(seur,
                       split.by="group.ident",
                       order=T,
                       features="TBX21",
                       ncol=5,
                       viridis_color_map = "B")
# ggsave("~/Projects/HumanThymusProject/data/human-thymus/HumanData_24_PlotsForFOCISconf/umap_tbx21.jpeg", width=15, height=8)


SCpubr::do_DimPlot(seur,
                   cells.highlight = rownames(seur@meta.data[seur@meta.data$TCR_Beta_Delta_V_gene_Dominant == "TRDV1*01" &
                                                               !is.na(seur@meta.data$TCR_Beta_Delta_V_gene_Dominant),]))



# GEP USAGE GD

# Import GEP usage
gep_usage <- read.table("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/cNMF_output/imputed_cNMF.usages.k_12.dt_0_02.consensus.txt", header=T)
dim(gep_usage)
nrow(gep_usage)==ncol(seur) # rows are cells
colnames(gep_usage) <- paste0("gep", 1:12, "_usage")
table(rownames(gep_usage)==rownames(seur@meta.data), useNA="ifany")
seur@meta.data <- cbind(seur@meta.data, gep_usage)

# VD2 VG9
# SCpubr::do_DimPlot(seur,
#                    cells.highlight = rownames(seur@meta.data[seur@meta.data$TCR_Beta_Delta_V_gene_Dominant %in% c("TRDV2*01", "TRDV2*02", "TRDV2*03") &
#                                                                !is.na(seur@meta.data$TCR_Beta_Delta_V_gene_Dominant),]))
SCpubr::do_FeaturePlot(seur,
                   cells.highlight = rownames(seur@meta.data[seur@meta.data$TCR_Beta_Delta_V_gene_Dominant %in% c("TRDV2*01", "TRDV2*02", "TRDV2*03") &
                                                               !is.na(seur@meta.data$TCR_Beta_Delta_V_gene_Dominant),]),
                   features="gep4_usage", order=T, viridis_color_map = "B", plot.title="TRDV2", min.cutoff = 0, max.cutoff = 0.8)
# ggsave("~/Projects/HumanThymusProject/data/human-thymus/HumanData_24_PlotsForFOCISconf/umap_trdv2_gep4usage.jpeg", width=6, height=7)


SCpubr::do_FeaturePlot(seur,
                   cells.highlight = rownames(seur@meta.data[seur@meta.data$TCR_Alpha_Gamma_V_gene_Dominant %in% c("TRGV9*01", "TRGV9*02") &
                                                               !is.na(seur@meta.data$TCR_Alpha_Gamma_V_gene_Dominant),]),
                   features="gep4_usage", order=T, viridis_color_map = "B", plot.title="TRGV9", min.cutoff = 0, max.cutoff = 0.8)
# ggsave("~/Projects/HumanThymusProject/data/human-thymus/HumanData_24_PlotsForFOCISconf/umap_trgv9_gep4usage.jpeg", width=6, height=7)

SCpubr::do_FeaturePlot(seur,
                       cells.highlight = rownames(seur@meta.data[seur@meta.data$TCR_Beta_Delta_V_gene_Dominant %in% c("TRDV2*01", "TRDV2*02", "TRDV2*03") &
                                                                   !is.na(seur@meta.data$TCR_Beta_Delta_V_gene_Dominant) &
                                                                 seur@meta.data$TCR_Alpha_Gamma_V_gene_Dominant %in% c("TRGV9*01", "TRGV9*02") &
                                                                   !is.na(seur@meta.data$TCR_Alpha_Gamma_V_gene_Dominant),]),
                       features="gep4_usage", order=T, viridis_color_map = "B", plot.title="VD2-VG9", min.cutoff = 0, max.cutoff = 0.8)
# ggsave("~/Projects/HumanThymusProject/data/human-thymus/HumanData_24_PlotsForFOCISconf/umap_trdv2_trvg9_gep4usage.jpeg", width=6, height=7)


# VD1
SCpubr::do_FeaturePlot(seur,
                       cells.highlight = rownames(seur@meta.data[seur@meta.data$TCR_Beta_Delta_V_gene_Dominant == "TRDV1*01" &
                                                                   !is.na(seur@meta.data$TCR_Beta_Delta_V_gene_Dominant),]),
                       features="gep1_usage", order=T, viridis_color_map = "B", plot.title="TRDV1", min.cutoff = 0, max.cutoff = 0.8)
# ggsave("~/Projects/HumanThymusProject/data/human-thymus/HumanData_24_PlotsForFOCISconf/umap_trdv1_gep1usage.jpeg", width=6, height=7)




