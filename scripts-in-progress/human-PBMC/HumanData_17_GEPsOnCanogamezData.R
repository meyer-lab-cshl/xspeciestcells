###
# Purpose: Score GEPs on Cano-Gamez data
# Date: May 3rd 2023
# Author: Salomé Carcy
###


#### IMPORT ####

# Libraries
library(ggplot2)
library(Seurat)
library(cowplot)
library(tidyverse)

# Data
seur.cano <- readRDS("./data/raw_data/human_data/canogamez-pbmcT/ALL_cells_no_erythro_v3.rds")
print(seur.cano) # 43,112 cells
DimPlot(seur.cano, reduction="umap", group.by="cytokine.condition")


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




# ************************
# 2. RUN MODULE SCORE ####
# ************************

# Make list of GEPs
gep_list <- list()
for (i in 1:ncol(gep_topgenes)){
  print(paste0("GEP", i))
  gep_allgenes <- gep_topgenes[1:200,i]
  gep_genesincano <- gep_allgenes[gep_allgenes %in% rownames(seur.cano)]
  print(length(gep_genesincano))
  gep_list[[colnames(gep_topgenes)[i]]] <- gep_genesincano
}
# names(gep_list) # sanity check

# Check for each GEP how many genes are not present in canogamez seurat object
# gep_list_genes_notincano <- lapply(gep_list, function(x) x[!x %in% rownames(seur.cano)])
# lengths(gep_list_genes_notincano)

# Run module score on canogamez data and our data
seur.cano.gep <- AddModuleScore(seur.cano, name = "GEP", features=gep_list)
seur.geps     <- AddModuleScore(seur, name = "GEP", features=gep_list)




# *******************************
# 3. VISUALIZE GEPs ON UMAPs ####
# *******************************

# Cano-Gamez
SCpubr::do_FeaturePlot(seur.cano.gep, reduction="umap", features=paste0("GEP", 1:length(gep_list)), ncol=6,
                       viridis_color_map = "B")
# ggsave("./data/human-thymus/HumanData_17_GEPsOnCanogamezData/cano_cNMF_geps.jpeg", width=40, height=20)
DimPlot(seur.cano, reduction="umap", group.by="cytokine.condition", label=T) |
  DimPlot(seur.cano, reduction="umap", group.by="cluster_id", label=T) + scale_color_manual(values=colorRampPalette(brewer.pal(8, "Paired"))(22))
# ggsave("./data/human-thymus/HumanData_17_GEPsOnCanogamezData/cano_umap.jpeg", width=15, height=8)


# Our data
SCpubr::do_FeaturePlot(seur.geps, reduction="UMAP_50", features=paste0("GEP", 1:length(gep_list)), ncol=6,
                       viridis_color_map = "B")
ggsave("./data/human-thymus/HumanData_17_GEPsOnCanogamezData/gapin_cNMF_geps_canogenes.jpeg", width=40, height=20)




# ***************************************
# 3. VISUALIZE GEPs ON CELL CLUSTERS ####
# ***************************************

VlnPlot(seur.cano.gep, group.by="cluster_id", features=c("GEP1", "GEP4", "GEP5", "GEP6"), ncol=2)
# ggsave("./data/human-thymus/HumanData_17_GEPsOnCanogamezData/cano_cNMF_geps_vlnplot.jpeg", width=20, height=15)


