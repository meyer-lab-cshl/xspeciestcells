# Purpose:  import Lee mouse GD data
# Author: Salomé Carcy
# Date: July 2023


# **************
# 1. IMPORT ####
# **************

# library(RaceID)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)

# Import data
load("/Volumes/CrucialX8/Projects/HumanThymusProject/data/raw_data/mouse_data/GDT_Lee/seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT_no_TCR_gene_in_hvg_no_c7.RData")
seu <- seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT_no_TCR_gene_in_hvg_no_c7
seu <- UpdateSeuratObject(seu) # 2667 cells (more likely to be the right one)

# load("./data/raw_data/mouse_data/GDT_Lee/seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT_no_TCR_gene_in_hvg.RData")
# seu2 <- seurat.qc.norm.rm.doublet.NA.11.14.TCR_RNA_gdT_no_TCR_gene_in_hvg
# seu2 <- UpdateSeuratObject(seu2) # 2681 cells


# Quick visualization
DimPlot(seu, reduction = "umap")
seu@meta.data

# Weird "newcelltype" annotation
table(seu@meta.data$newcelltype, useNA="ifany")
DimPlot(seu, reduction = "umap", group.by="newcelltype")




# ***************
# 2. CLEANUP ####
# ***************

# Give more clear names to clusters
DimPlot(seu, reduction = "umap", group.by="res.0.8", label=T)+ xlim(c(-5, 6)) + ylim(c(-6, 5))

# Try further subclustering of GD1
seu.gdt1 <- subset(seu, subset=res.0.8==1)
DimPlot(seu.gdt1, reduction = "umap", group.by="res.0.8", label=T) + xlim(c(-5, 6)) + ylim(c(-6, 5))
seu.gdt1 <- FindNeighbors(seu.gdt1, reduction="pca", dims=1:25)
seu.gdt1 <- FindClusters(seu.gdt1, resolution=0.8)
DimPlot(seu.gdt1, reduction = "umap", group.by="seurat_clusters", label=T) + xlim(c(-5, 6)) + ylim(c(-6, 5))
# add these subclusters to seu@meta.data
# table(rownames(seu@meta.data[rownames(seu@meta.data) %in% rownames(seu.gdt1@meta.data),]) == rownames(seu.gdt1@meta.data))
seu@meta.data[rownames(seu@meta.data) %in% rownames(seu.gdt1@meta.data),"subclusters_gdt1"] <- seu.gdt1@meta.data$RNA_snn_res.0.8

# Try further subclustering of GD17
seu.gdt17 <- subset(seu, subset=res.0.8==6)
DimPlot(seu.gdt17, reduction = "umap", group.by="res.0.8", label=T) + xlim(c(-5, 6)) + ylim(c(-6, 5))
seu.gdt17 <- FindNeighbors(seu.gdt17, reduction="pca", dims=1:25)
seu.gdt17 <- FindClusters(seu.gdt17, resolution=0.4)
DimPlot(seu.gdt17, reduction = "umap", group.by="seurat_clusters", label=T) + xlim(c(-5, 6)) + ylim(c(-6, 5))
# add these subclusters to seu@meta.data
# table(rownames(seu@meta.data[rownames(seu@meta.data) %in% rownames(seu.gdt17@meta.data),]) == rownames(seu.gdt17@meta.data))
seu@meta.data[rownames(seu@meta.data) %in% rownames(seu.gdt17@meta.data),"subclusters_gdt17"] <- seu.gdt17@meta.data$RNA_snn_res.0.4

seu@meta.data$gd_clusters <- case_when(
  seu@meta.data$res.0.8 == 0 ~ "G4 (immature Tγδ17)",
  seu@meta.data$res.0.8 == 1 & seu@meta.data$subclusters_gdt1 %in% 0:2 ~ "G7-1 (Tγδ1)",
  seu@meta.data$res.0.8 == 1 & seu@meta.data$subclusters_gdt1 == 3 ~ "G7-2 (Tγδ1)",
  seu@meta.data$res.0.8 == 2 ~ "G2 (immature Tγδ1/17)",
  seu@meta.data$res.0.8 == 3 ~ "G1 (Tγδp)",
  seu@meta.data$res.0.8 == 4 ~ "G3 (immature Tγδ1/17)",
  seu@meta.data$res.0.8 == 5 ~ "G5 (immature Tγδ17)",
  seu@meta.data$res.0.8 == 6 & seu@meta.data$subclusters_gdt17 == 0 ~ "G6-1 (Tγδ17)",
  seu@meta.data$res.0.8 == 6 & seu@meta.data$subclusters_gdt17 == 1 ~ "G6-2 (Tγδ17)"
)

DimPlot(seu, reduction = "umap", group.by="gd_clusters", label=T)

# Verify we have raw count data
seu@assays$RNA@counts[6:10,1:5] # perfect!
# seu@assays$RNA@data[6:10,1:5]

# Keep only metadata of interest
# colnames(seu@meta.data)
cols_of_interest <- c("nCount_RNA", "nFeature_RNA", "pct_counts_Mt", "gd_clusters")
seu@meta.data <- seu@meta.data[,cols_of_interest]




# ************
# 3. SAVE ####
# ************

saveRDS(seu, "./data/cross-species/00_Reproduce_UMAPs/ms_gdt_seurobj_lee_subclusters.rds")
