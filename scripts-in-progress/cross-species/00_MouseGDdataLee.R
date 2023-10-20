# Purpose:  import Lee mouse GD data
# Author: Salomé Carcy
# Date: July 2023


# **************
# 1. IMPORT ####
# **************

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

seu@meta.data$gd_subclusters <- case_when(
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

DimPlot(seu, reduction = "umap", group.by="gd_subclusters", label=T)


# Make a column with higher level annotation
seu@meta.data$gd_clusters <- case_when(
  seu@meta.data$gd_subclusters == "G1 (Tγδp)"             ~ "Tγδp",
  seu@meta.data$gd_subclusters == "G2 (immature Tγδ1/17)" ~ "immature Tγδ1/17",
  seu@meta.data$gd_subclusters == "G3 (immature Tγδ1/17)" ~ "immature Tγδ1/17",
  seu@meta.data$gd_subclusters == "G4 (immature Tγδ17)"   ~ "immature Tγδ17",
  seu@meta.data$gd_subclusters == "G5 (immature Tγδ17)"   ~ "immature Tγδ17",
  seu@meta.data$gd_subclusters == "G6-1 (Tγδ17)"          ~ "Tγδ17",
  seu@meta.data$gd_subclusters == "G6-2 (Tγδ17)"          ~ "Tγδ17",
  seu@meta.data$gd_subclusters == "G7-1 (Tγδ1)"           ~ "Tγδ1",
  seu@meta.data$gd_subclusters == "G7-2 (Tγδ1)"           ~ "Tγδ1"
)
DimPlot(seu, reduction = "umap", group.by="gd_clusters", label=T)

# sanity check
# table(seu@meta.data[,c("gd_subclusters", "gd_clusters")], useNA="ifany")

# Verify we have raw count data
seu@assays$RNA@counts[6:10,1:5] # perfect!
# seu@assays$RNA@data[6:10,1:5]

# Keep only metadata of interest
# colnames(seu@meta.data)
cols_of_interest <- c("celltype", "newcelltype", "VA14JA18", "VB8", "VB7", "VA19JA33", "VB6", "VG123456", "ident", "res.0.8",
                      "nCount_RNA", "nFeature_RNA", "pct_counts_Mt", "gd_clusters", "gd_subclusters")
seu@meta.data <- seu@meta.data[,cols_of_interest]




# ************
# 3. SAVE ####
# ************

saveRDS(seu, "./data/cross-species/00_Reproduce_UMAPs/ms_gdt_seurobj_lee_subclusters.rds")




# *******************************************
# 4. CONVERT ENSEMBL IDs TO GENE SYMBOLS ####
# *******************************************

library(SingleCellExperiment)
library(biomaRt)
sce <- as.SingleCellExperiment(seu)

# Get Mouse ENSEMBL - gene symbol correspondance from biomart
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

genemap <- getBM(attributes=c('ensembl_gene_id', "external_gene_name"),
                 filters = 'ensembl_gene_id',
                 values = rownames(sce),
                 mart = ensembl) %>%
  dplyr::rename(gene_id = ensembl_gene_id,
                symbol = external_gene_name) %>%
  distinct()

# how many gene names available?...
table(rownames(sce) %in% genemap$gene_id, useNA="ifany") # 25,645

# Change EMSEMBL IDs to gene symbols in SCE object
featureData <- tibble(gene_id=rownames(sce)) %>%
  left_join(genemap, by="gene_id") %>%
  drop_na() %>%
  filter(duplicated(symbol) == FALSE) %>%
  filter(symbol != "")

sce_final <- subset(sce, rownames(sce) %in% gene_id$gene_id)

mcols(sce_final) <- DataFrame(mcols(sce_final), featureData)
rownames(sce_final) <- rowData(sce_final)$symbol

# convert back to seurat
seu_final <- as.Seurat(sce_final)

# add the dimension reductions that got lost
seu_final@reductions$pca <- seu@reductions$pca
seu_final@reductions$tsne <- seu@reductions$tsne
seu_final@reductions$umap <- seu@reductions$umap
seu_final@meta.data <- seu_final@meta.data[,2:16]

DimPlot(seu_final, reduction="umap", group.by="gd_subclusters")
rownames(seu_final)[1:5]

# save for Laurent
saveRDS(seu_final, "./data/cross-species/00_Reproduce_UMAPs/ms_gdt_seurobj_lee_genesymbols.rds")




