###
# Purpose: Import & plot human iNKT data (curtesy of Laurent Gapin) with SCTransform and Harmony
# Date: Oct 3rd 2022
# Author: Salomé Carcy
###


#### IMPORT ####

# Libraries
library(ggplot2)
library(cowplot)
library(Seurat)
library(SeuratWrappers)
library(harmony)
library(Matrix)
library(tidyverse)

# Data (sorted human iNKT cells from 3 subjects)
path.plots <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/00_Reproduce_UMAPs"
path.data <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/human_data"
human5 <- read.csv(file.path(path.data, "CUThy13_220225_SampleTag05_hs_NKT_RSEC_MolsPerCell.csv"), sep=",", header=T, skip=7) # 404 cells
human8 <- read.csv(file.path(path.data, "CUTHY11BDRscRNA_seq_091621_SampleTag08_hs_NKT_RSEC_MolsPerCell.csv"), sep=",", header=T, skip=8) # 1913 cells
human12 <- read.csv(file.path(path.data, "CUTHY12BDRscRNA_seq_211101_SampleTag12_hs_NKT_RSEC_MolsPerCell.csv"), sep=",", header=T, skip=8) # 344 cells

# Any cluster ID in the colnames?
# table(stringr::str_detect(colnames(human5),"[[:lower:]]"))
# colnames(human5)[stringr::str_detect(colnames(human5),"[[:lower:]]")]


# Current df format:
# Rownames || Cell_Index | Gene1 | Gene2 | ...
#   -      ||   421953   |   0   |   0   | ...
#   -      ||   459121   |   0   |   3   | ...

# Invert the dataframes (cells as columns and genes as rows) and transform to dgCMatrix
convert_to_matrix <- function(df){
  rownames(df) <- df$Cell_Index
  df$Cell_Index <- NULL
  df <- t(df)
  df <- Matrix(df, sparse=T)
  return(df)
}

mat.hu5  <- convert_to_matrix(human5) # 404 cells
mat.hu8  <- convert_to_matrix(human8) # 1913 cells
mat.hu12 <- convert_to_matrix(human12) # 344 cells




#### PRE-PROCESSING ####

# Create Seurat Object
seur.h5  <- CreateSeuratObject(mat.hu5, project="Hu_Thymus_NKT_5")
seur.h8  <- CreateSeuratObject(mat.hu8, project="Hu_Thymus_NKT_8")
seur.h12 <- CreateSeuratObject(mat.hu12, project="Hu_Thymus_NKT_12")

# Combine into one seurat object
seur.combined <- merge(seur.h5, y=c(seur.h8,seur.h12), add.cell.ids=c('Hu5', 'Hu8', 'Hu12'), project="HuNKT") # 2661 cells

# QC
seur.combined[["percent.mt"]] <- PercentageFeatureSet(seur.combined, pattern = "^MT\\.")
head(seur.combined@meta.data, 5)
FeatureScatter(seur.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seur.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(seur.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=.01)
# ggsave(file.path(path.plots, "hu_qc.jpeg"), width=10, height=6)
seur.combined <- subset(seur.combined, subset = nFeature_RNA >= 500 & nFeature_RNA < 4000 & nCount_RNA >=500 & percent.mt < 25)
table(seur.combined@meta.data$orig.ident) # 399 cells hu5; 1828 cells hu8; 331 cells hu12

# Normalize with SCTransform and find variable features
# based on https://github.com/immunogenomics/harmony/issues/41 and https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/harmony.html
seur.list <- SplitObject(seur.combined, split.by = "orig.ident")
seur.list <- lapply(X = seur.list, FUN = SCTransform) # weird, SCTransform adds features?!
nkt.features <- SelectIntegrationFeatures(object.list = seur.list, nfeatures = 2000)

# Merge back seurat objects to run Harmony afterwards
huNKT_merged <- merge(seur.list$Hu_Thymus_NKT_5, y = c(seur.list$Hu_Thymus_NKT_8, seur.list$Hu_Thymus_NKT_12),
                      project = "HuNKT", merge.data = TRUE) # weird, now there are 40555 genes
VariableFeatures(huNKT_merged) <- nkt.features

# Run PCA and Harmony for integration (don't scale data when using SCTransform)
huNKT_merged <- RunPCA(object = huNKT_merged, assay = "SCT", npcs = 50)
# DimPlot(huNKT_merged, reduction = "pca", group.by = "orig.ident")
# ElbowPlot(huNKT_merged)
huNKT_merged <- RunHarmony(object = huNKT_merged,
                              assay.use = "SCT",
                              reduction = "pca",
                              dims.use = 1:50,
                              group.by.vars = "orig.ident",
                              plot_convergence = F) # converged only after 2 iterations
huNKT_merged <- RunUMAP(object = huNKT_merged, assay = "SCT", reduction = "harmony", dims = 1:30)
huNKT_merged <- FindNeighbors(object = huNKT_merged, assay = "SCT", reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.8)
# Display UMAP
p1 <- DimPlot(huNKT_merged, group.by = c("orig.ident"), pt.size = 0.1, label = F) + theme(legend.position="top", title = element_text(size=0))
p2 <- DimPlot(huNKT_merged, pt.size = 0.1, label = TRUE) + labs(title="Integrated (2558 cells)") + theme(legend.position="none")
p3 <- DimPlot(huNKT_merged, pt.size = 0.1, split.by = "orig.ident", ncol = 3)
ggdraw()+
  draw_plot(p1, x=0, y=0, width=.25, height=1)+
  draw_plot(p2, x=.25, y=0, width=.25, height=1)+
  draw_plot(p3, x=.5, y=0, width=.5, height=1)




