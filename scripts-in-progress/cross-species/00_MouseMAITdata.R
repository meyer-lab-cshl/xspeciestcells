###
# Purpose: Import & plot mouse MAIT data (Legoux paper)
# Date: Aug 9th 2022
# Author: Salomé Carcy
###


#### IMPORT ####

# Libraries
library(ggplot2)
library(Seurat)
# remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
library(cowplot)
library(tidyverse)

# Data
# path.plots <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/00_Reproduce_UMAPs"
path.data <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/mouse_data/MAIT_Legoux/Processed_data"
mouse1 <- read.table(file.path(path.data, "processed_ikura_tobiko_wt01.txt"), header=T)
mouse2 <- read.table(file.path(path.data, "processed_ikura_tobiko_wt02.txt"), header=T)

mouse1[1:5,1:5]
mouse2[1:5,1:5]

# Put gene names as rownames
rownames(mouse1) <- mouse1$symbol
rownames(mouse2) <- mouse2$symbol
mouse1$symbol <- NULL
mouse2$symbol <- NULL
# sanity check
mouse1[1:5,1:5]
mouse2[1:5,1:5]




#### PRE-PROCESSING ####

# Create Seurat Object
seur.ms1 <- CreateSeuratObject(counts=mouse1, project="B6_Thymus_MAIT_1", min.cells=3)
seur.ms2 <- CreateSeuratObject(counts=mouse2, project="B6_Thymus_MAIT_2", min.cells=3)

# Combine into one seurat object
seur.combined <- merge(seur.ms1, y = seur.ms2, add.cell.ids = c('B61', 'B62'), project = "MAIT")


# QC
seur.combined[["percent.mt"]] <- PercentageFeatureSet(seur.combined, pattern = "^mt-")
head(seur.combined@meta.data, 5)
FeatureScatter(seur.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seur.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(seur.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=.01)
# ggsave(file.path(path.plots, "ms_qc.jpeg"), width=10, height=6)
seur.combined <- subset(seur.combined, subset = nFeature_RNA >= 500 & percent.mt < 5)
table(seur.combined@meta.data$orig.ident) # 3,824 cells mouse 1 and 3,061 cells mouse 2


# Normalize and find variable features
set.seed(123)
seur.combined <- NormalizeData(seur.combined) # normalized data is saved in seur.combined[["RNA"]]@data
seur.combined <- FindVariableFeatures(seur.combined) # return 2,000 HVGs
# plot variable features with labels
LabelPoints(plot = VariableFeaturePlot(seur.combined),
            points = head(VariableFeatures(seur.combined), 10), repel = TRUE)


# Run integration & dimension reduction
seur.combined <- RunFastMNN(object.list = SplitObject(seur.combined, split.by = "orig.ident"), k = 15)
seur.combined <- RunUMAP(seur.combined, reduction = "mnn", dims = 1:10, return.model=T, min.dist=0.3, spread=1)
seur.combined <- FindNeighbors(seur.combined, reduction = "mnn", dims = 1:10, k.param=15, compute.SNN = TRUE)
seur.combined <- FindClusters(seur.combined, resolution = 0.5)

# Cluster numbers are not in the same order as in the paper, so we'll just replace them
DimPlot(seur.combined, group.by = c("seurat_clusters"), pt.size = 0.1, label = T)




#### REPRODUCE FIGURES ####

# Figure 1B
FeaturePlot(seur.combined, features = c('Cd24a', 'Zbtb16', 'Cd44', 'Tbx21', 'Cxcr3', 'Itga1'),
            cols = c('#deebf7', 'red'), pt.size = 0.7, ncol = 3, order=T)
FeaturePlot(seur.combined, features = c('Rorc', 'Ccr6', 'Itga5', 'Itgb3', 'Hells', 'Cdk1', 'Mki67', 'Sell'),
            cols = c('#deebf7', 'red'), pt.size = 0.7, ncol = 4, order=T)




#### SIGNATURES WITH VISION ####
library(VISION)

# Get the MAIT1 and MAIT17 signatures
mait_signatures <- read.csv("~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/00_Reproduce_UMAPs/MAIT_Legoux_signatures.csv", header=T)

# Create VISION gene signature objects
mait1 <- unique(mait_signatures$MicroArray_thyMAIT1[!mait_signatures$MicroArray_thyMAIT1 == ""])
mait1_vector <- rep(1, length(mait1))
names(mait1_vector) <- mait1
mait1_sig <- createGeneSignature(name="MAIT1", sigData=mait1_vector)

mait17 <- unique(mait_signatures$MicroArray_thyMAIT17[!mait_signatures$MicroArray_thyMAIT17 == ""])
mait17_vector <- rep(1, length(mait17))
names(mait17_vector) <- mait17
mait17_sig <- createGeneSignature(name="MAIT17", sigData=mait17_vector)

mait_all_sig <- c(mait1_sig, mait17_sig)

# Run Vision
vision.obj.seurat <- Vision(GetAssayData(seur.combined, slot="data"), signatures = mait_all_sig,
                            projection_methods = "UMAP", meta = seur.combined@meta.data)

vision.obj.seurat <- analyze(vision.obj.seurat)

viewResults(vision.obj.seurat)

