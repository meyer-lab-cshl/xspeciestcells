# Purpose:  import Sagar mouse GD data
# Author: Salomé Carcy
# Date: June 2023


# **************
# 1. IMPORT ####
# **************

library(RaceID)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)

# Import data
sc <- readRDS("./data/raw_data/mouse_data/GDT_Sagar/adult_raceID_object.rds")
plotmap(sc, um=T)
plotsymbolsmap(sc,sub("(\\_|\\.).+","", colnames(sc@ndata)),fr=F,
               samples_col=RColorBrewer::brewer.pal(11, "Paired"))




# *************************
# 2. CONVERT TO SEURAT ####
# *************************

# Explore structure of RaceID object
sc@umap # UMAP embedding
sub("(\\_|\\.).+","", colnames(sc@ndata))[1:5] # sorted cell types
sc@expdata[1:5,1:5] # raw data


# Convert to seurat
counts <- sc@expdata
genes_with_pipes <- grep("\\|", rownames(sc@expdata), value=T)
# data <- sc@ndata
metadata <- data.frame("cell_id" = colnames(sc@expdata),
                       "cell_sort" = sub("(\\_|\\.).+","", colnames(sc@expdata)),
                       row.names = colnames(sc@expdata))
metadata$cell_kept_sagar <- ifelse(metadata$cell_id %in% colnames(sc@ndata), T, F)
metadata %>%
  mutate(cell.ident = replace(cell_sort, cell_sort=="DN1", "cKIT+ DN1"),
         cell.ident = replace(cell_sort, cell_sort=="DN3GDLO", "Pre-selected GD"),
         cell.ident = replace(cell_sort, cell_sort=="DN3GDHI", "Post-selected GD"),
         cell.ident = replace(cell_sort, cell_sort=="GDCD24NEG", "CD24- GD"),
         cell.ident = replace(cell_sort, cell_sort=="GDCD24POS", "CD24+ GD"),
         cell.ident = replace(cell_sort, cell_sort=="CD24NEG", "CD24- GD"),)
  

# Create seurat object
seur.raw <- CreateSeuratObject(counts=counts[!rownames(counts) %in% genes_with_pipes,],
                               meta.data=metadata)
print(seur.raw)

# add umap embedding
# seur.raw[["umap_sagar"]] <- CreateDimReducObject(embeddings=as.matrix(sc@umap), key="UMAP_", assay=DefaultAssay(seur.raw))




# *********************
# 2. PREPROCESSING ####
# *********************


# QC
# seur.raw[["percent.mt"]] <- PercentageFeatureSet(seur.raw, pattern = "^mt-") # no mitochondrial genes
FeatureScatter(seur.raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(seur.raw, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size=.01)
seur.filt <- subset(seur.raw, subset = nCount_RNA >= 2500) # 4122 cells left
table(seur.filt$cell_kept_sagar) # they are all kept from authors

# add umap embedding
seur.filt[["umap_sagar"]] <- CreateDimReducObject(embeddings=as.matrix(sc@umap[colnames(seur.filt),]),
                                                  key="UMAP_", assay=DefaultAssay(seur.filt))
seur.filt[["tsne_sagar"]] <- CreateDimReducObject(embeddings=as.matrix(sc@tsne[colnames(seur.filt),]),
                                                  key="TSNE_", assay=DefaultAssay(seur.filt))
DimPlot(seur.filt, reduction="umap_sagar")
