# Purpose: make Park seurat object presentable to share with Laurent
# Author: Salomé Carcy
# Date: March 2023


# **************
# 1. IMPORT ####
# **************

# Import libraries
library(Seurat)
library(ggplot2)
# library(ggpointdensity)
# library(cowplot)
library(tidyverse)
library(dplyr)
# library(RColorBrewer)
# library(pals)
library(scater)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(Seurat)
library(scMerge)
library(scran)
# library(harmony)
library(scales)
library(biomaRt)


# Import Park seurat object
# Downloading the seurat object from [here](https://cellxgene.cziscience.com/collections/de13e3e2-23b6-40ed-a413-e9e12d7d3910)
park_thymus_seurat <- readRDS("data/human-thymus/HumanData_16_ParkIntegration/park_thymus_seurat.rds")




# *******************************************
# 2. CONVERT ENSEMBL IDs TO GENE SYMBOLS ####
# *******************************************

# Convert to sce object
park_sce <- as.SingleCellExperiment(park_thymus_seurat)

# Get Human ENSEMBL - gene symbol correspondance from biomart
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

genemap <- getBM(attributes=c('ensembl_gene_id', "external_gene_name"),
                 filters = 'ensembl_gene_id',
                 values = rownames(park_sce),
                 mart = ensembl) %>%
  dplyr::rename(gene_id = ensembl_gene_id,
                symbol = external_gene_name) %>%
  distinct()

# Change EMSEMBL IDs to gene symbols in SCE object
featureData <- tibble(gene_id=rownames(park_sce)) %>%
  left_join(genemap, by="gene_id") %>%
  drop_na() %>%
  filter(duplicated(symbol) == FALSE) %>%
  filter(symbol != "")
symbol <- dplyr::select(featureData, symbol) 
gene_id <- dplyr::select(featureData, gene_id) 

park_sce_final <- 
  subset(park_sce,
         rownames(park_sce) %in% gene_id$gene_id,
  )

mcols(park_sce_final) <- 
  DataFrame(mcols(park_sce_final), featureData)
rownames(park_sce_final) <- rowData(park_sce_final)$symbol

# convert back to seurat
park_seu_final <- as.Seurat(park_sce_final)




# *************************
# 3. ADD PARK CLUSTERS ####
# *************************

# The clusters saved in this seurat object are not as precise as the ones in the paper
DimPlot(park_seu_final, group.by="cell_type", label=T)

# Add metadata column (will be useful for integration)
park_seu_final@meta.data$study <- "park_data"

# I downloaded the AnnData object from the [developmental cell atlas](https://developmental.cellatlas.io/thymus-development)
# And then I exported the metadata from the anndata object into a .csv file (to get the more precise clusters)
# Import metadata
park_metadata <- read.csv("~/Projects/HumanThymusProject/data/human-thymus/HumanData_16_ParkIntegration/park_metadata.csv")

# Check that the row names are the same
table(park_metadata$index == rownames((park_seu_final@meta.data))) # yessss nice

# Add the cell types (cell clusters) to seurat object
park_seu_final@meta.data$cell_type <- park_metadata$cell.types

# Visualize
cols_park <- c("DN(early)" = "#78c679",
               "DN(P)" = "#41ab5d",
               "DN(Q)" = "#238443",
               "γδT" = "#92c051",
               "DP(P)" = "#b75347",
               "DP(Q)" = "#d8443c",
               "αβT(entry)" = "#e09351",
               "CD8+T"= "#5a97c1",
               "CD8αα(I)" = "#421401",
               "CD8αα(II)" = "#0a2e57",
               "CD4+T"= "gold",
               "T(agonist)" = "#9f5691",
               "Treg(diff)" = "#9f5691",
               "Treg" = "blueviolet",
               "Th17" = "#a40000",
               "NKT" = "#72bcd5")
DimPlot(park_seu_final, group.by="cell_type", label=T)+
  scale_color_manual(values=cols_park) # much better!!

# SCpubr::do_DimPlot(seur.park.thymocytes,
#                    group.by = "cell_type",
#                    pt.size=0.5,
#                    label=T,
#                    label.color="black",
#                    legend.position = "none",
#                    repel=T,
#                    colors.use=cols_park,
#                    font.size = 24)
# ggsave("data/human-thymus/HumanData_16_ParkIntegration/park_umap2.jpeg", width=10, height=9)




# **************************
# 4. SAVE SEURAT OBJECT ####
# **************************

saveRDS(park_seu_final, "data/human-thymus/HumanData_16_ParkIntegration/park_thymus_seu_gene_names.rds")
