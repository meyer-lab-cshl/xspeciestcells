# Purpose: transform ENSEMBL id in seurat object to gene symbols
# Author: Sarah Chapin
# Date: March 2023

library(DropletUtils)
library(scater)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(Seurat)
# library(scMerge)
library(scran)
# library(harmony)
library(cowplot)
library(scales)
library(biomaRt)

#Load full seurat data
park_thymus_seurat <- readRDS("~/Projects/HumanThymusProject/data/human-thymus/HumanData_16_ParkIntegration/park_thymus_seurat.rds")

#Convert to sce
park_sce <- as.SingleCellExperiment(park_thymus_seurat)

##Annotate by gene name
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

genemap <- getBM(attributes=c('ensembl_gene_id', "external_gene_name"),
                 filters = 'ensembl_gene_id',
                 values = rownames(park_sce),
                 mart = ensembl) %>%
  dplyr::rename(gene_id = ensembl_gene_id,
                symbol = external_gene_name) %>%
  distinct()

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

# save sce
# saveRDS(park_sce_final, 
#         "~/20221005_human_multiome/data/park_full_sce_gene_names.rds")
#                              reductions = reducedDim(park_sce_final))

#convert back to seurat
park_seu_final <- as.Seurat(park_sce_final)

#save seu
saveRDS(park_seu_final,
        "~/Projects/HumanThymusProject/data/human-thymus/HumanData_16_ParkIntegration/park_thymus_seu_gene_names.rds")
