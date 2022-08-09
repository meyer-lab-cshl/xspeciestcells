###
# Purpose: Import & plot pig iNKT data (curtesy of ...)
# Date: Aug 9th 2022
# Author: Salomé Carcy
###


#### IMPORT ####

# Libraries
library(ggplot2)
library(Seurat)

# Data
path.data <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/pig_data/GSE192520_RAW"
pig <- Read10X(path.data)


# cellular barcodes present in dataset
readr::read_tsv(file.path(path.data, "GSM5750559_thymus_barcodes.tsv.gz"), col_names = FALSE) # 9,112 cells
# IDs of quantified genes
readr::read_tsv(file.path(path.data, "GSM5750559_thymus_features.tsv.gz"), col_names = FALSE) # 21,301 genes
