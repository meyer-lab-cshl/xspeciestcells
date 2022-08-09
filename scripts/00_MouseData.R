###
# Purpose: Import & plot mouse iNKT data (curtesy of Laurent Gapin)
# Date: Aug 9th 2022
# Author: Salomé Carcy
###


#### IMPORT ####

# Libraries
library(ggplot2)
library(Seurat)

# Data
path.data <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/mouse_data"
mouse1 <- Read10X(file.path(path.data, "B6_1"))


# cellular barcodes present in dataset
readr::read_tsv(file.path(path.data, "B6_1/barcodes.tsv.gz"), col_names = FALSE) # 8,373 cells
# IDs of quantified genes
readr::read_tsv(file.path(path.data, "B6_1/features.tsv.gz"), col_names = FALSE) # 31,053 genes
