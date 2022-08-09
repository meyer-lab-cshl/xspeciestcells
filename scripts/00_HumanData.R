###
# Purpose: Import & plot human iNKT data (curtesy of Laurent Gapin)
# Date: Aug 9th 2022
# Author: Salomé Carcy
###


#### IMPORT ####

# Libraries
library(ggplot2)

# Data
path.data <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/human_data"
human5 <- read.csv(file.path(path.data, "CUThy13_220225_SampleTag05_hs_NKT_RSEC_MolsPerCell.csv"), sep=",", header=T, skip=7)
human8 <- read.csv(file.path(path.data, "CUTHY11BDRscRNA_seq_091621_SampleTag08_hs_NKT_RSEC_MolsPerCell.csv"), sep=",", header=T, skip=8)
human12 <- read.csv(file.path(path.data, "CUTHY12BDRscRNA_seq_211101_SampleTag12_hs_NKT_RSEC_MolsPerCell.csv"), sep=",", header=T, skip=8)

# Any cluster ID in the colnames?
# table(stringr::str_detect(colnames(human1),"[[:lower:]]"))
# colnames(human1)[stringr::str_detect(colnames(human1),"[[:lower:]]")]
