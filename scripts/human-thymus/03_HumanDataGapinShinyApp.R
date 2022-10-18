###
# Purpose: Create shiny interface on human (Gapin) data
# Date: Thur Oct 13th
# Author: Salomé Carcy
###


## INSTALL SHINY CELL ##
# code below for installation is from https://github.com/SGDDNB/ShinyCell

# Install required packages for ShinyCell
# reqPkg = c("data.table", "Matrix", "hdf5r", "reticulate", "ggplot2", 
#            "gridExtra", "glue", "readr", "RColorBrewer", "R.utils", "Seurat")
# newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
# if(length(newPkg)){install.packages(newPkg)}

# Install required packages for the Shiny app interface to work
# reqPkg = c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "hdf5r", 
#            "reticulate", "ggplot2", "gridExtra", "magrittr", "ggdendro")
# newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
# if(length(newPkg)){install.packages(newPkg)}

# Install Shiny cell
# devtools::install_github("SGDDNB/ShinyCell")



## MAKE SHINY APP FROM SEURAT OBJECT ##

library(Seurat)
library(ShinyCell)

# Import data
seur.human <- readRDS("~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/human_data/filtered_seurat_Harmony_07-22-22.RDS")

# Take quick peek at the data
# print(seur.human) # 79,801 cells and 34,778 genes
# seur.human@meta.data

# Create Shiny app
scConf = createConfig(seur.human)
makeShinyApp(seur.human, scConf,
             gex.assay="SCT", gex.slot="data",
             default.gene1="CD4",
             default.gene2="CD8A",
             shiny.title = "ShinyCell Gapin Human data",
             shiny.dir="shinyAppHumanGapinData/")

# Then, to run the shiny app, go into the "shinyAppHumanGapinData/" directory, open the "server.R" file into the Rstudio, and click on "Run App"