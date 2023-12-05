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
seur.human <- readRDS("~/Projects/HumanThymusProject/data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")

# Take quick peek at the data
# print(seur.human) # 78,607 cells and 17,204 genes
# seur.human@meta.data

# Create Shiny app
scConf = createConfig(seur.human)
makeShinyApp(seur, createConfig(seur),
             gex.assay="RNA", gex.slot="data",
             default.gene1="CD4",
             default.gene2="CD8A",
             default.dimred="UMAP_50",
             shiny.title = "ShinyCell Gapin Human data",
             shiny.dir="shinyAppHumanGapinData/")

# Then, to run the shiny app, go into the "shinyAppHumanGapinData/" directory, open the "server.R" file into the Rstudio, and click on "Run App"