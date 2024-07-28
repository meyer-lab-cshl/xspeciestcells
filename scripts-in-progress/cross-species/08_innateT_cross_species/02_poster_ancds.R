# Purpose:
# Author: Salomé Carcy
# Date:




# **************
# 1. IMPORT ####
# **************

# Import librairies
library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)
library(Seurat)
setwd("~/Projects/HumanThymusProject/")

# Import data
seur.ms <- readRDS("./data/cross-species/08_innateT_cross_species/Analysis_all_mouse-Tinn_filtered_seurat_MNN.rds")




# *****************
# 2. FUNCTIONS ####
# *****************





# ****************
# 3. ANALYSIS ####
# ****************

#___________________________
## 3.1. UMAP per dataset ####
head(seur.ms@meta.data)
table(seur.ms@meta.data$orig.ident)
table(seur.ms@meta.data$New.ident)

# rename some clusters
seur.ms@meta.data$datasets <- case_when(
  seur.ms@meta.data$New.ident=="Chandra_MAIT" ~ "(MAIT) Chandra et al.",
  seur.ms@meta.data$New.ident=="KMB_NKT" ~ "(iNKT) Maas-Bauer et al.",
  seur.ms@meta.data$New.ident=="Koay_MAIT" ~ "(MAIT) Koay et al.",
  seur.ms@meta.data$New.ident=="Krovi_NKT" ~ "(iNKT) Krovi et al.",
  seur.ms@meta.data$New.ident=="Lee_Gd" ~ "(GD) Lee et al.",
  seur.ms@meta.data$New.ident=="Lee_MAIT" ~ "(MAIT) Lee et al.",
  seur.ms@meta.data$New.ident=="Lee_NKT" ~ "(iNKT) Lee et al.",
  seur.ms@meta.data$New.ident=="Legoux_MAIT" ~ "(MAIT) Legoux et al.",
  seur.ms@meta.data$New.ident=="Li_Gd" ~ "(GD) Li et al.",
  seur.ms@meta.data$New.ident=="Paget_NKT" ~ "(iNKT) Paget et al.",
  seur.ms@meta.data$New.ident=="SAP_MAIT" ~ "(MAIT, SAP-/-) Legoux et al.",
  seur.ms@meta.data$New.ident=="Stage0_NKT" ~ "(iNKT) Wang et al.",
  seur.ms@meta.data$New.ident=="Stage1_NKT" ~ "(iNKT) Wang et al.",
  seur.ms@meta.data$New.ident=="Stage2_NKT" ~ "(iNKT) Wang et al.",
  seur.ms@meta.data$New.ident=="Stage3_NKT" ~ "(iNKT) Wang et al."
)
table(seur.ms@meta.data[,c("New.ident", "datasets")], useNA="ifany")

# color palette
colors_pal <- RColorBrewer::brewer.pal(12, "Paired")
# colors_pal <- c(
#   # GD colors
#   "#08306b",
#   "#2171b5",
#   # iNKT colors
#   "#3f007d",
#   "#6a51a3",
#   "#9e9ac8",
#   "#dadaeb",
#   "#fcfbfd",
#   # MAIT colors
#   # "#41b6c4",
#   # "#7fcdbb",
#   # "#c7e9b4",
#   # "#edf8b1",
#   # "#ffffd9"
#   "#6baed6",
#   "#9ecae1",
#   "#c6dbef",
#   "#deebf7",
#   "#f7fbff"
# )
names(colors_pal) <- sort(unique(seur.ms@meta.data$datasets))

# order MAIT & GD on top, because there are a lot of iNKT
# table(seur.ms@meta.data$Cell.type, useNA="ifany")
# cells_order <- c(
#   sample(rownames(seur.ms@meta.data[seur.ms@meta.data$Cell.type=="GD",]), size=6804, replace=FALSE),
#   sample(rownames(seur.ms@meta.data[seur.ms@meta.data$Cell.type=="MAIT",]), size=6151, replace=FALSE),
#   sample(rownames(seur.ms@meta.data[seur.ms@meta.data$Cell.type=="NKT",]), size=30535, replace=FALSE)
# )

ggrastr::rasterise(
  SCpubr::do_DimPlot(
    seur.ms,
    group.by="datasets",
    # order=rev(cells_order),
    shuffle=TRUE,
    legend.position="right",
    colors.use=colors_pal
  ),
  layers="Point",
  dpi=300
)
# ggsave("~/Desktop/Meyer-lab/Conferences/2024-07_ANDCS_Paris/fig3_crosspecies_ms_umap_datasets1.pdf", width=24, height=20, units="cm")
## /end ####


#___________________________
## 3.2. UMAP per cluster ####
seur.ms@meta.data$clusters_annotated <- case_when(
  seur.ms@meta.data$New_clusters == "0" ~ "post-selection",
  seur.ms@meta.data$New_clusters == "1" ~ "immature CD24+ GD",
  seur.ms@meta.data$New_clusters == "2" ~ "immature CD24+ GD",
  seur.ms@meta.data$New_clusters == "3" ~ "signaling",
  seur.ms@meta.data$New_clusters == "4" ~ "signaling",
  seur.ms@meta.data$New_clusters == "5" ~ "cycling",
  seur.ms@meta.data$New_clusters == "6" ~ "transition",
  seur.ms@meta.data$New_clusters == "7" ~ "typeII",
  seur.ms@meta.data$New_clusters == "8" ~ "typeI",
  seur.ms@meta.data$New_clusters == "9" ~ "typeI",
  seur.ms@meta.data$New_clusters == "10" ~ "typeIII",
  seur.ms@meta.data$New_clusters == "11" ~ "typeIII",
  seur.ms@meta.data$New_clusters == "12" ~ "12"
)


# colors
colors_pal_clusters <- RColorBrewer::brewer.pal(9, "Pastel1")
names(colors_pal_clusters) <- sort(unique(seur.ms@meta.data$clusters_annotated))

ggrastr::rasterise(
  SCpubr::do_DimPlot(
    seur.ms,
    group.by="clusters_annotated",
    legend.position="right",
    colors.use=colors_pal_clusters
  ),
  layers="Point",
  dpi=300
)
# ggsave("~/Desktop/Meyer-lab/Conferences/2024-07_ANDCS_Paris/fig3_crosspecies_ms_umap_annotation.pdf", width=24, height=20, units="cm")


## /end ####

