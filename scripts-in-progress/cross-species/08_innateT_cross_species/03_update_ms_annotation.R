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
source("~/Projects/phd/scripts/colors_universal.R")

# Import data
seur.ms <- readRDS("./data/cross-species/08_innateT_cross_species/Analysis_all_mouse-Tinn_filtered_seurat_MNN.rds")





# *****************
# 2. FUNCTIONS ####
# *****************

plot_DotPlot <- function(seurobj, group, features, scaling=T){
  # get plot data
  p <- Seurat::DotPlot(
    seurobj,
    group.by=group,
    features=features,
    scale=scaling
  )
  # plot in personalized way
  p <- ggplot(p$data, aes(x=id, y=features.plot, fill=avg.exp.scaled, size=pct.exp))+
    geom_point(color="black", shape=21)+
    # scale_fill_gradient2(low=scales::muted("blue"), high=scales::muted("red"), name="z-score\nnormalized\navg expression")+
    scale_size_continuous(range=c(0,6), limits=c(0,100), name="%cells\nexpressing\ngene")+
    theme_bw()+
    theme(axis.text.y=element_text(face="italic"),
          axis.text.x=element_text(angle=45, hjust=1))+
    labs(y="", x="Clusters")
  # different color scale if scaled or not
  if(scaling==T){
    p <- p + scale_fill_gradient2(low=scales::muted("blue"), high=scales::muted("red"), name="z-score\naverage\nnorm expression")
  } else{
    p <- p + viridis::scale_fill_viridis(option="B", direction=-1, name="avg norm.\nexpression")
  }
  return(p)
}



# ****************
# 3. ANALYSIS ####
# ****************

#___________________________
## 3.1. Clean-up metadata a bit ####

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
## /end ####


#___________________________
## 3.2. Lineage composition per cluster ####
table(seur.ms@meta.data$Cell.type,useNA="ifany")

seur.ms@meta.data %>%
  as_tibble() %>%
  group_by(New_clusters, Cell.type) %>%
  count() %>%
  mutate(Cell.type=str_replace(Cell.type, "NKT", "iNKT")) %>%
  ggplot(aes(x=New_clusters, y=n, fill=Cell.type)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=cols_lineages)
# ggsave("~/Projects/HumanThymusProject/data/cross-species/08_innateT_cross_species/ms_integrated_cell_annotation/lineage_composition_barplot.jpeg")

## /end ####


#___________________________
## 3.3. Dotplot genes of interest ####

plot_DotPlot(
  seur.ms,
  group="New_clusters",
  features=c("Il2ra", "Clec12a", "Cd24a", "Sox13", "Cd200", "Nt5e", "Rorc", "Tbx21", "Cd44"),
  scaling=F
)
# ggsave("~/Projects/HumanThymusProject/data/cross-species/08_innateT_cross_species/ms_integrated_cell_annotation/genes_of_interest_dotplot.jpeg", width=5, height=4)

## /end ####


#___________________________
## 3.4. Fourth analysis ####

## /end ####

