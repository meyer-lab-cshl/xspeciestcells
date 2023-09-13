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
library(scico)
source("./scripts-final/colors_universal.R")

# Import data
seur.nkt  <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.nkt.RDS")
seur.mait <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.mait.RDS")
seur.gdt  <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.gd.RDS")
seur.cd4  <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.cd4.RDS")
seur.cd8  <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.cd8.RDS")

DimPlot(seur.nkt,  group.by="cell_annot", label=T)
DimPlot(seur.mait, group.by="cell_annot", label=T)
DimPlot(seur.gdt,  group.by="cell_annot", label=T)
DimPlot(seur.cd4,  group.by="cell_annot", label=T)
DimPlot(seur.cd8,  group.by="cell_annot", label=T)

# Cleanup a bit
seur.nkt@meta.data[,33:75]  <- NULL
seur.mait@meta.data[,33:75] <- NULL
seur.gdt@meta.data[,33:75]  <- NULL
seur.cd4@meta.data[,33:75]  <- NULL
seur.cd8@meta.data[,33:75]  <- NULL

# Import GEPs
cnmf_gapin  <- read.csv("./data/human-PBMC/HumanData_27_ComparisonGEPsGarner/limited_nonimputed_genes_per_gep_post_rank_threshold_k12.csv", row.names=1)
head(cnmf_gapin)



# *****************
# 2. FUNCTIONS ####
# *****************





# ****************
# 3. ANALYSIS ####
# ****************

#___________________________
## 3.1. Score GEPs on seurat objects ####

dflong_gapin <- gather(cnmf_gapin, key=geneprogram, value=gene, colnames(cnmf_gapin)) %>%
  filter(!is.na(gene)) %>%
  filter(geneprogram != "GEP12")
head(dflong_gapin)
table(dflong_gapin$geneprogram, useNA="ifany")
table(is.na(dflong_gapin$gene))

geneprograms_list <- list()
for(gp in unique(dflong_gapin$geneprogram)){
  print(gp)
  geneprograms_list[[gp]] <- dflong_gapin %>%
    filter(geneprogram==gp & gene %in% rownames(seur.nkt)) %>%
    pull(gene)
}
print(lengths(geneprograms_list))

# Compute cell scores
seur.nkt   <- AddModuleScore(seur.nkt,  name = names(geneprograms_list), features=geneprograms_list, seed=1)
seur.mait  <- AddModuleScore(seur.mait,  name = names(geneprograms_list), features=geneprograms_list, seed=1)
seur.gdt   <- AddModuleScore(seur.gdt,  name = names(geneprograms_list), features=geneprograms_list, seed=1)
seur.cd4   <- AddModuleScore(seur.cd4,  name = names(geneprograms_list), features=geneprograms_list, seed=1)
seur.cd8   <- AddModuleScore(seur.cd8,  name = names(geneprograms_list), features=geneprograms_list, seed=1)


# Remove the annoying numbers that are being added
colnames(seur.nkt@meta.data)[35:45] <- paste0("GEP", 1:11)
colnames(seur.mait@meta.data)[35:45] <- paste0("GEP", 1:11)
colnames(seur.gdt@meta.data)[35:45] <- paste0("GEP", 1:11)
colnames(seur.cd4@meta.data)[35:45] <- paste0("GEP", 1:11)
colnames(seur.cd8@meta.data)[35:45] <- paste0("GEP", 1:11)


## /end ####


#___________________________
## 3.2. Plot all GEPs per lineage ####
ggsave(filename = "./scripts-in-progress/human-thymus/HumanThymus_23_PlotThymicGEPs/plots/allgeps_umaps_nkt.jpeg",
       SCpubr::do_FeaturePlot(seur.nkt, features=paste0("GEP", 1:11), ncol=3),
       width=12, height=18)
ggsave(filename = "./scripts-in-progress/human-thymus/HumanThymus_23_PlotThymicGEPs/plots/allgeps_umaps_mait.jpeg",
       SCpubr::do_FeaturePlot(seur.mait, features=paste0("GEP", 1:11), ncol=3),
       width=12, height=18)
ggsave(filename = "./scripts-in-progress/human-thymus/HumanThymus_23_PlotThymicGEPs/plots/allgeps_umaps_gdt.jpeg",
       SCpubr::do_FeaturePlot(seur.gdt, features=paste0("GEP", 1:11), ncol=3),
       width=12, height=18)
ggsave(filename = "./scripts-in-progress/human-thymus/HumanThymus_23_PlotThymicGEPs/plots/allgeps_umaps_cd4.jpeg",
       SCpubr::do_FeaturePlot(seur.cd4, features=paste0("GEP", 1:11), ncol=3),
       width=12, height=18)
ggsave(filename = "./scripts-in-progress/human-thymus/HumanThymus_23_PlotThymicGEPs/plots/allgeps_umaps_cd8.jpeg",
       SCpubr::do_FeaturePlot(seur.cd8, features=paste0("GEP", 1:11), ncol=3),
       width=12, height=18)
## /end ####


#___________________________
## 3.3. Plot thymic GEPs per lineage in grid ####
geps_of_interest <- c("GEP9", "GEP1", "GEP2", "GEP3", "GEP4", "GEP5", "GEP6", "GEP7")
p_nkt  <- SCpubr::do_FeaturePlot(seur.nkt, features=geps_of_interest, ncol=1, legend.position="right", border.size=1, pt.size=3) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0, end = 1, direction = -1)
p_mait <- SCpubr::do_FeaturePlot(seur.mait, features=geps_of_interest, ncol=1, legend.position="right", border.size=1, pt.size=3) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0, end = 1, direction = -1)
p_gdt  <- SCpubr::do_FeaturePlot(seur.gdt, features=geps_of_interest, ncol=1, legend.position="right", border.size=1, pt.size=3) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0, end = 1, direction = -1)

p_nkt | p_mait | p_gdt
ggsave("./scripts-in-progress/human-thymus/HumanThymus_23_PlotThymicGEPs/plots/combined_umaps_innate.jpeg", width=18, height=35, limitsize=F)



## /end ####


#___________________________
## 3.4. Fourth analysis ####

## /end ####

