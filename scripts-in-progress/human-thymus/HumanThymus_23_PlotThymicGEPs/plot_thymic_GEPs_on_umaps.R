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

# Import GEP usage
gep_usage <- read.table("./data/human-PBMC/HumanData_20_RibbonPlotCellStateToID/cNMF_output/non_imputed_cNMF_allcells.usages.k_12.dt_0_02.consensus.txt", header=T)
dim(gep_usage)
colnames(gep_usage) <- paste0("GEP", c(2,5,3,1,4,12,6,7,8,10,9,11), "_usage") # rename column names
head(gep_usage)
gep_usage$GEP12_usage <- NULL



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
ggsave("./scripts-in-progress/human-thymus/HumanThymus_23_PlotThymicGEPs/plots/combined_umaps_innate.jpeg", width=20, height=35, limitsize=F)



## /end ####


#___________________________
## 3.4. Plot GEP usage per lineage in grid ####

# Split gep_usage to cells of interest
gep_usage_nkt  <- gep_usage[rownames(gep_usage) %in% rownames(seur.nkt@meta.data), ]
gep_usage_mait <- gep_usage[rownames(gep_usage) %in% rownames(seur.mait@meta.data), ]
gep_usage_gdt  <- gep_usage[rownames(gep_usage) %in% rownames(seur.gdt@meta.data), ]
gep_usage_cd4  <- gep_usage[rownames(gep_usage) %in% rownames(seur.cd4@meta.data), ]
gep_usage_cd8  <- gep_usage[rownames(gep_usage) %in% rownames(seur.cd8@meta.data), ]

# sanity checks
table(rownames(seur.nkt@meta.data) ==rownames(gep_usage_nkt),  useNA="ifany")
table(rownames(seur.mait@meta.data)==rownames(gep_usage_mait), useNA="ifany")
table(rownames(seur.gdt@meta.data) ==rownames(gep_usage_gdt),  useNA="ifany")
table(rownames(seur.cd4@meta.data) ==rownames(gep_usage_cd4),  useNA="ifany")
table(rownames(seur.cd8@meta.data) ==rownames(gep_usage_cd8),  useNA="ifany")

# Add GEP usage to seurat objects
seur.nkt@meta.data  <- cbind(seur.nkt@meta.data, gep_usage_nkt)
seur.mait@meta.data <- cbind(seur.mait@meta.data, gep_usage_mait)
seur.gdt@meta.data  <- cbind(seur.gdt@meta.data, gep_usage_gdt)
seur.cd4@meta.data  <- cbind(seur.cd4@meta.data, gep_usage_cd4)
seur.cd8@meta.data  <- cbind(seur.cd8@meta.data, gep_usage_cd8)

# Plot
gepsusage_of_interest <- paste0("GEP", 1:11, "_usage")
seur_vector <- list("NKT"=seur.nkt, "MAIT"=seur.mait, "GDT"=seur.gdt, "CD4"=seur.cd4, "CD8"=seur.cd8)

plist_massive <- list()
for(gep in gepsusage_of_interest){
  print(gep)
  # get row of plots (for one GEP)
  plist <- list()
  for(i_seur in names(seur_vector)){
    print(i_seur)
    seur <- seur_vector[[i_seur]]
    # get legend once
    if(gep=="GEP1_usage" & i_seur=="NKT"){
      plegend <- ggpubr::get_legend(
        SCpubr::do_FeaturePlot(seur,  features=gep, ncol=1, legend.position="right", legend.title="GEP usage", border.size=1, pt.size=2, order=T)+
          scale_color_viridis_c(limits=c(0,1.25), option="B")+
          theme(plot.background = element_rect(color = "black"))
      )
    }
    # make plot
    legendpos <- "none"
    plt_title <- ""
    if(gep=="GEP1_usage"){plt_title=i_seur}
    p <- SCpubr::do_FeaturePlot(seur,  features=gep, ncol=1, legend.position=legendpos, border.size=1, pt.size=2, order=T, plot.title=plt_title)+
      scale_color_viridis_c(limits=c(0,1.25), option="B")+
      theme(plot.background = element_rect(color = "black", linewidth = 2),
            plot.margin=unit(c(0,5.5,0,5.5), "points"))
    plist[[i_seur]] <- ggrastr::rasterise(p, layers="Point", dpi=300)
  }
  prow_gep <- plot_grid(plotlist=plist, nrow=1, scale=0.9)
  plist_massive[[gep]] <- prow_gep 
}

plot_grid(plot_grid(plotlist=plist_massive, ncol=1, align="h"),
          ggpubr::as_ggplot(plegend), ncol=2, rel_widths = c(5, .5), scale=c(0.95, 5))
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_23_PlotThymicGEPs/plots/combined_umaps_gepusage3.pdf", width=25, height=40, limitsize=F, dpi=80)


## /end ####

