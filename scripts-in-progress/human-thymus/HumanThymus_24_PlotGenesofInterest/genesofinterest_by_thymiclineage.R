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

# Import human data
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


# Import mouse data
seur.nkt.ms <- readRDS("./data/cross-species/00_Reproduce_UMAPs/ms_nkt_seurobj.rds")
seur.mait.ms <- readRDS("./data/cross-species/00_Reproduce_UMAPs/ms_mait_seurobj.rds")
ortholog.df <- read.csv("./data/cross-species/03_BiomartTable/big_ass_ortholog_table.csv")
DimPlot(seur.nkt.ms, group.by="cell_type", label=T)
DimPlot(seur.mait.ms, group.by="cell_type", label=T)



# *****************
# 2. FUNCTIONS ####
# *****************

plot_four_genes <- function(seur, genelist, ordercells=F){
  plist <- list()
  # identify which gene has highest max count
  max_per_gene <- lapply(genelist, function(x) max(seur@assays$RNA@data[x,]))
  names(max_per_gene) <- genelist
  gene_with_max_value <- names(which.max(max_per_gene))

  for (gene in genelist){
    print(gene)
    if(gene=="ZBTB16"){ordercells=T} # put PLZF at front
    # get legend once
    if(gene==gene_with_max_value){
      plegend <- ggpubr::get_legend(
        SCpubr::do_FeaturePlot(seur, features=gene, ncol=1, legend.position="right", legend.title="", border.size=1, pt.size=2, order=T) &
          scale_colour_scico(palette = "lapaz", alpha = 0.8, direction = -1)
      )
    }
    # featureplot
    p <- SCpubr::do_FeaturePlot(seur, features=gene, order=ordercells, ncol=1, legend.position="none", border.size=1, pt.size=3) +
      labs(title=gene)+
      theme(plot.background = element_rect(color = "black"),
            plot.title = element_text(hjust = 0.5, vjust=0.1, size=20)) &
      scale_colour_scico(palette = "lapaz", alpha = 0.8, direction = -1)
    plist[[gene]] <- ggrastr::rasterise(p, layers="Point", dpi=300)
  }
  
  
  # combine featureplots
  pcombined <- plot_grid(
    plot_grid(plotlist=plist, ncol=2, scale=1),
    ggpubr::as_ggplot(plegend),
    ncol=2, rel_widths = c(5, .5), scale=c(0.95, 0.5)
    )
  
  return(pcombined)
}

plot_genesignature <- function(seur, genesignature, ordercells=F){
  # get legend
  plegend <- ggpubr::get_legend(
    SCpubr::do_FeaturePlot(seur, features=genesignature, legend.position="right", border.size=1, pt.size=2, order=ordercells)+
      scale_color_viridis_c(limits=c(0,max(seur@meta.data[,genesignature])), option="B")
  )
  # get featureplot
  p <- SCpubr::do_FeaturePlot(seur, features=genesignature, legend.position="none", border.size=1, pt.size=2, order=ordercells)+
    scale_color_viridis_c(limits=c(0,max(seur@meta.data[,genesignature])), option="B")+
    theme(plot.background = element_rect(color = "black", linewidth = 2))
  p <- ggrastr::rasterise(p, layers="Point", dpi=300)
  # combine
  pcombined <- plot_grid(
    p,
    ggpubr::as_ggplot(plegend) + theme(plot.margin=margin(0,50,0,30)),
    ncol=2, rel_widths = c(5, 2), scale=c(0.95, 0.5)
  )
  return(pcombined)
}




# ****************
# 3. ANALYSIS ####
# ****************

#___________________________
## 3.1. Score gene signatures ####
gene_signatures <- list("effector_new"=c("HOPX", "GZMB", "NKG7", "TBX21", "PRF1", "GZMA", "KLRD1", "KLF6",
                                     "CCR6", "RORC", "TMEM176A", "TMEM176B", "JUNB", "FOS", "RORA", "FOSB"),
                        "effector"= c("HOPX", "GZMB", "GZMK", "ZEB2", "NKG7", "GNLY", "TBX21", "EOMES", "TYROBP", "PRF1",
                                      "CCL4", "CCL5", "KLRB1", "GZMH", "GZMA", "KLRD1", "CST7", "KLF6", "CXCR4"),
                        "naive"=c("SATB1", "TCF7", "LEF1", "CCR7", "SELL", "MYC", "EIF3E", "SOX4", "ID3", "BACH2"),
                        "egress"=c("KLF2", "CORO1A", "CCR7", "CXCR4", "CXCR6", "FOXO1", "CXCR3", "S1PR1", "S1PR4",
                                   "S100A4", "S100A6", "EMP3"))

seur.nkt   <- AddModuleScore(seur.nkt,  name = names(gene_signatures), features=gene_signatures, seed=1)
seur.mait  <- AddModuleScore(seur.mait, name = names(gene_signatures), features=gene_signatures, seed=1)
seur.gdt   <- AddModuleScore(seur.gdt,  name = names(gene_signatures), features=gene_signatures, seed=1)
seur.cd4   <- AddModuleScore(seur.cd4,  name = names(gene_signatures), features=gene_signatures, seed=1)
seur.cd8   <- AddModuleScore(seur.cd8,  name = names(gene_signatures), features=gene_signatures, seed=1)

colnames(seur.nkt@meta.data)[35:37]  <- names(gene_signatures)
colnames(seur.mait@meta.data)[35:37] <- names(gene_signatures)
colnames(seur.gdt@meta.data)[35:37]  <- names(gene_signatures)
colnames(seur.cd4@meta.data)[35:37]  <- names(gene_signatures)
colnames(seur.cd8@meta.data)[35:37]  <- names(gene_signatures)

# Score same signatures in mouse nkt/mait
gene_signatures_ms <- list("effector_new"=c("Hopx", "Gzmb", "Nkg7", "Tbx21", "Prf1", "Gzma", "Klrd1", "Klf6",
                                        "Ccr6", "Rorc", "Tmem176a", "Tmem176b", "Junb", "Fos", "Rora", "Fosb"),
                           "effector"=c("Hopx", "Gzmb", "Gzmk", "Zeb2", "Nkg7", "Tbx21", "Eomes", "Tyrobp", "Prf1",
                                        "Ccl4", "Ccl5", "Klrb1a", "Gzmg", "Gzma", "Klrd1", "Cst7", "Klf6", "Cxcr4"), # manual curating
                           "naive"=ortholog.df[ortholog.df$hu_symbol %in% gene_signatures$naive, "ms_symbol_data"],
                           "egress"=ortholog.df[ortholog.df$hu_symbol %in% gene_signatures$egress, "ms_symbol_data"])

seur.nkt.ms   <- AddModuleScore(seur.nkt.ms,  name = names(gene_signatures_ms), features=gene_signatures_ms, seed=1)
seur.mait.ms  <- AddModuleScore(seur.mait.ms, name = names(gene_signatures_ms), features=gene_signatures_ms, seed=1)

colnames(seur.nkt.ms@meta.data)[8:10]  <- names(gene_signatures_ms)
colnames(seur.mait.ms@meta.data)[10:12] <- names(gene_signatures_ms)
## /end ####


#___________________________
## 3.2. Plot gene signatures ####

# human
plot_genesignature(seur.nkt, genesignature = "naive", ordercells=F)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_naivesig.jpeg", width=7, height=5)
plot_genesignature(seur.mait, genesignature = "naive", ordercells=F)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_naivesig.jpeg", width=7, height=5)

plot_genesignature(seur.nkt, genesignature = "effector_new", ordercells=F)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_effectorsigSC.jpeg", width=7, height=5)
plot_genesignature(seur.mait, genesignature = "effector_new", ordercells=F)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_effectorsigSC.jpeg", width=7, height=5)
plot_genesignature(seur.gdt, genesignature = "effector_new", ordercells=F)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_gdt_effectorsigSC.jpeg", width=7, height=5)



# mouse
plot_genesignature(seur.nkt.ms, genesignature = "naive")
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_naivesig_mouse.jpeg", width=7, height=5)
plot_genesignature(seur.mait.ms, genesignature = "naive")
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_naivesig_mouse.jpeg", width=7, height=5)

plot_genesignature(seur.nkt.ms, genesignature = "effector")
plot_genesignature(seur.mait.ms, genesignature = "effector")


## /end ####



#___________________________
## 3.3. Plot genes of interest ####

# naive genes
plot_four_genes(seur.nkt, c("CCR9", "CCR7", "SELL", "TCF7"))
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_naivegenes.jpeg", width=6, height=5)
plot_four_genes(seur.mait, c("CCR9", "CCR7", "SELL", "TCF7"))
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_naivegenes.jpeg", width=6, height=5)


# effector genes
plot_four_genes(seur.nkt, c("ZBTB16", "EOMES", "KLRB1", "GZMK"))
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_effectorgenes.jpeg", width=6, height=5)
plot_four_genes(seur.mait, c("ZBTB16", "EOMES", "KLRB1", "GZMK"))
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_effectorgenes.jpeg", width=6, height=5)
plot_four_genes(seur.gdt, c("KLRD1", "EOMES", "KLRB1", "GZMK"), ordercells=T)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_gdt_effectorgenes.jpeg", width=6, height=5)


# CD4 and CD8
FeaturePlot(seur.nkt, features=c("CD4", "CD8A"), blend=T, blend.threshold=0.01, cols=c("lightgrey", "red", "blue"), order=T, pt.size=2)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_cd4cd8.jpeg", width=18, height=5)
FeaturePlot(seur.mait, features=c("CD4", "CD8A"), blend=T, cols=c("lightgrey", "red", "blue"), order=T, pt.size=2)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_cd4cd8.jpeg", width=18, height=5)
FeaturePlot(seur.gdt, features=c("CD4", "CD8A"), blend=T, cols=c("lightgrey", "red", "blue"), order=T, pt.size=2)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_gdt_cd4cd8.jpeg", width=18, height=5)


plot_four_genes(seur.nkt, c("CD4", "CD8A"))
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_effectorgenes.jpeg", width=6, height=3)
plot_four_genes(seur.mait, c("ZBTB16", "EOMES", "GZMK", "KLRB1"))
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_effectorgenes.jpeg", width=6, height=5)
plot_four_genes(seur.gdt, c("ZBTB16", "EOMES", "GZMK", "KLRB1"))
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_gdt_effectorgenes.jpeg", width=6, height=5)




## /end ####


#___________________________
## 3.2. Second analysis ####

## /end ####


#___________________________
## 3.3. Third analysis ####

## /end ####


#___________________________
## 3.4. Fourth analysis ####

## /end ####

