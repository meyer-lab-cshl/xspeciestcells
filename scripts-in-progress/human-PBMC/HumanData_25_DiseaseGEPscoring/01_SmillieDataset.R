###
# Purpose: Score GEPs on disease datasets
# Date: July 2023
# Author: Salomé Carcy
###


# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)
library(Seurat)
library(RColorBrewer)

## 1.2. Data ####
seu <- readRDS("./data/raw_data/human_data/smillie-ulcerativecolitis/train.Imm.seur")
DimPlot(seu, reduction="tsne")
head(seu@meta.data)




# *******************
# 2. GET T CELLS ####
# *******************

table(seu@meta.data$Cluster, useNA="ifany")
tcells <- c("NKs", "ILCs", "CD4+ Activated Fos-hi", "CD4+ Activated Fos-lo", "CD4+ Memory", "Tregs",
            "CD4+ PD1+", "CD8+ IELs", "CD8+ LP", "CD8+ IL17+", "MT-hi", "Cycling T")

# Subset to t-cells
seuT <- subset(seu, subset=Cluster %in% tcells)
seuT@meta.data$Cluster <- factor(seuT@meta.data$Cluster, levels=unique(as.character(seuT$Cluster)))
DimPlot(seuT, reduction="tsne", group.by="Cluster")+
  scale_color_manual(values=RColorBrewer::brewer.pal(12, "Paired"))

# Re-run Tsne
# test <- RunTSNE(seuT, perplexity=25, dims=1:20, reduction.name="tsne_tcells", reduction.key="tSNEtcells_", max.iter=1000)
# DimPlot(test, reduction="tsne_tcells", group.by="Cluster")+
#   scale_color_manual(values=RColorBrewer::brewer.pal(12, "Paired"))
# Impossible to reproduce their T cell T-SNE




# ******************
# 3. SCORE GEPs ####
# ******************

# Import GEPs
gep_topgenes <- read.csv("./data/human-thymus/HumanData_17_GEPsOnParkData/genes_per_GEP_df_2023-04-07.csv", row.names=1)
dim(gep_topgenes) # 12 GEPs

# Keep only GEPs of interest
gep_pbmc <- c("GEP_1", "GEP_4", "GEP_5", "GEP_6", "GEP_8", "GEP_11", "GEP_12")
gep_topgenes <- gep_topgenes[,gep_pbmc]

# Get the genes into a list
gep_list <- list()
for (i in colnames(gep_topgenes)){
  print(i)
  # gep_allgenes <- gep_topgenes[1:200,i][!is.na(gep_topgenes[1:200,i])]
  gep_allgenes <- gep_topgenes[,i][!is.na(gep_topgenes[,i])]
  print(length(gep_allgenes))
  # Remove genes not present in seurat object
  gep_allgenes <- gep_allgenes[gep_allgenes %in% rownames(seuT)]
  print(length(gep_allgenes))
  gep_list[[i]] <- gep_allgenes
}


# Score these GEPs on our data & the UC data
seur <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
seur.geps   <- AddModuleScore(seur, name = names(gep_list), features=gep_list)
seurUC.geps <- AddModuleScore(seuT, name = names(gep_list), features=gep_list)

# remove the number added at the end of the score name
colnames(seur.geps@meta.data)[76:82]   <- str_sub(colnames(seur.geps@meta.data)[76:82], end=-2)
colnames(seurUC.geps@meta.data)[10:16] <- str_sub(colnames(seurUC.geps@meta.data)[10:16], end=-2)


# Plot on dimred
SCpubr::do_FeaturePlot(subset(seur.geps, subset=Tissue=="PBMC"), reduction="UMAP_50", features=gep_pbmc,
                       ncol=4, viridis_color_map = "B", order=T)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/smillie_umapgapin_scores.jpeg", width=35, height=20)
SCpubr::do_FeaturePlot(seurUC.geps, reduction="tsne", features=gep_pbmc,
                       ncol=4, viridis_color_map = "B", order=T)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/smillie_tsnesmillie_scores.jpeg", width=35, height=20)
cols_smillie <- RColorBrewer::brewer.pal(12, "Paired")
names(cols_smillie) <- unique(seurUC.geps$Cluster)
SCpubr::do_DimPlot(seurUC.geps, reduction="tsne", group.by = "Cluster", colors.use = cols_smillie)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/smillie_tsnesmillie_clusters.jpeg", width=7, height=7)


# Plot on VlnPlot
VlnPlot(seurUC.geps, group.by="Cluster", features=gep_pbmc, cols = cols_smillie, ncol = 4)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/smillie_vlnplotsmillie_scores.jpeg", width=15, height=10)
plot_grid(VlnPlot(seurUC.geps, group.by="Cluster", split.by="Health", features="GEP_1"),
          VlnPlot(seurUC.geps, group.by="Cluster", split.by="Health", features="GEP_4"),
          VlnPlot(seurUC.geps, group.by="Cluster", split.by="Health", features="GEP_5"),
          VlnPlot(seurUC.geps, group.by="Cluster", split.by="Health", features="GEP_6"),
          VlnPlot(seurUC.geps, group.by="Cluster", split.by="Health", features="GEP_11"),
          VlnPlot(seurUC.geps, group.by="Cluster", split.by="Health", features="GEP_12"),
          ncol=1)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/smillie_vlnplotsmillie_scores2.jpeg", width=8, height=30)

# Better comparison of healthy-uninflammed-inflammed
dfplot <- seurUC.geps@meta.data[,c("Cluster", "Health")]
dfplot <- dfplot %>%
  rownames_to_column("cellid") %>%
  as_tibble() %>%
  # get nb cells per group
  group_by(Cluster, Health) %>%
  count() %>%
  ungroup() %>%
  # get proportions
  group_by(Health) %>%
  mutate(totcells=sum(n)) %>%
  ungroup() %>%
  mutate(freq=n*100/totcells)

ggplot(dfplot, aes(x=Health, y=freq, fill=Cluster))+
  geom_bar(stat="identity", position="stack")+
  scale_fill_manual(values=cols_smillie)
ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/smillie_cellproportions.jpeg", width=8, height=8)
