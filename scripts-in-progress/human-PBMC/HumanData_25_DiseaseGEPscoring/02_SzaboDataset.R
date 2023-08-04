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
library(SingleCellExperiment)
library(RColorBrewer)
library(biomaRt)


## 1.2. Data ####
seu.szabo <- readRDS("./data/raw_data/human_data/szabo-Tstim/local.rds")

# Quick visualization
DimPlot(seu.szabo, reduction="umap", group.by="donor_id")
# Donor 1 and 2 have lung/BM, lymph node (from deceased donors)
# Donor A and B have blood samples (from alive donors)
DimPlot(seu.szabo, reduction="umap", group.by="stimulation_status")
DimPlot(seu.szabo, reduction="umap", group.by="tissue")
DimPlot(seu.szabo, reduction="umap", group.by="cell_type")
DimPlot(seu.szabo, reduction="umap", group.by="cd4cd8_status")

# head(seu.szabo@meta.data)

# Change some metadata
seu.szabo@meta.data$tissue <- as.character(seu.szabo@meta.data$tissue)
seu.szabo@meta.data[seu.szabo@meta.data$tissue=="bone marrow", "tissue"] <- "BM"
seu.szabo@meta.data[seu.szabo@meta.data$tissue=="lower lobe of left lung", "tissue"] <- "LG"
seu.szabo@meta.data[seu.szabo@meta.data$tissue=="venous blood", "tissue"] <- "BL"
seu.szabo@meta.data[seu.szabo@meta.data$tissue=="bronchopulmonary lymph node", "tissue"] <- "LN"

seu.szabo@meta.data$stimulation_status <- as.character(seu.szabo@meta.data$stimulation_status)
seu.szabo@meta.data[seu.szabo@meta.data$stimulation_status=="act",  "stimulation_status"] <- "activated"
seu.szabo@meta.data[seu.szabo@meta.data$stimulation_status=="rest", "stimulation_status"] <- "resting"

temp <- as.data.frame(seu.szabo@meta.data[,c("tissue", "stimulation_status")]) %>%
  mutate(tissue_source=paste(tissue, stimulation_status))
seu.szabo@meta.data$tissue_source <- temp$tissue_source


## 1.3. Transform ENSBL to gene symbols ####
rownames(seu.szabo)[1:5]
# Get ENSEMBL-gene symbol dict
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
genemap <- getBM(attributes=c('ensembl_gene_id', "external_gene_name"),
                 filters = 'ensembl_gene_id',
                 values = rownames(seu.szabo),
                 mart = ensembl) %>%
  dplyr::rename(gene_id = ensembl_gene_id,
                symbol = external_gene_name) %>%
  distinct()

# Verify that these correspond to rownames(seu.szabo)
table(genemap$gene_id %in% rownames(seu.szabo), useNA="ifany")
table(rownames(seu.szabo) %in% genemap$gene_id, useNA="ifany") # 126 ENSEMBL were not found

# Check if there are any duplicated ENSEMBL or gene symbols
table(duplicated(genemap$gene_id), useNA="ifany") # none
table(duplicated(genemap$symbol), useNA="ifany") # 39,716 symbols are unique (the rest are duplicated)

# Check if any GEP gene is present in the duplicated symbols (otherwise we'll just remove them)
# Import GEPs
gep_topgenes <- read.csv("./data/human-thymus/HumanData_17_GEPsOnParkData/genes_per_GEP_df_2023-04-07.csv", row.names=1)
dim(gep_topgenes) # 12 GEPs
# Keep only GEPs of interest
gep_pbmc <- c("GEP_1", "GEP_4", "GEP_5", "GEP_6", "GEP_8", "GEP_11", "GEP_12")
gep_topgenes <- gep_topgenes[,gep_pbmc]
# Check if any GEP gene is present in the duplicated symbols
duplicatedgenes <- genemap$symbol[duplicated(genemap$symbol)==T]
gepgenes <- unlist(as.vector(gep_topgenes)[!is.na(as.vector(gep_topgenes))])
table(duplicatedgenes %in% gepgenes) # only 2 of them are present, so let's just get rid of them

# Keep only unique symbols
genemap <- genemap %>%
  filter(duplicated(symbol)==FALSE) # 39,716 genes with unique symbol

# Subset seurat object to the genes of interest
seu.szabo <- seu.szabo[genemap$gene_id,]
table(rownames(seu.szabo) == genemap$gene_id, useNA="ifany")

# Replace gene names
sce <- as.SingleCellExperiment(seu.szabo)
rowData(sce) <- genemap
rownames(sce) <- rowData(sce)$symbol
seu.szabo.new <- seu.szabo
seu.szabo.new@assays$RNA@counts <- assays(sce)[[1]]
seu.szabo.new@assays$RNA@data <- assays(sce)[[2]]




# ******************
# 2. SCORE GEPs ####
# ******************


# Get the genes into a list
gep_list <- list()
for (i in colnames(gep_topgenes)){
  print(i)
  # gep_allgenes <- gep_topgenes[1:200,i][!is.na(gep_topgenes[1:200,i])]
  gep_allgenes <- gep_topgenes[,i][!is.na(gep_topgenes[,i])]
  print(length(gep_allgenes))
  # Remove genes not present in seurat object
  gep_allgenes <- gep_allgenes[gep_allgenes %in% rownames(seu.szabo.new)]
  print(length(gep_allgenes))
  gep_list[[i]] <- gep_allgenes
}


# Score these GEPs on our data & the UC data
seur <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
seur.geps       <- AddModuleScore(seur, name = names(gep_list), features=gep_list)
seur.szabo.geps <- AddModuleScore(seu.szabo.new, name = names(gep_list), features=gep_list)

# remove the number added at the end of the score name
colnames(seur.geps@meta.data)[76:82]   <- str_sub(colnames(seur.geps@meta.data)[76:82], end=-2)
colnames(seur.szabo.geps@meta.data)[27:33] <- str_sub(colnames(seur.szabo.geps@meta.data)[27:33], end=-2)


# Plot on dimred
SCpubr::do_FeaturePlot(subset(seur.geps, subset=Tissue=="PBMC"), reduction="UMAP_50", features=gep_pbmc,
                       ncol=4, viridis_color_map = "B", order=T)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/szabo_umapgapin_scores.jpeg", width=35, height=20)
SCpubr::do_FeaturePlot(seur.szabo.geps, reduction="umap", features=gep_pbmc,
                       ncol=4, viridis_color_map = "B", order=T)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/szabo_umapszabo_scores.jpeg", width=35, height=20)


# Plot on dimred, split by tissue & activation
ggsave(plot=SCpubr::do_FeaturePlot(seur.szabo.geps, split.by="tissue_source", ncol=2, reduction="umap", features="GEP_1", viridis_color_map = "B", order=T),
       filename="./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/szabo_umapszabo_gep1.jpeg", width=10, height=20)
ggsave(plot=SCpubr::do_FeaturePlot(seur.szabo.geps, split.by="tissue_source", ncol=2, reduction="umap", features="GEP_4", viridis_color_map = "B", order=T),
       filename="./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/szabo_umapszabo_gep4.jpeg", width=10, height=20)
ggsave(plot=SCpubr::do_FeaturePlot(seur.szabo.geps, split.by="tissue_source", ncol=2, reduction="umap", features="GEP_5", viridis_color_map = "B", order=T),
       filename="./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/szabo_umapszabo_gep5.jpeg", width=10, height=20)
ggsave(plot=SCpubr::do_FeaturePlot(seur.szabo.geps, split.by="tissue_source", ncol=2, reduction="umap", features="GEP_6", viridis_color_map = "B", order=T),
       filename="./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/szabo_umapszabo_gep6.jpeg", width=10, height=20)


# Plot umap with different metadata
p1 <- DimPlot(seu.szabo, reduction="umap", group.by="stimulation_status")+
  scale_color_manual(values=c("#cb181d", "#2171b5"))+
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.line = element_blank())
p2 <- DimPlot(seu.szabo, reduction="umap", group.by="cd4cd8_status")+
  scale_color_manual(values=c("#41ab5d", "#6a51a3", "black"))+
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.line = element_blank())
p3 <- DimPlot(seu.szabo, reduction="umap", group.by="tissue")+
  scale_color_manual(values=c("#6a51a3", "#ef3b2c", "#41ab5d", "#4292c6"))+
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.line = element_blank())
p4 <- DimPlot(seu.szabo, reduction="umap", group.by="tissue_source")+
  scale_color_manual(values=c("#6a51a3", "#dadaeb", "#ef3b2c", "#fcbba1", "#41ab5d", "#c7e9c0", "#4292c6", "#c6dbef"))+
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.line = element_blank())
plot_grid(p1, p2, p3, p4, ncol=2)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/szabo_umapszabo_clusters.jpeg", width=10, height=9)


# Plot on VlnPlot
VlnPlot(seur.szabo.geps, group.by="tissue_source", features=gep_pbmc,
        cols = c("#6a51a3", "#dadaeb", "#ef3b2c", "#fcbba1", "#41ab5d", "#c7e9c0", "#4292c6", "#c6dbef"), ncol = 4)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/szabo_vlnplotszabo_scores.jpeg", width=15, height=10)

VlnPlot(seur.szabo.geps, group.by="tissue_source", split.by="cd4cd8_status", features="GEP_1",
        cols = c("#6a51a3", "#dadaeb", "#ef3b2c", "#fcbba1", "#41ab5d", "#c7e9c0", "#4292c6", "#c6dbef"))
VlnPlot(seur.szabo.geps, group.by="tissue_source", split.by="cd4cd8_status", features="GEP_4",
        cols = c("#6a51a3", "#dadaeb", "#ef3b2c", "#fcbba1", "#41ab5d", "#c7e9c0", "#4292c6", "#c6dbef"))

dfplot <- seur.szabo.geps@meta.data[,c("cd4cd8_status", "tissue_source", "tissue", "stimulation_status", gep_pbmc)] %>%
  rownames_to_column("cellid") %>%
  as_tibble() %>%
  pivot_longer(cols=starts_with("GEP"), names_to="GEP", values_to="score")

ggplot(dfplot %>% filter(cd4cd8_status != "unassigned" & GEP %in% c("GEP_1", "GEP_4", "GEP_6")),
       aes(x=tissue_source, y=score))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(fill=tissue_source))+
  geom_jitter(width=0.01, size=0.1)+
  facet_grid(GEP ~ cd4cd8_status)+
  scale_fill_manual(values=c("#6a51a3", "#dadaeb", "#ef3b2c", "#fcbba1", "#41ab5d", "#c7e9c0", "#4292c6", "#c6dbef"))+
  theme(axis.text.x=element_text(angle=45, hjust=1))
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/szabo_vlnplotszabo_scores2.jpeg", width=8, height=8)
