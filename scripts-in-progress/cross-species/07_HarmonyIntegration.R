# **************
# 1. IMPORT ####
# **************

# Import librairies
library(Seurat)
library(ggplot2)
library(tidyverse)
library(harmony)
library(RColorBrewer)

# Import data
seur.ms <- readRDS("~/Projects/HumanThymusProject/data/cross-species/00_Reproduce_UMAPs/ms_seurobj.rds")
DimPlot(seur.ms)
cols_nkt <- c("#810f7c", "#8856a7", "#8c96c6", "#b3cde3", "#edf8fb")
names(cols_nkt) <- unique(seur.ms$cell_type)
SCpubr::do_DimPlot(seur.ms, 
                   group.by = "cell_type",
                   label=T,
                   label.color="black",
                   legend.position = "none",
                   repel=T,
                   plot.title = "Mouse",
                   colors.use=cols_nkt,
                   font.size = 24)

# seur.hu <- readRDS("~/Projects/HumanThymusProject/data/human-thymus/HumanData_12_AnalysisByLineage/seurat_thymus-NKT_2023-02-27.rds")
seur.hu2 <- readRDS("~/Projects/HumanThymusProject/data/human-thymus/HumanData_12_AnalysisByLineage/thymus.nkt_02_28_23.RDS")
Idents(seur.hu) <- seur.hu$cell_annot
DimPlot(seur.hu)
cols_nkt2 <- c("#fd8d3c", "#feedde", "#fdbe85", "#a63603", "#e6550d")
names(cols_nkt2) <- unique(seur.hu$cell_annot)
colors_clusters_NKT <- c("hu0" = "#d8443c", "hu1" = "#e09351", "hu2" = "gold", "hu3" = "#74c8c3", "hu4" = "#5a97c1",
                         "hu5" = "#a40000", "hu6" = "#72bcd5")
SCpubr::do_DimPlot(seur.hu, 
                   group.by = "cell_annot",
                   label=T,
                   label.color="black",
                   legend.position = "none",
                   repel=T,
                   plot.title = "Human",
                   colors.use=cols_nkt2,
                   font.size = 24)

ortholog.df <- read.csv("~/Projects/HumanThymusProject/data/cross-species/03_BiomartTable/big_ass_ortholog_table.csv")
length(unique(ortholog.df$ms_symbol_data)) # 17,085 ms genes
length(unique(ortholog.df$hu_symbol)) # 17,152 hu genes




# ****************************
# 2. MERGE SEURAT OBJECTS ####
# ****************************

# Get counts, HVGs and metadata
# ms.hvg <- VariableFeatures(FindVariableFeatures(seur.ms))
ms.counts <- seur.ms[["RNA"]]@counts
ms.metadata <- seur.ms@meta.data

# hu.hvg <- VariableFeatures(seur.hu)
hu.counts <- seur.hu[["RNA"]]@counts
hu.metadata <- seur.hu@meta.data


# First, check whether genes can all be found in the ortholog table (many can't be found because I removed genes with orthology confidence=0)
table(unique(rownames(ms.counts)) %in% ortholog.df$ms_symbol_data) # 13,969 not
table(unique(rownames(hu.counts)) %in% ortholog.df$hu_symbol) # 4747 not


# Subset the ortholog table to only genes that we can "translate"
dictionary <- ortholog.df %>%
  as_tibble() %>%
  ## Intersection
  filter(ms_symbol_data %in% unique(rownames(ms.counts)) & hu_symbol %in% unique(rownames(hu.counts))) %>%
  ## Union
  # filter(ms_symbol_data %in% unique(rownames(ms.counts)) | hu_symbol %in% unique(rownames(hu.counts))) %>%
  # filter(hu_orthology_confidence == 1) %>%
  ## Remove any symbols that are NAs
  filter(!is.na(ms_symbol_data)) %>%
  filter(!is.na(hu_symbol)) %>%
  ## Keep only 1:1 orthologs
  group_by(ms_symbol_data) %>% filter(n_distinct(hu_symbol) == 1) %>% ungroup() %>%
  group_by(hu_symbol) %>% filter(n_distinct(ms_symbol_data) == 1) %>% ungroup()


# Translate the mouse HVGs into "human gene" language
# ms.hvg.translated <- pull(dictionary %>% filter(ms_symbol_data %in% ms.hvg), # not all ms HVGs are found in ortholog.df
#                           hu_symbol)
# hu.hvg.translated <- pull(dictionary %>% filter(hu_symbol %in% hu.hvg), # not all hu HVGs are found in ortholog.df
#                           hu_symbol)
# total.hvg <- unique(union(ms.hvg.translated, hu.hvg.translated))
# total.hvg <- total.hvg[!total.hvg %in% test$hu_symbol] # remove duplicates


# Keep only ms and hu genes that have 1:1 orthologs
table(unique(rownames(ms.counts)) %in% dictionary$ms_symbol_data) # 12,201 genes should have a translation
table(unique(rownames(hu.counts)) %in% dictionary$hu_symbol) # 12,201 genes should have a translation
ms.counts <- ms.counts[rownames(ms.counts) %in% dictionary$ms_symbol_data,]
hu.counts <- hu.counts[rownames(hu.counts) %in% dictionary$hu_symbol,]

# Translate the mouse genes in count table into "human gene"
ms.dict <- dictionary %>%
  filter(ms_symbol_data %in% rownames(ms.counts)) %>%
  select(ms_symbol_data, hu_symbol, hu_orthology_confidence) %>%
  # distinct() %>%
  # group_by(ms_symbol_data) %>% filter(n_distinct(hu_symbol)>1)
  distinct(ms_symbol_data, .keep_all=T)
ms.dict <- ms.dict[match(rownames(ms.counts), ms.dict$ms_symbol_data),]
table(ms.dict$ms_symbol_data == rownames(ms.counts)) # all true
table(is.na(ms.dict$hu_symbol), useNA="ifany") # no NAs
# Translate
rownames(ms.counts) <- ms.dict$hu_symbol


# Check if some key genes are included
table(c("ZBTB16", "EGR2", "TBX21", "RORC", "GATA3") %in% rownames(ms.counts))


# Merge everything into one
ms.metadata$study <- "Mouse"
ms.metadata$Method <- "10X"
ms.metadata <- ms.metadata[,c("cell_type", "study", "orig.ident", "Method")]
colnames(ms.metadata)[1] <- "cell_annot"
colnames(ms.metadata)[3] <- "Batch"
head(ms.metadata)

hu.metadata$study <- "Human"
hu.metadata <- hu.metadata[,c("cell_annot", "new_clusters_NKT", "study", "Batch", "Method")]
head(hu.metadata)

ms.seur <- CreateSeuratObject(counts=ms.counts, meta.data=ms.metadata, project="mouseNKT")
hu.seur <- CreateSeuratObject(counts=hu.counts, meta.data=hu.metadata, project="humanNKT")
seur.total <- merge(ms.seur, hu.seur)

# Sanity checks
head(seur.total@meta.data)
table(seur.total$cell_annot, useNA="ifany")
table(seur.total$study, useNA="ifany")
table(seur.total$Batch, useNA="ifany")
table(seur.total$Method, useNA="ifany")




# ***************************************
# 2. INTEGRATE DATASETS WITH HARMONY ####
# ***************************************

# Do preprocessing (without integration)
seur.int <- NormalizeData(seur.total)
seur.int <- FindVariableFeatures(seur.int, selection.method="vst", nfeatures=5000)
seur.int <- ScaleData(seur.int)
seur.int <- RunPCA(seur.int, npcs=50, verbose=TRUE, approx=FALSE)
ElbowPlot(seur.int, ndims=50, reduction="pca")
seur.int <- RunUMAP(seur.int, reduction = "pca", dims = 1:10, n.neighbors=50, reduction.key = "UMAP50_", reduction.name="umap_bef_int")
# seur.int <- FindNeighbors(seur.int, k.param=20, reduction = "pca", dims = 1:10)
# seur.int <- FindClusters(seur.int, resolution = 0.8, random.seed = 0)
# Visualization
DimPlot(seur.int, reduction = "umap_bef_int", group.by="study", label = TRUE, repel = TRUE)
DimPlot(seur.int, reduction = "umap_bef_int", group.by="cell_annot", label = TRUE, repel = TRUE)


# Do preprocessing (with integration)
seur.int <- RunHarmony(seur.int, reduction = "pca", max.iter.harmony=30,
                       group.by.vars = c("study", "Method", "Batch"),
                       lambda=c(0.5,0.5,0.5))
ElbowPlot(seur.int, ndims=50, reduction="harmony")
DimPlot(seur.int, reduction="pca", group.by="study") | DimPlot(seur.int, reduction="harmony", group.by="study")
DimPlot(seur.int, reduction="pca", group.by="Batch") + theme(legend.position="none") | DimPlot(seur.int, reduction="harmony", group.by="Batch") + theme(legend.position="none")
DimPlot(seur.int, reduction="pca", group.by="Method") | DimPlot(seur.int, reduction="harmony", group.by="Method")

seur.int <- RunUMAP(seur.int, reduction = "harmony", dims = 1:10, n.neighbors=50, reduction.key = "UMAP50_")
# seur.int <- FindNeighbors(seur.int, k.param=20, reduction = "harmony", dims = 1:10)
# seur.int <- FindClusters(seur.int, resolution = 0.8, random.seed = 0)

# saveRDS(seur.int, "~/Projects/HumanThymusProject/data/human-thymus/HumanData_16_ParkIntegration/seurat_integrated_park_23-03-13.rds") # checkpoint save
# seur.int <- readRDS("~/Projects/HumanThymusProject/data/human-thymus/HumanData_16_ParkIntegration/seurat_integrated_park_23-03-13.rds")


# Add laurent's clusters
# table(rownames(seur.int@meta.data[seur.int@meta.data$study=="Human",])==rownames(seur.hu2@meta.data), useNA="ifany")
# seur.int@meta.data$clusters_NKT <- seur.int@meta.data$cell_annot
# seur.hu2@meta.data$new_clusters_NKT <- paste0("hu", seur.hu2@meta.data$new_clusters_NKT)
# seur.int@meta.data[seur.int@meta.data$study=="Human","clusters_NKT"] <- seur.hu2@meta.data$new_clusters_NKT
# table(seur.int$cell_annot, useNA="ifany")
# table(seur.int$clusters_NKT, useNA="ifany")


# Visualization sanity check
DimPlot(seur.int, reduction = "umap", group.by="Batch", label = T, shuffle=T)
DimPlot(seur.int, reduction = "umap", group.by="Method", label = T, shuffle=T)

DimPlot(seur.int, reduction = "umap", group.by="study", label = F, shuffle=T)+
  scale_color_manual(values=c("#d6604d", "#4393c3"))+
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_text(hjust=0),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) |
DimPlot(seur.int, reduction = "umap", group.by="clusters_NKT", label = TRUE, repel = TRUE)+
  # scale_color_manual(values=c(cols_nkt, cols_nkt2))+
  scale_color_manual(values=c(cols_nkt, colors_clusters_NKT))+
  labs(title="Cell annotation")+
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
# ggsave("~/Projects/HumanThymusProject/data/cross-species/07_HarmonyIntegration/umap_harmony_integration3.jpeg", width=12, height=5)




# **************************************
# 2. INTEGRATE DATASETS WITH SEURAT ####
# **************************************

# normalize and identify variable features for each dataset independently
seur.total.list <- SplitObject(seur.total, split.by = "study")
seur.total.list <- lapply(X = seur.total.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seur.total.list)

# perform integration
anchors <- FindIntegrationAnchors(object.list = seur.total.list, anchor.features = features)
seur.combined <- IntegrateData(anchorset = anchors)
DefaultAssay(seur.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
seur.combined <- ScaleData(seur.combined, verbose = FALSE)
seur.combined <- RunPCA(seur.combined, npcs = 50, verbose = FALSE, approx=FALSE)
ElbowPlot(seur.int, ndims=50, reduction="pca")
seur.combined <- RunUMAP(seur.combined, reduction = "pca", dims = 1:10)
# seur.combined <- FindNeighbors(seur.combined, reduction = "pca", dims = 1:30)
# seur.combined <- FindClusters(seur.combined, resolution = 0.5)

# Add laurent's clusters
# table(rownames(seur.combined@meta.data[seur.combined@meta.data$study=="Human",])==rownames(seur.hu2@meta.data), useNA="ifany")
# seur.combined@meta.data$clusters_NKT <- seur.combined@meta.data$cell_annot
# seur.hu2@meta.data$new_clusters_NKT <- paste0("hu", seur.hu2@meta.data$new_clusters_NKT)
# seur.combined@meta.data[seur.combined@meta.data$study=="Human","clusters_NKT"] <- seur.hu2@meta.data$new_clusters_NKT
# table(seur.combined$cell_annot, useNA="ifany")
# table(seur.combined$clusters_NKT, useNA="ifany")


# Visualize
DimPlot(seur.combined, reduction = "umap", group.by="Batch", label = T, shuffle=T)
DimPlot(seur.combined, reduction = "umap", group.by="Method", label = T, shuffle=T)

DimPlot(seur.combined, reduction = "umap", group.by="study", label = F, shuffle=T)+
  scale_color_manual(values=c("#d6604d", "#4393c3"))+
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_text(hjust=0),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) |
DimPlot(seur.combined, reduction = "umap", group.by="clusters_NKT", label = TRUE, repel = TRUE)+
  # scale_color_manual(values=c(cols_nkt, cols_nkt2))+
  scale_color_manual(values=c(cols_nkt, colors_clusters_NKT))+
  labs(title="Cell annotation")+
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
# ggsave("~/Projects/HumanThymusProject/data/cross-species/07_HarmonyIntegration/umap_seurat_integration3.jpeg", width=12, height=5)

DimPlot(seur.combined, reduction = "umap", group.by="clusters_NKT", split.by="study", label = TRUE, repel = TRUE)+
  # scale_color_manual(values=c(cols_nkt, cols_nkt2))+
  scale_color_manual(values=c(cols_nkt, colors_clusters_NKT))+
  labs(title="Cell annotation")+
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_text(hjust=0),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
# ggsave("~/Projects/HumanThymusProject/data/cross-species/07_HarmonyIntegration/umap_seurat_integration4.jpeg", width=8, height=5)
