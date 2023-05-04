###
# Purpose: Import & plot mouse MAIT data (Legoux paper)
# Date: Oct 31st 2022
# Author: Salomé Carcy
###


#### IMPORT ####

# Libraries
library(ggplot2)
library(Seurat)
library(SeuratWrappers)
library(cowplot)
library(tidyverse)

# Data
path.data <- "./data/raw_data/mouse_data/MAIT_Legoux/Processed_data"
mouse1 <- read.table(file.path(path.data, "processed_ikura_tobiko_wt01.txt"), header=T)
mouse2 <- read.table(file.path(path.data, "processed_ikura_tobiko_wt02.txt"), header=T)

mouse1[1:5,1:5]
mouse2[1:5,1:5]

# Put gene names as rownames
rownames(mouse1) <- mouse1$symbol
rownames(mouse2) <- mouse2$symbol
mouse1$symbol <- NULL
mouse2$symbol <- NULL
# sanity check
mouse1[1:5,1:5]
mouse2[1:5,1:5]

# Keep only common gene names
genes_common <- intersect(rownames(mouse1), rownames(mouse2))
length(genes_common) # 12,033 genes
mouse1 <- mouse1[genes_common,]
mouse2 <- mouse2[genes_common,]
dim(mouse1) # sanity check
dim(mouse2) # sanity check





#### PRE-PROCESSING ####

# Create Seurat Object
seur.ms1 <- CreateSeuratObject(counts=mouse1, project="B6_Thymus_MAIT_1", min.cells=3)
seur.ms2 <- CreateSeuratObject(counts=mouse2, project="B6_Thymus_MAIT_2", min.cells=3)

# Combine into one seurat object
seur.combined <- merge(seur.ms1, y = seur.ms2, add.cell.ids = c('B61', 'B62'), project = "MAIT")


# QC
seur.combined[["percent.mt"]] <- PercentageFeatureSet(seur.combined, pattern = "^mt-")
head(seur.combined@meta.data, 5)
FeatureScatter(seur.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seur.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(seur.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=.01)
# ggsave(file.path(path.plots, "ms_qc.jpeg"), width=10, height=6)
seur.combined <- subset(seur.combined, subset = nFeature_RNA >= 500 & percent.mt < 5)
table(seur.combined@meta.data$orig.ident) # 3,823 cells mouse 1 and 3,061 cells mouse 2


# Normalize and find variable features
set.seed(123)
seur.combined <- NormalizeData(seur.combined) # normalized data is saved in seur.combined[["RNA"]]@data
seur.combined <- FindVariableFeatures(seur.combined) # return 2,000 HVGs
# plot variable features with labels
# LabelPoints(plot = VariableFeaturePlot(seur.combined),
#             points = head(VariableFeatures(seur.combined), 10), repel = TRUE)


# Run integration & dimension reduction
seur.combined <- RunFastMNN(object.list = SplitObject(seur.combined, split.by = "orig.ident"), cos.norm=TRUE, k = 15)
seur.combined <- RunUMAP(seur.combined, reduction = "mnn", dims = 1:20, return.model=T, min.dist=0.3, spread=1)
seur.combined <- FindNeighbors(seur.combined, reduction = "mnn", dims = 1:10, k.param=15, compute.SNN = TRUE)
seur.combined <- FindClusters(seur.combined, resolution = 0.5)

# Cluster numbers are not in the same order as in the paper, so we'll just replace them
DimPlot(seur.combined, group.by = c("seurat_clusters"), pt.size = 0.1, label = T)
DimPlot(seur.combined, group.by = c("orig.ident"), pt.size = 0.1, label = T)




#### REPRODUCE FIGURES ####

# Figure 1B
FeaturePlot(seur.combined, features = c('Cd24a', 'Zbtb16', 'Cd44', 'Tbx21', 'Cxcr3', 'Itga1'),
            cols = c('#deebf7', 'red'), pt.size = 0.7, ncol = 3, order=T)
FeaturePlot(seur.combined, features = c('Rorc', 'Ccr6', 'Itga5', 'Itgb3', 'Hells', 'Cdk1', 'Mki67', 'Sell'),
            cols = c('#deebf7', 'red'), pt.size = 0.7, ncol = 4, order=T)


# Figure 1A (MAIT1 and MAIT17 signatures)
# Get the MAIT1 and MAIT17 signatures
mait_signatures <- read.csv("./data/cross-species/00_Reproduce_UMAPs/MAIT_Legoux_signatures.csv", header=T)
# Compute gene scores
seur.combined <- AddModuleScore(object = seur.combined, assay = "RNA",
                                       features = list("mait1"=mait_signatures$MicroArray_thyMAIT1[mait_signatures$MicroArray_thyMAIT1!=""],
                                                       "mait17"=mait_signatures$MicroArray_thyMAIT17),
                                       name=c("mait1", "mait17"))
FeaturePlot(seur.combined, features="mait11", cols = c('#deebf7', 'red'))+
  labs(title="MAIT1 microarray signature")
FeaturePlot(seur.combined, features="mait172", cols = c('#deebf7', 'red'))+
  labs(title="MAIT17 microarray signature")




#### PERFORM CELL ANNOTATION ####

seur.combined@meta.data$cell_type <- case_when(
  seur.combined@meta.data$seurat_clusters %in% c(0,2) ~ "MAIT17b",
  seur.combined@meta.data$seurat_clusters == 1 ~ "Cluster 7",
  seur.combined@meta.data$seurat_clusters == 3 ~ "CyclingS",
  seur.combined@meta.data$seurat_clusters == 4 ~ "MAIT1",
  seur.combined@meta.data$seurat_clusters == 5 ~ "MAIT17a",
  seur.combined@meta.data$seurat_clusters == 6 ~ "CyclingG2M",
  seur.combined@meta.data$seurat_clusters == 7 ~ "MAIT0",
)

cols <- c("MAIT0"     = "#fee391",
          "Cluster 7" = "#bf812d",
          "MAIT1"     = "#80cdc1",
          "MAIT17a"   = "#084594",
          "MAIT17b"   = "#6baed6",
          "CyclingS"  = "#f768a1",
          "CyclingG2M"= "#8c6bb1")

DimPlot(seur.combined, group.by = c("cell_type"), pt.size = 0.7, label = T)+
  scale_color_manual(values=cols)+
  labs(title="Legoux dataset")
# ggsave("./data/cross-species/00_Reproduce_UMAPs/ms_mait_umap.jpeg", width=7, height=6)

# Dotplot with markers
Idents(seur.combined) <- "cell_type"
legoux.markers <- FindAllMarkers(seur.combined, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25)
test <- legoux.markers %>%
  group_by(cluster) %>%
  slice_max(n=5, order_by=avg_log2FC) %>%
  mutate(cluster=factor(cluster, levels=c("MAIT17b", "MAIT17a", "MAIT1", "MAIT0", "CyclingS", "CyclingG2M", "Cluster 7"))) %>%
  arrange(cluster)

DotPlot(seur.combined,
        features=unique(test %>% pull(gene)),
        group.by="cell_type")+
  theme(axis.text.x=element_text(angle=45, hjust=1))
# ggsave("./data/cross-species/00_Reproduce_UMAPs/ms_mait_dotplot.jpeg", width=10, height=5)

# Save
seur.combined@meta.data$Batch <- case_when(
  seur.combined@meta.data$orig.ident == "B6_Thymus_MAIT_1" ~ "ms_wt_1",
  seur.combined@meta.data$orig.ident == "B6_Thymus_MAIT_2" ~ "ms_wt_2"
)
seur.combined@meta.data$RNA_snn_res.0.5 <- NULL
colnames(seur.combined@meta.data)[6:7] <- c("score_mait1", "score_mait17")
# saveRDS(seur.combined, "./data/cross-species/00_Reproduce_UMAPs/ms_mait_seurobj.rds")




#### DECONTAMINATION CHECK ####
library(SingleCellExperiment)
library(celda)
library(scales)

counts   <- seur.combined@assays$RNA@counts
metadata <- seur.combined@meta.data
umap <- seur.combined[["umap"]]@cell.embeddings
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata,
                            reducedDims = list(umap = umap)) # 6,884 cells
# scater::plotReducedDim(sce, dimred="umap", colour_by="cell_type", point_size=0.5) # sanity check

# Run decontX with new_clusters
sce <- decontX(sce, z=sce$cell_type, batch=sce$orig.ident)
# scater::plotDimReduceCluster(x = sce$decontX_clusters,
#     dim1 = reducedDim(sce, "umap")[, 1],
#     dim2 = reducedDim(sce, "umap")[, 2])
scater::plotReducedDim(sce, dimred="umap", colour_by="decontX_contamination", point_size=0.1)+
  # scale_colour_gradient2(low=muted("blue"), mid="white", high=muted("red"), midpoint=0.25)+
  labs(title="decontX", color="Contamination")




