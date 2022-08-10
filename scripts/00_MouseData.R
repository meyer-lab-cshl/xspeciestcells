###
# Purpose: Import & plot mouse iNKT data (curtesy of Laurent Gapin)
# Date: Aug 9th 2022
# Author: Salomé Carcy
###


#### IMPORT ####

# Libraries
library(ggplot2)
library(Seurat)
# remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
library(cowplot)
library(tidyverse)

# Data
path.plots <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/00_Reproduce_UMAPs"
path.data <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/mouse_data"
mouse1 <- Read10X(file.path(path.data, "B6_1"))
mouse2 <- Read10X(file.path(path.data, "B6_2"))

# Look at individual files
# # cellular barcodes present in dataset
# readr::read_tsv(file.path(path.data, "B6_1/barcodes.tsv.gz"), col_names = FALSE) # 8,373 cells mouse 1
# readr::read_tsv(file.path(path.data, "B6_2/barcodes.tsv.gz"), col_names = FALSE) # 4,243 cells mouse 2
# # IDs of quantified genes
# readr::read_tsv(file.path(path.data, "B6_1/features.tsv.gz"), col_names = FALSE) # 31,053 genes
# readr::read_tsv(file.path(path.data, "B6_2/features.tsv.gz"), col_names = FALSE) # 31,053 genes




#### PRE-PROCESSING ####

# Create Seurat Object
seur.ms1 <- CreateSeuratObject(mouse1, project="B6_Thymus_NKT_1")
seur.ms2 <- CreateSeuratObject(mouse2, project="B6_Thymus_NKT_2")

# Combine into one seurat object
seur.combined <- merge(seur.ms1, y = seur.ms2, add.cell.ids = c('B61', 'B62'), project = "NKT")


# QC
seur.combined[["percent.mt"]] <- PercentageFeatureSet(seur.combined, pattern = "^mt-")
head(seur.combined@meta.data, 5)
FeatureScatter(seur.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seur.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(seur.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=.01)
# ggsave(file.path(path.plots, "ms_qc.jpeg"), width=10, height=6)
seur.combined <- subset(seur.combined, subset = nFeature_RNA >= 500 & nFeature_RNA < 4100 & nCount_RNA >=1000 & percent.mt < 6)
table(seur.combined@meta.data$orig.ident) # 7,231 cells mouse 1 and 3,651 cells mouse 2


# Normalize and find variable features
set.seed(123)
seur.combined <- NormalizeData(seur.combined) # normalized data is saved in seur.combined[["RNA"]]@data
seur.combined <- FindVariableFeatures(seur.combined) # return 2,000 HVGs
# plot variable features with labels
LabelPoints(plot = VariableFeaturePlot(seur.combined),
            points = head(VariableFeatures(seur.combined), 10), repel = TRUE)
# ggsave(file.path(path.plots, "ms_hvg.jpeg"), width=8, height=6)


# Run integration & dimension reduction
seur.combined <- RunFastMNN(object.list = SplitObject(seur.combined, split.by = "orig.ident"), k = 20)
seur.combined <- RunUMAP(seur.combined, reduction = "mnn", dims = 1:30)
seur.combined <- FindNeighbors(seur.combined, reduction = "mnn", dims = 1:30, compute.SNN = TRUE)
seur.combined <- FindClusters(seur.combined, resolution = 0.8)

# Cluster numbers are not in the same order as in the paper, so we'll just replace them
DimPlot(seur.combined, group.by = c("seurat_clusters"), pt.size = 0.1, label = T)
new.cluster.ids <- c("c8", #0
                     "c4", #1
                     "c6", #2
                     "c1", #3
                     "c10",#4
                     "c9", #5
                     "c0", #6
                     "c7", #7
                     "c3", #8
                     "c5", #9
                     "c2", #10
                     "c11") # 11
names(new.cluster.ids) <- levels(seur.combined)
seur.combined <- RenameIdents(seur.combined, new.cluster.ids)
levels(seur.combined) <- paste0("c", 0:11)




#### REPRODUCE FIGURES ####

# Display UMAP
p1 <- DimPlot(seur.combined, group.by = c("orig.ident"), pt.size = 0.1, label = F) + theme(legend.position="top", title = element_text(size=0))
p2 <- DimPlot(seur.combined, pt.size = 0.1, label = TRUE) + labs(title="Integrated (10882 cells)") + theme(legend.position="none")
p3 <- DimPlot(seur.combined, pt.size = 0.1, split.by = "orig.ident", ncol = 4)

jpeg(filename = file.path(path.plots, "ms_umap.jpeg"), width=4500, height=1500, res=300)
ggdraw()+
  draw_plot(p1, x=0, y=0, width=.25, height=1)+
  draw_plot(p2, x=.25, y=0, width=.25, height=1)+
  draw_plot(p3, x=.5, y=0, width=.5, height=1)
dev.off()
# ggsave(file.path(path.plots, "ms_umap.jpeg"), width=15, height=6) # weird error message


# Display differentially expressed genes on UMAP
FeaturePlot(seur.combined, features = c('Cd24a', 'Mki67', 'Zbtb16', 'Rorc', 'Tbx21'),
            cols = c('#deebf7', 'red'), pt.size = 0.7, ncol = 5)
# ggsave(file.path(path.plots, "ms_umap_deg.jpeg"), width=25, height=6)


# Create a dotplot of genes highly expressed in different iNKT subsets
markers.to.plot <- c("Cd24a", "Cd69", "Egr2", "Il2rb", "Ifng", "Tbx21", "Ccr6", "Rorc", "Il17rb", 
                     "Ccr7", "Il4", "Gata3", "Zbtb16", "Cd4", "Cd44")
DotPlot(seur.combined, features = markers.to.plot, cols = "RdBu", dot.scale = 8, col.min=-1, col.max=1) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(path.plots, "ms_dotplot.jpeg"), width=8, height=6, bg="white")


#Create a heatmap of differentially expressed genes by cluster
diff.markers <- FindAllMarkers(seur.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- diff.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# DoHeatmap(seur.combined, features = top5$gene, slot="data") + NoLegend()




#### WHAT IS CLUSTER 11? ####

VlnPlot(seur.combined, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), pt.size = 0.5, ncol=1)
ggsave(file.path(path.plots, "ms_c11_vlnplot.jpeg"), width=8, height=10, bg="white")
