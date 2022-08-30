###
# Purpose: Import & plot human iNKT data (curtesy of Laurent Gapin)
# Date: Aug 9th 2022
# Author: Salomé Carcy
###


#### IMPORT ####

# Libraries
library(ggplot2)
library(cowplot)
library(Seurat)
library(SeuratWrappers)
library(Matrix)
library(tidyverse)

# Data (sorted human iNKT cells from 3 subjects)
path.plots <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/00_Reproduce_UMAPs"
path.data <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/human_data"
human5 <- read.csv(file.path(path.data, "CUThy13_220225_SampleTag05_hs_NKT_RSEC_MolsPerCell.csv"), sep=",", header=T, skip=7) # 404 cells
human8 <- read.csv(file.path(path.data, "CUTHY11BDRscRNA_seq_091621_SampleTag08_hs_NKT_RSEC_MolsPerCell.csv"), sep=",", header=T, skip=8) # 1913 cells
human12 <- read.csv(file.path(path.data, "CUTHY12BDRscRNA_seq_211101_SampleTag12_hs_NKT_RSEC_MolsPerCell.csv"), sep=",", header=T, skip=8) # 344 cells

# Any cluster ID in the colnames?
# table(stringr::str_detect(colnames(human5),"[[:lower:]]"))
# colnames(human5)[stringr::str_detect(colnames(human5),"[[:lower:]]")]


# Current df format:
# Rownames || Cell_Index | Gene1 | Gene2 | ...
#   -      ||   421953   |   0   |   0   | ...
#   -      ||   459121   |   0   |   3   | ...

# Invert the dataframes (cells as columns and genes as rows) and transform to dgCMatrix
convert_to_matrix <- function(df){
  rownames(df) <- df$Cell_Index
  df$Cell_Index <- NULL
  df <- t(df)
  df <- Matrix(df, sparse=T)
  return(df)
}

mat.hu5  <- convert_to_matrix(human5) # 404 cells
mat.hu8  <- convert_to_matrix(human8) # 1913 cells
mat.hu12 <- convert_to_matrix(human12) # 344 cells




#### PRE-PROCESSING ####

# Create Seurat Object
seur.h5  <- CreateSeuratObject(mat.hu5, project="Hu_Thymus_NKT_5")
seur.h8  <- CreateSeuratObject(mat.hu8, project="Hu_Thymus_NKT_8")
seur.h12 <- CreateSeuratObject(mat.hu12, project="Hu_Thymus_NKT_12")

# Combine into one seurat object
seur.combined <- merge(seur.h5, y=c(seur.h8,seur.h12), add.cell.ids=c('Hu5', 'Hu8', 'Hu12'), project="HuNKT") # 2661 cells


# Identify mitochondrial genes
# Use AnnotationHub to retrieve chromosomal locations given the Ensembl IDs (species: susscrofa)
# library(AnnotationHub)
# ah <- AnnotationHub()
# display(ah) # find the reference for "the pig reference genome Sscrofa 11.1 with gene transfer format file 11.1.98"
# ens.hu.v29 <- AnnotationHub()[["AH75174"]] # GRCh38 (Gencode v29 annotated genes)
# # Can we find our human genes in the GRCh38 reference?
# table(rownames(seur.combined) %in% ens.hu.v29$Gene) # not all
# # Let's get the list of mitochondrial genes from the GRCh38 reference
# hu.mt <- subset(ens.hu.v29, seqnames(ens.hu.v29)=="chrM")
# mt.genes <- as.character(hu.mt$Gene)
# mt.genes <- mt.genes[!is.na(mt.genes)]
# print(mt.genes)
# # Compare the reference list of mitochondrial genes with the human genes from seur.combined that start with MT.
# mt.genes <- paste0("MT.", mt.genes)
# sort(mt.genes)
# sort(rownames(seur.combined)[stringr::str_detect(rownames(seur.combined), "^MT\\.")])
# # same genes!! => genes starting with MT. in seur.combined are mitochondrial


# QC
seur.combined[["percent.mt"]] <- PercentageFeatureSet(seur.combined, pattern = "^MT\\.")
head(seur.combined@meta.data, 5)
FeatureScatter(seur.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seur.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(seur.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=.01)
# ggsave(file.path(path.plots, "hu_qc.jpeg"), width=10, height=6)
seur.combined <- subset(seur.combined, subset = nFeature_RNA >= 500 & nFeature_RNA < 4000 & nCount_RNA >=500 & percent.mt < 25)
table(seur.combined@meta.data$orig.ident) # 399 cells hu5; 1828 cells hu8; 331 cells hu12


# Normalize and find variable features
set.seed(123)
seur.combined <- NormalizeData(seur.combined) # normalized data is saved in seur.combined[["RNA"]]@data
seur.combined <- FindVariableFeatures(seur.combined) # return 2,000 HVGs
# plot variable features with labels
LabelPoints(plot = VariableFeaturePlot(seur.combined),
            points = head(VariableFeatures(seur.combined), 10), repel = TRUE)
# ggsave(file.path(path.plots, "hu_hvg.jpeg"), width=8, height=6, bg="white")


# Run integration & dimension reduction
seur.combined <- RunFastMNN(object.list = SplitObject(seur.combined, split.by = "orig.ident"), k = 20)
seur.combined <- RunUMAP(seur.combined, reduction = "mnn", dims = 1:30)
seur.combined <- FindNeighbors(seur.combined, reduction = "mnn", dims = 1:30, compute.SNN = TRUE)
seur.combined <- FindClusters(seur.combined, resolution = 0.8)




#### FIGURES ####

# Display UMAP
p1 <- DimPlot(seur.combined, group.by = c("orig.ident"), pt.size = 0.1, label = F) + theme(legend.position="top", title = element_text(size=0))
p2 <- DimPlot(seur.combined, pt.size = 0.1, label = TRUE) + labs(title="Integrated (2558 cells)") + theme(legend.position="none")
p3 <- DimPlot(seur.combined, pt.size = 0.1, split.by = "orig.ident", ncol = 3)

jpeg(filename = file.path(path.plots, "hu_umap.jpeg"), width=5000, height=1500, res=300)
ggdraw()+
  draw_plot(p1, x=0, y=0, width=.25, height=1)+
  draw_plot(p2, x=.25, y=0, width=.25, height=1)+
  draw_plot(p3, x=.5, y=0, width=.5, height=1)
dev.off()


# Display differentially expressed genes on UMAP
rownames(seur.combined)[stringr::str_detect(rownames(seur.combined),"CD4")]
FeaturePlot(seur.combined, features = c('CD24', 'MKI67', 'ZBTB16', 'RORC', 'TBX21'),
            cols = c('#deebf7', 'red'), pt.size = 0.7, ncol = 5)
ggsave(file.path(path.plots, "hu_umap_deg.jpeg"), width=25, height=6)


# Create a dotplot of genes highly expressed in different iNKT subsets
levels(seur.combined) <- c("0", "4", "1", "5", "6", "3", "2")
markers.to.plot <- c("CD24", "CD69", "EGR2", "IL2RB", "IFNG", "TBX21", "CCR6", "RORC", "IL17RB", 
                     "CCR7", "IL4", "GATA3", "ZBTB16", "CD4", "CD44")
DotPlot(seur.combined, features = markers.to.plot, cols = "RdBu", dot.scale = 8, col.min=-1, col.max=1) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
# ggsave(file.path(path.plots, "hu_dotplot.jpeg"), width=8, height=6, bg="white")

# Create a heatmap of differentially expressed genes by cluster
levels(seur.combined) <- as.character(0:6)
diff.markers <- FindAllMarkers(seur.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- diff.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(seur.combined, features = top5$gene, slot="data") + scale_fill_gradientn(colors = c("#f0f0f0", "#b2182b"))
# ggsave(file.path(path.plots, "hu_heatmap.jpeg"), width=8, height=6, bg="white")
