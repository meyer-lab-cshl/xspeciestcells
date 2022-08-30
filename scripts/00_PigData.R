###
# Purpose: Import & plot pig iNKT data (curtesy of ...)
# Date: Aug 9th 2022
# Author: Salomé Carcy
###


#### IMPORT ####

# Libraries
library(ggplot2)
library(Seurat)

# Data
path.plots <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/00_Reproduce_UMAPs"
path.data <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/pig_data/GSE192520_RAW/iNKT"
pig <- Read10X(path.data)

# Look at individual files
# # cellular barcodes present in dataset
# readr::read_tsv(file.path(path.data, "barcodes.tsv.gz"), col_names = FALSE) # 12,097 cells
# # IDs of quantified genes
# readr::read_tsv(file.path(path.data, "features.tsv.gz"), col_names = FALSE) # 21,301 genes




#### PRE-PROCESSING ####

# Create Seurat Object
seur.pig <- CreateSeuratObject(counts = pig, project = "pig_Thymus_NKT", min.cells = 3, min.features = 200)
seur.pig # 12,097 cells


# Identify mitochondrial genes
# Use AnnotationHub to retrieve chromosomal locations given the Ensembl IDs (species: susscrofa)
# library(AnnotationHub)
# ah <- AnnotationHub()
# display(ah) # find the reference for "the pig reference genome Sscrofa 11.1 with gene transfer format file 11.1.98"
# ens.ss.v98 <- AnnotationHub()[["AH75101"]] # Sscrofa11.1 (v98)
# # Get the list of genes (ENSEMBL ID)
# features.df <- readr::read_tsv(file.path(path.data, "features.tsv.gz"), col_names = FALSE)
# # Get the list of genes & their chromosomal location
# chr.loc <- mapIds(ens.ss.v98, keys=features.df$X1,
#                   keytype="GENEID", column="SEQNAME")
# # Keep the mitochondrial genes
# is.mito <- features.df[features.df$X1 %in% names(chr.loc[chr.loc=="MT"]),]
# print(is.mito) # the genes correspond to the one that the authors identified as mitochondrial (see below PercentageFeatureSet)
# # check the chromosomal location of genes that start with "MT"
# table(chr.loc[features.df$X1[stringr::str_detect(test$X2, "^MT")]])


# QC
seur.pig[["percent.mt"]] <- PercentageFeatureSet(seur.pig, features=c("ND1","ND2", "ND3", "ND4",
                                                                      "ND4L", "ND5","ND6","COX1", "COX2", "COX3", "CYTB", "ATP6", "ATP8"))
head(seur.pig@meta.data, 5)
FeatureScatter(seur.pig, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seur.pig, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(seur.pig, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=.01)
# ggsave(file.path(path.plots, "pig_qc.jpeg"), width=10, height=6)
seur.pig <- subset(seur.pig, subset = nFeature_RNA >= 800 & nFeature_RNA < 5000 & percent.mt < 11)
seur.pig # 11,123 cells


# Normalize and find variable features
seur.pig <- NormalizeData(seur.pig)
seur.pig <- FindVariableFeatures(seur.pig, selection.method = "vst", nfeatures = 2000) # default parameters
# plot variable features with labels
LabelPoints(plot = VariableFeaturePlot(seur.pig),
            points = head(VariableFeatures(seur.pig), 10), repel = TRUE)
# ggsave(file.path(path.plots, "pig_hvg.jpeg"), width=8, height=6, bg="white")
seur.pig <- ScaleData(seur.pig, features = rownames(seur.pig))


# Regress out cell cyle 
s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seur.pig  <- CellCycleScoring(seur.pig, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# head(seur.pig@meta.data, 5) # cell cycle scoring is incorporated in the metadata
seur.pig <- ScaleData(seur.pig, vars.to.regress = c("S.Score","G2M.Score"), features = rownames(seur.pig)) # take 20-30min to run


# Run dimension reduction
seur.pig <- RunPCA(seur.pig, npcs = 30, verbose = FALSE)
# ElbowPlot(seur.pig)
seur.pig <- FindNeighbors(seur.pig, dims = 1:19)
set.seed(123)
seur.pig <- FindClusters(seur.pig, resolution = 0.5) # random seed!!
seur.pig <- RunUMAP(seur.pig, dims = 1:19)

# Cluster numbers are not in the same order as in the paper, so we'll just replace them
DimPlot(seur.pig, pt.size = 0.1, label = T)
# new.cluster.ids <- c("c6", #0
#                      "c4", #1
#                      "c3-5", #2
#                      "c1", #3
#                      "c1-4",#4
#                      "c7", #5
#                      "c8", #6
#                      "c2", #7
#                      "c8") #8
new.cluster.ids <- c("iNKT2.3",  #0
                     "iNKT2.1",  #1
                     "iNKTt", #2
                     "iNKT2.0",  #3
                     "iNKT2.2",  #4
                     "iNKT-ISG",    #5
                     "iNKTc",    #6
                     "iNKT-sw")  #7
names(new.cluster.ids) <- levels(seur.pig)
seur.pig <- RenameIdents(seur.pig, new.cluster.ids)
# levels(seur.pig) <- c("c1", "c2", "c3-5", "c4", "c1-4", "c6", "c7", "c8")
levels(seur.pig) <- c("iNKT2.0", "iNKTc", "iNKTt", "iNKT2.1", "iNKT2.2", "iNKT2.3", "iNKT-ISG", "iNKT-sw")




#### REPRODUCE FIGURES ####

# Display UMAP
DimPlot(seur.pig, reduction = "umap", label = T)
# ggsave(file.path(path.plots, "pig_umap.jpeg"), width=7, height=5)


# Dotplot
genes<-c("XCL1", "STAT4", "CXCR3", "NKG7", "TBX21",
         "IFIT1", "MX1", "ISG15",
         "S1PR1", "CD9", "S100A6", "S100A4", "KLF2",
         "LEF1", "CCR9", "SATB1", "TOX2",
         "MKI67", "HMGB2", "PCNA", "PCLAF", # MKI67 ensembl id not provided in the supplementaries!
         "ID3", "EGR2", "CD69", "ENSSSCG00000033465",
         "CD247", "EIF3I", "CD200", "IL4", "IL6R", "ICOS", "FOXO1", "ENSSSCG00000009240", "CCR7", "ZBTB16", "GATA3")
DotPlot(seur.pig, features = rev(genes)) + coord_flip()+theme(text = element_text(face="bold"))
ggsave(file.path(path.plots, "pig_dotplot.jpeg"), width=7, height=7, bg="white")




#### DE Genes & Average Exp data for cross-species ####

# Get differentially expressed genes by cluster
diff.markers <- FindAllMarkers(seur.pig, assay="RNA", only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.25)
dim(diff.markers) # 3457 genes
# sanity check
# top5 <- diff.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# DoHeatmap(seur.pig, features = top5$gene, slot="data") + NoLegend()
saveRDS(diff.markers, "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/02_CorrelationComparion/pig_DEG.rds")


# Get average expression data per cluster
avgexp <- AverageExpression(seur.pig)
avgexp <- as.data.frame(avgexp$RNA) # get it into a df
head(avgexp, 10)
dim(avgexp) # 8 clusters and 13,922 genes
saveRDS(avgexp, "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/02_CorrelationComparion/pig_avgexp.rds")





#### TRY A DIFFERENT PROCESSING (personal) ####

seur.pig <- CreateSeuratObject(counts = pig, project = "pig_Thymus_NKT", min.cells = 3, min.features = 200)
seur.pig # 12,097 cells

# QC
seur.pig[["percent.mt"]] <- PercentageFeatureSet(seur.pig, features=c("ND1","ND2", "ND3", "ND4",
                                                                      "ND4L", "ND5","ND6","COX1", "COX2", "COX3", "CYTB", "ATP6", "ATP8"))
head(seur.pig@meta.data, 5)
seur.pig <- subset(seur.pig, subset = nFeature_RNA > 800 & nFeature_RNA < 5000 & nCount_RNA >=1000 & percent.mt < 8)
seur.pig # 10,759 cells


# Normalize and find variable features
set.seed(123)
seur.pig <- NormalizeData(seur.pig) # normalized data is saved in seur.combined[["RNA"]]@data
seur.pig <- FindVariableFeatures(seur.pig) # return 2,000 HVGs
seur.pig <- ScaleData(seur.pig, features = rownames(seur.pig))

# Run dimension reduction
seur.pig <- RunPCA(seur.pig, features = VariableFeatures(object = seur.pig))
ElbowPlot(seur.pig)
seur.pig <- FindNeighbors(seur.pig, dims = 1:15)
seur.pig <- FindClusters(seur.pig, resolution = 0.8)
seur.pig <- RunUMAP(seur.pig, dims = 1:15)

DimPlot(seur.pig, pt.size = 0.1, label = T)
