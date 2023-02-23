library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)

# https://satijalab.org/seurat/articles/spatial_vignette.html#dataset-2
directory <- "data/GSE207205_RAW"
sample_ids <- c("GSM6281327", "GSM6281320", "GSM6281322",  "GSM6281321",
                "GSM6281323", "GSM6281324", "GSM6281325", "GSM6281326"
)

samples <- file.path(directory, sample_ids, "outs")

# barcodes
file_barcodes <- file.path(samples[1], "raw_feature_bc_matrix", "barcodes.tsv")
df_barcodes <- read.csv(file_barcodes, sep = "\t", header = FALSE, 
                        col.names = c("barcode_id"))

# counts
file_counts <- file.path(samples[1], "raw_feature_bc_matrix", "matrix.mtx")
counts <- readMM(file = file_counts)

#  features
file_features <- file.path(samples[1], "raw_feature_bc_matrix", "features.tsv")
df_features <- read.csv(file_features, sep = "\t", header = FALSE, 
                        col.names = c("gene_id", "gene_name", "feature_type"))

# tissue map
file_tisspos <- file.path(samples[1], "spatial", "tissue_positions_list.csv")
tisspos <- read.csv(file_tisspos, header = FALSE, 
                    col.names=c("barcode_id", "in_tissue", "array_row", "array_col",
                                "pxl_row_in_fullres", "pxl_col_in_fullres"))

image10x <- Read10X_Image(
  as.character(dirname(file_tisspos)),
  image.name = "tissue_lowres_image.png",
  filter.matrix = TRUE
)

# Load into seurat 
so <- CreateSeuratObject(
  counts = Read10X(data.dir=file.path(samples[3], "raw_feature_bc_matrix")),
  assay = 'Spatial',
  min.cells = 3,
  min.features = 200
)


image <- image10x[Cells(x = so)]
DefaultAssay(object = image) <- 'Spatial'
so[['Slice']] <- image


generateSeuratSpatial <- function(basedir) {
  file_tisspos <- file.path(basedir, "spatial", "tissue_positions_list.csv")
  image10x <- Read10X_Image(
    as.character(dirname(file_tisspos)),
    image.name = "tissue_lowres_image.png",
    filter.matrix = TRUE
  )

  so <- CreateSeuratObject(
    counts = Read10X(data.dir=file.path(basedir, "raw_feature_bc_matrix")),
    assay = 'Spatial',
    min.cells = 3,
    min.features = 200
  )

  image <- image10x[Cells(x = so)]
  DefaultAssay(object = image) <- 'Spatial'
  so[['Slice']] <- image
  
  return(so)
}

thymus.slices <- lapply(samples, generateSeuratSpatial)
thymus.slices <- lapply(thymus.slices, function(x) {
  x <- SCTransform(x, assay = "Spatial", verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = 2000, verbose = FALSE)
  return(x)
})

thymus.anchors    <- FindIntegrationAnchors(object.list = thymus.slices,
                                            dims = 1:30)

thymus <- IntegrateData(anchorset = thymus.anchors, dims = 1:30)
#thymus.merge <- Reduce(merge, thymus.slices)

#plot1 <- VlnPlot(so, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
#plot2 <- SpatialFeaturePlot(so, features = "nCount_Spatial") + theme(legend.position = "right")
#wrap_plots(plot1, plot2)

#so <- SCTransform(so, assay = "Spatial", verbose = FALSE)
#SpatialFeaturePlot(so, features = c("CD1D"), pt.size.factor = 1)

DefaultAssay(thymus) <- "integrated"
thymus <- ScaleData(thymus, verbose = FALSE)
thymus <- RunPCA(thymus, npcs = 30, verbose = FALSE)
thymus <- RunUMAP(thymus, reduction = "pca", dims = 1:30, verbose = FALSE)

#so <- RunPCA(so, assay = "SCT", verbose = FALSE)
thymus <- FindNeighbors(thymus, reduction = "pca", dims = 1:30)
thymus <- FindClusters(thymus, verbose = FALSE, algorithm = 1, resolution = 0.4)
#so <- RunUMAP(so, reduction = "pca", dims = 1:30)

p1 <- DimPlot(thymus, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(thymus, images="Slice.5", label = TRUE, label.size = 3,
                     pt.size.factor = 1)

p1/p2

so <- FindSpatiallyVariableFeatures(so, assay = "SCT",
                                    features = VariableFeatures(so)[1:1000],
                                    selection.method = "markvariogram")

top.features <- head(SpatiallyVariableFeatures(so, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(so, features = top.features, ncol = 3, alpha = c(0.1, 1))
