###
# Purpose: Quality filtering of seurat object
# Date: Dec 8th 2022
# Author: Salomé Carcy
###


## 1. IMPORT ####

library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)
library(Seurat)
library(RColorBrewer)
# library(harmony)
# library(ggalluvial)


# Import seurat object
seur.human <- readRDS("~/Projects/HumanThymusProject/data/raw_data/human_data/filtered_seurat_harmony.11.22.22.RDS")
print(seur.human) # 81,378 cells (it's the whole seurat object)

# Quick visualization
DimPlot(seur.human, reduction="UMAP_50", group.by="new_clusters")+
  scale_color_manual(values=colorRampPalette(brewer.pal(12,"Paired"))(18))




## 2. GET RAW SEURAT OBJECT (optional) ####

# # Get raw counts
# counts <- seur.human@assays$RNA@counts
# metadata <- seur.human@meta.data
# 
# # Remove clustering information in metadata, except for relevant clusters
# metadata[,c("RNA_snn_res.0.7", "seurat_clusters", "RNA_snn_res.1.2")] <- NULL
# colnames(metadata)[14] <- "clusters_laurent" # rename Laurent's clusters to "clusters_laurent" (to keep as reference)
# 
# # Create seurat object
# seur.human <- CreateSeuratObject(counts=counts, meta.data=metadata, min.cells=10) # keep genes expressed in at least 10 cells




## 3. REMOVE CELLS WITH LOW MITOCHONDRIAL CONTENT ####

# Look at qc measures per batch & donor
VlnPlot(seur.human, features="percent.mt",   group.by="Batch", split.by="Donor") + labs(x="Batch", fill="Donor")
# VlnPlot(seur.human, features="nFeature_RNA", group.by="Batch", split.by="Donor") + labs(x="Batch", fill="Donor")
# VlnPlot(seur.human, features="nCount_RNA",   group.by="Batch", split.by="Donor") + labs(x="Batch", fill="Donor")

# Split seurat object to find percent.mt outliers per batch (& subset seurat object)
seur.list <- SplitObject(seur.human, split.by = "Batch")
seur.list <- lapply(X = seur.list, FUN = function(x) {
  threshold <- attr(scuttle::isOutlier(x@meta.data$percent.mt, type="higher"), "thresholds")["higher"]
  print(threshold)
  x <- subset(x, subset= percent.mt<threshold)
  # x@meta.data$highmt <- ifelse(x@meta.data$percent.mt > threshold, TRUE, FALSE)
  # return(x)
})
# Combine seurat object back together
seur.filt <- purrr::reduce(seur.list, function(x, y){
  merge(x = x, y = y, add.cell.ids = NULL, merge.data = T)
})

# Sanity check
VlnPlot(seur.filt, features="percent.mt", group.by="Batch", split.by="Donor") + labs(x="Batch", fill="Donor", title="post-filtering")
# ggsave("~/Projects/HumanThymusProject/data/human-thymus/HumanData_08_QualityControl/mt_filtered.jpeg", width=10, height=6)




## 4. REMOVE DOUBLETS ####

library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(scater)

# Transform to SCE object
counts   <- seur.filt@assays$RNA@counts
metadata <- seur.filt@meta.data
umap <- seur.human[["UMAP_50"]]@cell.embeddings
umap <- umap[rownames(umap) %in% rownames(seur.filt@meta.data),] # keep only cells that passed MT filt
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata,
                            reducedDims = list(umap = umap)) # 80,614 cells


# note: expected doublet rate is 5% (for <10k cells)
# Run scDblFinder (with seurat clusters)
set.seed(123)
sce.srtclusters <- scDblFinder(sce, clusters="new_clusters", samples="Batch", BPPARAM=MulticoreParam(4))

# Run scDblFinder (with cell lineage)
set.seed(123)
sce.cellident <- scDblFinder(sce, clusters="cell.ident", samples="Batch", BPPARAM=MulticoreParam(4))

plotReducedDim(sce.srtclusters, dimred="umap", colour_by="scDblFinder.class", point_size=0.1) + labs(title="scDblFinder on seurat clusters")+
  scale_color_manual(values=c("#bdbdbd", "#de2d26"), name="") |
  plotReducedDim(sce.cellident, dimred="umap", colour_by="scDblFinder.class", point_size=0.1) + labs(title="scDblFinder on cell identity")+
  scale_color_manual(values=c("#bdbdbd", "#de2d26"), name="")
# ggsave("~/Projects/HumanThymusProject/data/human-thymus/HumanData_08_QualityControl/scDblFinder.jpeg", width=10, height=5)


# Check how much the two methods agree
table(rownames(colData(sce.srtclusters)) == rownames(colData(sce.cellident))) # sanity check
dbl <- data.frame(srtclusters=colData(sce.srtclusters)$scDblFinder.class,
                  cellident=colData(sce.cellident)$scDblFinder.class)
table(dbl)

# Get vector of cells that are considered as doublets in both cases
doublets <- intersect(rownames(colData(sce.srtclusters)[colData(sce.srtclusters)$scDblFinder.class=="doublet",]),
                      rownames(colData(sce.cellident)[colData(sce.cellident)$scDblFinder.class=="doublet",]))
length(doublets) # 3640 cells

# Remove doublets
sce.filt <- sce[,!colnames(sce) %in% doublets]




## 5. GET SEURAT OBJECT AND SAVE ####
seur.filt2 <- CreateSeuratObject(counts=assays(sce.filt)$counts, meta.data=as.data.frame(colData(sce.filt)))
# seur.filt2@meta.data[,c("RNA_snn_res.0.7", "seurat_clusters", "RNA_snn_res.1.2")] <- NULL
# colnames(seur.filt2@meta.data)[14] <- "clusters_laurent" # rename Laurent's clusters to "clusters_laurent" (to keep as reference)

# Save
saveRDS(seur.filt2, "~/Projects/HumanThymusProject/data/raw_data/human_data/seur_filt_23-02-08.rds")
