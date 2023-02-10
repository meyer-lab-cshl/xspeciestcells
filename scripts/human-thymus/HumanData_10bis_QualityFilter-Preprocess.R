###
# Purpose: Quality filtering of seurat object
# Date: Feb 9th 2023
# Author: Salomé Carcy
###


## 1. IMPORT ####

library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(harmony)
# library(ggalluvial)


# Import seurat object
seur.human <- readRDS("~/Projects/HumanThymusProject/data/raw_data/human_data/filtered_seurat_harmony.11.22.22.RDS")
print(seur.human) # 81,378 cells (it's the whole seurat object)

# Quick visualization
DimPlot(seur.human, reduction="UMAP_50", group.by="new_clusters")+
  scale_color_manual(values=colorRampPalette(brewer.pal(12,"Paired"))(18))




## 2. REMOVE CELLS WITH LOW MITOCHONDRIAL CONTENT ####

# Look at qc measures per batch & donor
VlnPlot(seur.human, features="percent.mt",   group.by="Batch", split.by="Donor") + labs(x="Batch", fill="Donor")

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

seur.filt@meta.data[,c("RNA_snn_res.0.7", "seurat_clusters", "RNA_snn_res.1.2")] <- NULL
colnames(seur.filt@meta.data)[14] <- "clusters_laurent" # rename Laurent's clusters to "clusters_laurent" (to keep just in case)

# Save object (but there isn't UMAP, PCA, etc. saved in it)
saveRDS(seur.filt, "~/Projects/HumanThymusProject/data/CLUSTER/seurat-integration/input/seur_filt_23-02-09.rds")




## 3. PREPROCESSING INTEGRATED OBJECT ####

# Preprocessing function
preprocess <- function(seurobj, ndim=20, res=0.8, approxPCA=TRUE){
  print(seurobj)
  cat("Cells:", unique(seurobj$group.ident), "\n")
  
  print("Normalize")
  seurobj <- NormalizeData(seurobj)
  print("HVG")
  seurobj <- FindVariableFeatures(seurobj, selection.method="vst", nfeatures=2000)
  print(seurobj)
  print("Scale")
  seurobj <- ScaleData(seurobj)
  print("PCA")
  seurobj <- RunPCA(seurobj, npcs=50, verbose=FALSE, approx=approxPCA)
  print("Harmony")
  seurobj <- RunHarmony(seurobj, reduction = "pca", max.iter.harmony=30, group.by.vars = c("Method", "Batch"))
  # DimPlot(seurobj, reduction="pca", group.by="Batch") | DimPlot(seurobj, reduction="harmony", group.by="Batch")
  # DimPlot(seurobj, reduction="pca", group.by="Method") | DimPlot(seurobj, reduction="harmony", group.by="Method")
  print(ElbowPlot(seurobj, ndims=50, reduction="harmony"))
  print("UMAP")
  seurobj <- RunUMAP(seurobj, reduction = "harmony", dims = 1:ndim, n.neighbors=50, reduction.key = "UMAP50_")
  print("Cluster")
  seurobj <- FindNeighbors(seurobj, reduction = "harmony", dims = 1:ndim)
  seurobj <- FindClusters(seurobj, resolution = res, random.seed = 0)
  
  return(seurobj)
}

# Preprocessing on integrated object
seur.total <- preprocess(seur.filt, res=1.2, ndim=20)

# Re-order clusters
seur.total@meta.data$seurat_clusters <- case_when(
  seur.total@meta.data$RNA_snn_res.1.2 == 0 ~ 14,
  seur.total@meta.data$RNA_snn_res.1.2 == 1 ~ 15,
  seur.total@meta.data$RNA_snn_res.1.2 == 2 ~ 5,
  seur.total@meta.data$RNA_snn_res.1.2 == 3 ~ 10,
  seur.total@meta.data$RNA_snn_res.1.2 == 4 ~ 8,
  seur.total@meta.data$RNA_snn_res.1.2 == 5 ~ 13,
  seur.total@meta.data$RNA_snn_res.1.2 == 6 ~ 11,
  seur.total@meta.data$RNA_snn_res.1.2 == 7 ~ 12,
  seur.total@meta.data$RNA_snn_res.1.2 == 8 ~ 4,
  seur.total@meta.data$RNA_snn_res.1.2 == 9 ~ 7,
  seur.total@meta.data$RNA_snn_res.1.2 == 10 ~ 9,
  seur.total@meta.data$RNA_snn_res.1.2 == 11 ~ 16,
  seur.total@meta.data$RNA_snn_res.1.2 == 12 ~ 6,
  seur.total@meta.data$RNA_snn_res.1.2 == 13 ~ 0,
  seur.total@meta.data$RNA_snn_res.1.2 == 14 ~ 3,
  seur.total@meta.data$RNA_snn_res.1.2 == 15 ~ 1,
  seur.total@meta.data$RNA_snn_res.1.2 == 16 ~ 17,
  seur.total@meta.data$RNA_snn_res.1.2 == 17 ~ 2,
  seur.total@meta.data$RNA_snn_res.1.2 == 18 ~ 18,
  seur.total@meta.data$RNA_snn_res.1.2 == 19 ~ 1
)

# Quick visualization
# seur.total@meta.data$clusters_laurent <- as.numeric(seur.total@meta.data$clusters_laurent) # if clusters_laurent appears as character
DimPlot(seur.total, reduction="umap", group.by="clusters_laurent")+
  scale_color_manual(values=colorRampPalette(brewer.pal(12,"Paired"))(20)) |
DimPlot(seur.total, reduction="umap", group.by="seurat_clusters")+
  scale_color_manual(values=colorRampPalette(brewer.pal(12,"Paired"))(20))

# Save object
Idents(seur.total) <- "seurat_clusters"
saveRDS(seur.total, "~/Projects/HumanThymusProject/data/human-thymus/HumanData_10bis_QualityFilter-Preprocessing/seurat_filtered_harmony_02-09-23.RDS")




## 4. PREPROCESSING THYMIC LINEAGES SEPARATELY ####


### 4.1. THYMIC CD4 ####

seur.cd4 <- preprocess(seurobj = subset(seur.filt, subset=group.ident=="CD4_Thymus"), ndim=10, res=0.2)

# Annotate
# FindAllMarkers(seur.cd4, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.5) %>%
#   group_by(cluster) %>%
#   slice_max(n=5, order_by=avg_log2FC)
seur.cd4@meta.data$cell_annot <- case_when(
  seur.cd4@meta.data$seurat_clusters == 0 ~ "thyCD4_ccr9",
  seur.cd4@meta.data$seurat_clusters == 1 ~ "thyCD4_ccr7",
  seur.cd4@meta.data$seurat_clusters == 2 ~ "thyCD4_ISP",
  seur.cd4@meta.data$seurat_clusters == 3 ~ "thyCD4_DPp",
  seur.cd4@meta.data$seurat_clusters == 4 ~ "thyCD4_Treg",
  seur.cd4@meta.data$seurat_clusters == 5 ~ "thyCD4_Tagonist",
  seur.cd4@meta.data$seurat_clusters == 6 ~ "thyCD4_DPq"
)

# Visualize
DimPlot(seur.cd4, group.by="cell_annot", label=T, repel=T) + 
  theme(legend.position="none") + labs(title="Thymic CD4")


### 4.2. THYMIC CD8 ####

seur.cd8 <- preprocess(seurobj = subset(seur.filt, subset=group.ident=="CD8_Thymus"), ndim=10, res=0.2)

# Annotate
# markers <- FindAllMarkers(seur.cd8, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.5) %>%
#   group_by(cluster) %>%
#   slice_max(n=10, order_by=avg_log2FC)
seur.cd8@meta.data$cell_annot <- case_when(
  seur.cd8@meta.data$seurat_clusters == 0 ~ "thyCD8_ccr7",
  seur.cd8@meta.data$seurat_clusters == 1 ~ "thyCD8_ccr9", # ccr9?
  seur.cd8@meta.data$seurat_clusters == 2 ~ "thyCD8_cd8aa1",
  seur.cd8@meta.data$seurat_clusters == 3 ~ "thyCD8_idk",
  seur.cd8@meta.data$seurat_clusters == 4 ~ "thyCD8_cd8aa2",
  seur.cd8@meta.data$seurat_clusters == 5 ~ "thyCD8_DP"
)

# Visualize
DimPlot(seur.cd8, group.by="cell_annot", label=T, repel=T) +
                    theme(legend.position="none") + labs(title="Thymic CD8")


### 4.3. THYMIC MAIT ####

seur.mait <- preprocess(seurobj = subset(seur.filt, subset=group.ident=="MAIT_Thymus"), ndim=10, res=0.2)

# Annotate
# markers <- FindAllMarkers(seur.mait, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.5) %>%
#   group_by(cluster) %>%
#   slice_max(n=5, order_by=avg_log2FC)
seur.mait@meta.data$cell_annot <- case_when(
  seur.mait@meta.data$seurat_clusters == 0 ~ "thyMAIT_ccr7",
  seur.mait@meta.data$seurat_clusters == 1 ~ "thyMAIT_ccr9",
  seur.mait@meta.data$seurat_clusters == 2 ~ "thyMAIT_effector",
  seur.mait@meta.data$seurat_clusters == 3 ~ "thyMAIT_idk",
  seur.mait@meta.data$seurat_clusters == 4 ~ "thyMAIT_cd8aa",
  seur.mait@meta.data$seurat_clusters == 5 ~ "thyMAIT_DP",
  seur.mait@meta.data$seurat_clusters == 6 ~ "thyMAIT_IFNsig"
)

# Visualize
DimPlot(seur.mait, group.by="cell_annot", label=T, repel=T) +
  theme(legend.position="none") + labs(title="Thymic MAIT")


### 4.4. THYMIC NKT ####

seur.nkt <- preprocess(seurobj = subset(seur.filt, subset=group.ident=="NKT_Thymus"), ndim=10, res=0.3)

# Annotate
# markers <- FindAllMarkers(seur.nkt, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.5) %>%
#   group_by(cluster) %>%
#   slice_max(n=5, order_by=avg_log2FC)
seur.nkt@meta.data$cell_annot <- case_when(
  seur.nkt@meta.data$seurat_clusters == 0 ~ "thyNKT_ccr7",
  seur.nkt@meta.data$seurat_clusters == 1 ~ "thyNKT_effFOSJUN",
  seur.nkt@meta.data$seurat_clusters == 2 ~ "thyNKT_ccr9",
  seur.nkt@meta.data$seurat_clusters == 3 ~ "thyNKT_cd8aa",
  seur.nkt@meta.data$seurat_clusters == 4 ~ "thyNKT_effector",
  seur.nkt@meta.data$seurat_clusters == 5 ~ "thyNKT_IFNsig"
)

# Visualize
DimPlot(seur.nkt, group.by="cell_annot", label=T, repel=T) +
  theme(legend.position="none") + labs(title="Thymic NKT")


### 4.5. THYMIC GDT ####

seur.gd <- preprocess(seurobj = subset(seur.filt, subset=group.ident=="GD_Thymus"), ndim=10, res=0.2)

# Annotate
# markers <- FindAllMarkers(seur.gd, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.5) %>%
#   group_by(cluster) %>%
#   slice_max(n=10, order_by=avg_log2FC)
seur.gd@meta.data$cell_annot <- case_when(
  seur.gd@meta.data$seurat_clusters == 0 ~ "thyGDT_IFNsig",
  seur.gd@meta.data$seurat_clusters == 1 ~ "thyGDT_ccr9",
  seur.gd@meta.data$seurat_clusters == 2 ~ "thyGDT_immat",
  seur.gd@meta.data$seurat_clusters == 3 ~ "thyGDT_effector",
  seur.gd@meta.data$seurat_clusters == 4 ~ "thyGDT_immat_cycl",
  seur.gd@meta.data$seurat_clusters == 5 ~ "thyGDT_DP"
)

# Visualize
DimPlot(seur.gd, group.by="cell_annot", label=T, repel=T) +
  theme(legend.position="none") + labs(title="Thymic GDT")
