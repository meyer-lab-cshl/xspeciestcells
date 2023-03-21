#################
## libraries ####
#################
library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(harmony)
library(pals)

#################
## functions ####
#################

# Preprocessing function
preprocess <- function(seurobj, ndim=20, res=0.8, colvalues, celltype, 
                       approxPCA=FALSE){
  print(seurobj)
  cat("Cells:", unique(seurobj$group.ident), "\n")
  
  print("Normalize")
  seurobj <- NormalizeData(seurobj)
  print("HVG")
  seurobj <- FindVariableFeatures(seurobj, selection.method="vst",
                                  nfeatures=2000)
  print("Scale")
  seurobj <- ScaleData(seurobj)
  
  print("PCA")
  seurobj <- RunPCA(seurobj, npcs=50, verbose=FALSE, approx=approxPCA)
  
  print("Harmony")
  seurobj <- RunHarmony(seurobj, reduction = "pca", max.iter.harmony=30,
                        group.by.vars = c("Method", "Batch"))
  print(ElbowPlot(seurobj, ndims=50, reduction="harmony"))
  
  print("UMAP")
  seurobj <- RunUMAP(seurobj, reduction = "harmony", dims = 1:ndim,
                     n.neighbors=50, reduction.key = "UMAP50_")
  
  print("Cluster")
  seurobj <- FindNeighbors(seurobj, reduction = "harmony", dims = 1:ndim)
  seurobj <- FindClusters(seurobj, resolution = res, random.seed = 0)
  
  # Visualization sanity check
  print(DimPlot(seurobj, reduction = "umap", group.by="seurat_clusters",
                label = TRUE, repel = TRUE) +
          #scale_color_manual(values=colvalues) +
          labs(title=celltype))
  
  return(seurobj)
}

############
## data ####
############

colvalues <- colorRampPalette(brewer.pal(12,"Paired"))(18)
all_cells <- readRDS('data/seurat_filtered_harmony_02_15_23.RDS')
print(all_cells) # 78,607 cells

# filter for thymus only:
thymus <- subset(all_cells,
                      subset=Tissue=="Thymus")
print(thymus) # 37,369 cells

# Make raw object
thymus.filt <- CreateSeuratObject(counts=thymus@assays$RNA@counts,
                                  meta.data = thymus@meta.data)

print(thymus.filt) # 37,369 cells

################
## analysis ####
################

## find marker genes for cell type specific clusters ####

## NKT CELLS ####
thymus.nkt <- preprocess(seurobj = subset(thymus.filt,
                                          subset=group.ident=="NKT_Thymus"),
                         celltype = "NKTthy", colvalues=colvalues, ndim=10, res=0.4)

# Annotate
thymus.nkt@meta.data$cell_annot <- case_when(
  thymus.nkt@meta.data$seurat_clusters == 0 ~ "thyNKT_ccr7",
  thymus.nkt@meta.data$seurat_clusters == 1 ~ "thyNKT_effFOSJUN",
  thymus.nkt@meta.data$seurat_clusters == 2 ~ "thyNKT_ccr9",
  thymus.nkt@meta.data$seurat_clusters == 3 ~ "thyNKT_cd8aa",
  thymus.nkt@meta.data$seurat_clusters == 4 ~ "thyNKT_effector",
  #seur.nkt@meta.data$seurat_clusters == 5 ~ "thyNKT_IFNsig"
)

p_nkt <- DimPlot(thymus.nkt, group.by="cell_annot", label=T, repel=T) +
  theme(legend.position="none") + 
  scale_color_manual(values=colvalues) +
  labs(title="Thymic NKT")

# Markers NKT
markers.nkt <- FindAllMarkers(thymus.nkt, only.pos=TRUE,
                              min.pct=0.25, logfc.threshold=0.5) %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)

markers.nkt <- markers.nkt %>%
  select(avg_log2FC, p_val_adj, cluster, gene) %>%
  left_join(thymus.nkt@meta.data %>%
              select(cell_annot, RNA_snn_res.0.3) %>%
              dplyr::rename(cluster=RNA_snn_res.0.3) %>%
              distinct(),
            by="cluster") %>%
  relocate(gene, cell_annot) %>%
  as.data.frame


## MAIT CELLS ####
thymus.mait <- preprocess(seurobj = subset(thymus.filt,
                                           subset=group.ident=="MAIT_Thymus"),
                          celltype = "MAITthy", colvalues=colvalues,
                          ndim=10, res=0.2)

thymus.mait@meta.data$cell_annot <- case_when(
  thymus.mait@meta.data$seurat_clusters == 0 ~ "thyMAIT_ccr9",
  thymus.mait@meta.data$seurat_clusters == 1 ~ "thyMAIT_ccr7",
  thymus.mait@meta.data$seurat_clusters == 2 ~ "thyMAIT_effector",
  thymus.mait@meta.data$seurat_clusters == 3 ~ "thyMAIT_tbc",
  thymus.mait@meta.data$seurat_clusters == 4 ~ "thyMAIT_cd8aa",
  thymus.mait@meta.data$seurat_clusters == 5 ~ "thyMAIT_DP",
  thymus.mait@meta.data$seurat_clusters == 6 ~ "thyMAIT_IFNsig"
)

p_mait <- DimPlot(thymus.mait, group.by="cell_annot", label=T, repel=T) +
  theme(legend.position="none") +
  scale_color_manual(values=colvalues) +
  labs(title="Thymic MAIT")

# Markers MAIT
markers.mait <- FindAllMarkers(thymus.mait, only.pos=TRUE,
                               min.pct=0.25, logfc.threshold=0.5) %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)

markers.mait <- markers.mait %>%
  select(avg_log2FC, p_val_adj, cluster, gene) %>%
  left_join(thymus.mait@meta.data %>%
              select(cell_annot, RNA_snn_res.0.2) %>%
              dplyr::rename(cluster=RNA_snn_res.0.2) %>%
              distinct(),
            by="cluster") %>%
  relocate(gene, cell_annot)



## CD4 CELLS ####
thymus.cd4 <- preprocess(seurobj = subset(thymus.filt,
                                          subset=group.ident=="CD4_Thymus"), 
                       celltype = "CD4thy", ndim=10, res=0.2)

# Annotate
thymus.cd4@meta.data$cell_annot <- case_when(
  thymus.cd4@meta.data$seurat_clusters == 0 ~ "thyCD4_ccr9",
  thymus.cd4@meta.data$seurat_clusters == 1 ~ "thyCD4_ccr7",
  thymus.cd4@meta.data$seurat_clusters == 2 ~ "thyCD4_ISP",
  thymus.cd4@meta.data$seurat_clusters == 3 ~ "thyCD4_Treg",
  thymus.cd4@meta.data$seurat_clusters == 4 ~ "thyCD4_DPp",
  thymus.cd4@meta.data$seurat_clusters == 5 ~ "thyCD4_Tagonist",
  thymus.cd4@meta.data$seurat_clusters == 6 ~ "thyCD4_DPq"
)

p_cd4 <- DimPlot(thymus.cd4, group.by="cell_annot", label=TRUE, repel=TRUE) +
  theme(legend.position="none") +
  scale_color_manual(values=colvalues) +
  labs(title="Thymic CD4")

# Markers CD4
markers.cd4 <- FindAllMarkers(thymus.cd4, only.pos=TRUE,
                              min.pct=0.25, logfc.threshold=0.5) %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)

markers.cd4 <- markers.cd4 %>%
  select(avg_log2FC, p_val_adj, cluster, gene) %>%
  left_join(thymus.cd4@meta.data %>%
              select(cell_annot, RNA_snn_res.0.2) %>%
              dplyr::rename(cluster=RNA_snn_res.0.2) %>%
              distinct(),
            by="cluster") %>%
  relocate(gene, cell_annot)


## CD8 CELLS ####
thymus.cd8 <- preprocess(seurobj = subset(thymus.filt,
                                          subset=group.ident=="CD8_Thymus"),
                       celltype = "CD8thy", ndim=10, res=0.2)

# Annotate
thymus.cd8@meta.data$cell_annot <- case_when(
  thymus.cd8@meta.data$seurat_clusters == 0 ~ "thyCD8_ccr7",
  thymus.cd8@meta.data$seurat_clusters == 1 ~ "thyCD8_ccr9", 
  thymus.cd8@meta.data$seurat_clusters == 2 ~ "thyCD8_cd8aa1",
  thymus.cd8@meta.data$seurat_clusters == 3 ~ "thyCD8_idk",
  thymus.cd8@meta.data$seurat_clusters == 4 ~ "thyCD8_cd8aa2",
  thymus.cd8@meta.data$seurat_clusters == 5 ~ "thyCD8_DP"
)

p_cd8 <- DimPlot(thymus.cd8, group.by="cell_annot", label=TRUE, repel=TRUE) +
  theme(legend.position="none") +
  scale_color_manual(values=colvalues) +
  labs(title="Thymic CD8")

# Markers CD8
markers.cd8 <- FindAllMarkers(thymus.cd8, only.pos=TRUE,
                              min.pct=0.25, logfc.threshold=0.5) %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)

markers.cd8 <- markers.cd8 %>%
  select(avg_log2FC, p_val_adj, cluster, gene) %>%
  left_join(thymus.cd8@meta.data %>%
              select(cell_annot, RNA_snn_res.0.2) %>%
              dplyr::rename(cluster=RNA_snn_res.0.2) %>%
              distinct(),
            by="cluster") %>%
  relocate(gene, cell_annot)

## GD CELLS ####
thymus.gd <- preprocess(seurobj = subset(thymus.filt,
                                       subset=group.ident=="GD_Thymus"),
                      celltype = "GDthy", ndim=10, res=0.2)

# Annotate
thymus.gd@meta.data$cell_annot <- case_when(
  thymus.gd@meta.data$seurat_clusters == 0 ~ "thyGDT_IFNsig",
  thymus.gd@meta.data$seurat_clusters == 1 ~ "thyGDT_immat",
  thymus.gd@meta.data$seurat_clusters == 2 ~ "thyGDT_ccr9",
  thymus.gd@meta.data$seurat_clusters == 3 ~ "thyGDT_effector",
  thymus.gd@meta.data$seurat_clusters == 4 ~ "thyGDT_immat_cycl",
  thymus.gd@meta.data$seurat_clusters == 5 ~ "thyGDT_DP"
)

p_gd <- DimPlot(thymus.gd, group.by="cell_annot", label=TRUE, repel=TRUE) +
  theme(legend.position="none") + 
  scale_color_manual(values=colvalues) +
  labs(title="Thymic GDT")

# Markers GDT
markers.gd <- FindAllMarkers(thymus.gd, only.pos=TRUE, min.pct=0.25,
                              logfc.threshold=0.5) %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)

markers.gd <- markers.gd %>%
  select(avg_log2FC, p_val_adj, cluster, gene) %>%
  left_join(thymus.gd@meta.data %>%
              select(cell_annot, RNA_snn_res.0.2) %>%
              dplyr::rename(cluster=RNA_snn_res.0.2) %>%
              distinct(),
            by="cluster") %>%
  relocate(gene, cell_annot)

# Annotate integrated  based on clusters derived here
# combine all metadata:
subset_annotations <- c(thymus.nkt$cell_annot, thymus.mait$cell_annot,
                thymus.cd4$cell_annot, thymus.cd8$cell_annot,
                thymus.gd$cell_annot)

subset_annotations <- subset_annotations[match(colnames(thymus),
                                               names(subset_annotations))]
# check: sum(names(subset_annotations) == colnames(thymus)) == ncol(thymus)

thymus@meta.data$subset_annotations <- subset_annotations
saveRDS(thymus, "data/seurat_filtered_harmony_02_15_23_thymus_subset_annotations.RDS")


# combine hvg:
length_hvg <- length(VariableFeatures(thymus.nkt))
subset_hvgs <- tibble(celltypes=c(rep("NKT", length_hvg),
                        rep("MAIT", length_hvg),
                        rep("CD4", length_hvg),
                        rep("CD8", length_hvg),
                        rep("GD", length_hvg)),
                      Features=c(VariableFeatures(thymus.nkt),
                        VariableFeatures(thymus.mait),
                        VariableFeatures(thymus.cd4),
                        VariableFeatures(thymus.cd8),
                        VariableFeatures(thymus.gd))
)

readr::write_csv(subset_hvgs, "results/hvgs_thymic_populations.csv")

# combine all plots
p_all <- cowplot::plot_grid(p_cd4, p_cd8, p_nkt, p_gd, p_mait)
ggsave("results/plots/thymic_subsets.pdf", plot=p_all, width=9, height=6)

# Combine MAIT and NKT markers and export as .csv file
markers.subsets <- rbind(markers.nkt, markers.mait, markers.cd4, markers.cd8,
                           markers.gd) %>%
  ungroup() %>%
  select(-cluster)
table(markers.subsets$cell_annot)
readr::write_csv(markers.subsets, "results/top20markers_thymic_populations.csv")
