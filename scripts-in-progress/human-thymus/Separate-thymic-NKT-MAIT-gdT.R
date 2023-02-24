###
# Purpose: Generate seurat objects for human thymic NKT, MAIT and gdT cells (with clusters)
# Date: Nov 7th 2022
###


#### IMPORT ####
# Librairies
library(Seurat)
library(harmony)
library(dplyr)

# Data
seur.human <- readRDS("~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/human_data/filtered_seurat_Harmony_07-22-22.RDS")
print(seur.human) # 79,801 cells and 34,778 genes




#### DEFINE FUNCTIONS ####

# This function will preprocess the data (normalization, integration, etc.)
Preprocess <- function(seuratobject, cluster_res=0.5, harmony_vars= c("Sex", "Method"), pcs=15, printComputation=T){
  
  # Extract the counts data (from RNA assay) and metadata
  counts <- GetAssayData(object = seuratobject, slot = "counts", assay="RNA")
  metadata <- seuratobject@meta.data
  metadata[, c("nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.7", "seurat_clusters", "SCT_snn_res.0.9", "SCT_snn_res.1")] <- NULL # remove some cluster columns to avoid any confusion
  colnames(metadata)[14] <- "clusters_bigUMAP" # rename "new_clusters" column
  # Create new seurat object with only the counts and the "cleaned" metadata
  initseurat <- CreateSeuratObject(counts = counts, meta.data = metadata, min.cells=20) # keep only genes expressed in at least 20 cells
  
  cat("Initial seurat object:\n")
  print(initseurat) # this is for sanity check
  
  
  cat("\n1-Perform SCTransform on data split by sample/batch \n")
  # Normalize with SCTransform by splitting object (per batch/sample)
  seur.list <- SplitObject(initseurat, split.by = "orig.ident")
  for(i in 1:length(seur.list)){
    cat("...SCTransform on", names(seur.list)[i], "\n")
    seur.list[[i]] <-  SCTransform(seur.list[[i]],
                                   vars.to.regress = 'percent.mt',
                                   assay="RNA",
                                   new.assay.name="SCT",
                                   do.correct.umi=TRUE,
                                   ncells=5000,
                                   seed.use = 1448145,
                                   verbose=printComputation)
  }
  seur.SCT <- merge(seur.list[[1]], y = unlist(seur.list[-1], use.names=FALSE), merge.data = TRUE) # merge back the data into one seurat object
  
  
  # Define shared HVGs between the samples & save them as VariableFeatures in the seur.SCT (combined data) object
  sct.hvg  <- SelectIntegrationFeatures(object.list = seur.list, nfeatures = 3000)
  VariableFeatures(seur.SCT) <- sct.hvg
  
  
  # Run PCA (don't scale data when using SCTransform)
  cat("\n2-Run PCA pre-integration")
  seur.SCT <- RunPCA(object = seur.SCT,
                     assay = "SCT",
                     npcs = 50,
                     seed.use=42,
                     weight.by.var = TRUE,
                     verbose=printComputation)
  p1 <- DimPlot(seur.SCT, reduction = "pca", group.by = "orig.ident") + ggtitle("PCA pre-integration")
  print(ElbowPlot(seur.SCT))
  cat("\nSeurat object pre-integration:\n")
  print(seur.SCT) # sanity check
  
  
  # Run Harmony for integration
  cat("\n3-Run Harmony\n")
  seur.harmony <- RunHarmony(object = seur.SCT,
                             assay.use = "SCT",
                             reduction = "pca",
                             dims.use = 1:50,
                             group.by.vars = harmony_vars,
                             max.iter.harmony=30,
                             plot_convergence = printComputation)
  cat("\nSeurat object post-integration:\n")
  print(seur.harmony) # sanity check
  
  
  # Run UMAP & clustering
  cat("\n4-Run UMAP & clustering\n")
  seur.harmony <- RunUMAP(object = seur.harmony,
                          dims = 1:pcs, # default 15
                          assay = "SCT",
                          reduction = "harmony",
                          seed.use=42,
                          verbose=printComputation)
  seur.harmony <- FindNeighbors(object = seur.harmony,
                                dims = 1:pcs, # default 15
                                k.param = 20,
                                assay = "SCT",
                                reduction = "harmony",
                                verbose=printComputation)
  seur.harmony <- FindClusters(object = seur.harmony,
                               resolution = cluster_res,
                               algorithm=1,
                               random.seed=42,
                               verbose=printComputation)
  
  
  # Show PCA after integration
  p2 <- DimPlot(seur.harmony, reduction = "harmony", group.by = "orig.ident") + ggtitle("PCA post-integration")
  print(p1 | p2)
  
  return(seur.harmony)
}




# __________________
#### NKT CELLS ####

# Isolate thymic NKT cells
NKT_Thymus <- subset(seur.human, ident = "NKT_Thymus")
print(NKT_Thymus) # 2551 cells

# Quick QC
# FeatureScatter(NKT_Thymus, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by="Batch")
# FeatureScatter(NKT_Thymus, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Batch")

# Preprocess
NKT_integrated <- Preprocess(seuratobject = NKT_Thymus, harmony_vars = c("Sex", "Method"), cluster_res = 0.4, printComputation = F)
DimPlot(NKT_integrated, dims = 1:2, reduction = "umap", group.by = "SCT_snn_res.0.4", label = T) + ggtitle("Thymic NKT cells")

# Re-name clusters
NKT_integrated$new_clusters_NKT <- case_when(
  NKT_integrated$SCT_snn_res.0.4 == '3' ~ '0',
  NKT_integrated$SCT_snn_res.0.4 == '2' ~ '1',
  NKT_integrated$SCT_snn_res.0.4 == '0' ~ '2',
  NKT_integrated$SCT_snn_res.0.4 == '1' ~ '3',
  NKT_integrated$SCT_snn_res.0.4 == '4' ~ '4',
)
DimPlot(NKT_integrated, dims = 1:2, reduction = "umap", group.by = "new_clusters_NKT", label = T) +
  ggtitle("New cluster numbers")

# Save seurat object as .rds file
# saveRDS(NKT_integrated, "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/08_SplitData/Thymic_NKT_22-11-07.rds")




# __________________
#### MAIT CELLS ####

# Isolate thymic MAIT cells
MAIT_Thymus <- subset(seur.human, ident = "MAIT_Thymus")
print(MAIT_Thymus) # 4689 cells

# Quick QC
# FeatureScatter(MAIT_Thymus, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by="orig.ident")
# FeatureScatter(MAIT_Thymus, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="orig.ident")

# Preprocess
MAIT_integrated <- Preprocess(seuratobject = MAIT_Thymus, harmony_vars = c("Sex", "Method"), cluster_res = 0.3, printComputation = F)
DimPlot(MAIT_integrated, dims = 1:2, reduction = "umap", group.by = "SCT_snn_res.0.3", label = T) + ggtitle("Thymic MAIT cells")

# Re-name clusters
MAIT_integrated$new_clusters_MAIT <- case_when(
  MAIT_integrated$SCT_snn_res.0.3 == '5' ~ '0', # DPs
  MAIT_integrated$SCT_snn_res.0.3 == '3' ~ '1',
  MAIT_integrated$SCT_snn_res.0.3 == '1' ~ '2',
  MAIT_integrated$SCT_snn_res.0.3 == '0' ~ '3',
  MAIT_integrated$SCT_snn_res.0.3 == '4' ~ '4',
  MAIT_integrated$SCT_snn_res.0.3 == '6' ~ '5',
  MAIT_integrated$SCT_snn_res.0.3 == '2' ~ '6',
)
DimPlot(MAIT_integrated, dims = 1:2, reduction = "umap", group.by = "new_clusters_MAIT", label = T) +
  ggtitle("New cluster numbers")

# Save seurat object as .rds file
# saveRDS(MAIT_integrated, "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/08_SplitData/Thymic_MAIT_22-11-07.rds")




# _________________
#### gdT CELLS ####

# Isolate thymic MAIT cells
gdT_Thymus <- subset(seur.human, ident = "GD_Thymus")
print(gdT_Thymus) # 2981 cells

# Quick QC
# FeatureScatter(gdT_Thymus, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by="orig.ident")
# FeatureScatter(gdT_Thymus, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="orig.ident")

# Preprocess
gdT_integrated <- Preprocess(seuratobject = gdT_Thymus, harmony_vars = c("Sex", "Method"), cluster_res = 0.4, printComputation = F)
DimPlot(gdT_integrated, dims = 1:2, reduction = "umap", group.by = "SCT_snn_res.0.4", label = T) + ggtitle("Thymic gdT cells")

# Re-name clusters
gdT_integrated$new_clusters_gdT <- case_when(
  gdT_integrated$SCT_snn_res.0.4 == '0' ~ '0',
  gdT_integrated$SCT_snn_res.0.4 == '2' ~ '1',
  gdT_integrated$SCT_snn_res.0.4 == '3' ~ '2',
  gdT_integrated$SCT_snn_res.0.4 == '1' ~ '3',
  gdT_integrated$SCT_snn_res.0.4 == '4' ~ '4',
  gdT_integrated$SCT_snn_res.0.4 == '5' ~ '5',
  gdT_integrated$SCT_snn_res.0.4 == '6' ~ '6',
  gdT_integrated$SCT_snn_res.0.4 == '7' ~ '7',
)
DimPlot(gdT_integrated, dims = 1:2, reduction = "umap", group.by = "new_clusters_gdT", label = T) +
  ggtitle("New cluster numbers")

# Save seurat object as .rds file
# saveRDS(gdT_integrated, "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/08_SplitData/Thymic_GDT_22-11-07.rds")