###
# Purpose: Get human raw data 
# Date: Oct 18th 2022
# Author: Salomé Carcy
###


#### IMPORT ####
library(Seurat)

# Data
seur.human <- readRDS("~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/human_data/filtered_seurat_Harmony_07-22-22.RDS")
print(seur.human) # 79,801 cells and 34,778 genes




#### EXTRACT RAW DATA ####

# Get the raw counts
rawcounts <- seur.human@assays$RNA@counts

# Take quick peek at the metadata
head(seur.human@meta.data)

table(seur.human$orig.ident)
table(seur.human$cell.ident)
table(seur.human$Sex)
table(seur.human$Age_in_weeks)
table(seur.human$Donor)
table(seur.human$Batch)
table(seur.human$Method)
table(seur.human$Tissue)
table(seur.human$group.ident)
table(seur.human$group.ident, seur.human$cell.ident) # sanity check
table(seur.human$TCRa_g_chain)
table(seur.human$TCRa_g_clonotype)
table(seur.human$TRAV10_TRAJ18, seur.human$cell.ident)

# Metadata of interest
cols <- colnames(seur.human@meta.data)[!colnames(seur.human@meta.data) %in% c("nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_SCT", "nFeature_SCT",
                                                                              "SCT_snn_res.0.7", "Test_cell_ident_for_cLISI", "SCT_snn_res.0.9",
                                                                              "SCT_snn_res.1", "new_clusters")]
metadata <- seur.human@meta.data[,cols]
head(metadata)
colnames(metadata)[10] <- "clusters_alldata"


# Create new seurat object
rawseur <- CreateSeuratObject(counts=rawcounts, meta.data=metadata, assay="RNA")
saveRDS(rawseur, "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/06_HumanData/seurat_raw_hu.rds")
