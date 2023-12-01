# Purpose:
# Author: Salomé Carcy
# Date:




# **************
# 1. IMPORT ####
# **************

# Import librairies
library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)
library(Seurat)
source("./scripts-final/colors_universal.R")

# Import integrated data
seur.human <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_08_28_23.RDS")
seur.ms <- readRDS("./data/cross-species/08_innateT_cross_species/Analysis_all_mouse-Tinn_filtered_seurat_MNN.rds")

# Import thymus data
seur.thym.nkt  <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.nkt.RDS")
seur.thym.mait <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.mait.RDS")
seur.thym.gdt  <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.gd.RDS")
seur.thym.cd4  <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.cd4.RDS")
seur.thym.cd8  <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.cd8.RDS")

# Import PBMC data
seur.pbmc.cd4  <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.CD4_03_16_23.RDS")
seur.pbmc.cd8  <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.CD8_03_16_23.RDS")
seur.pbmc.gdt  <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.GD_03_16_23.RDS")
seur.pbmc.mait <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.MAIT_03_16_23.RDS")
seur.pbmc.nkt  <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.NKT_03_16_23.RDS")



# **********************************************
# 2. PREPARE HUMAN INTEGRATED SEURAT OBJECT ####
# **********************************************

#___________________________________
## 2.1. Keep columns of interest ####

colnames(seur.human@meta.data)
seur_human_integrated_columns_to_keep <- c(
  "nCount_RNA", "nFeature_RNA", "percent.mt",
  "cell.ident", "group.ident",
  "Sex", "Age_in_weeks", "Donor", "Batch", "Method", "Tissue",
  "new_clusters",
  "TCR_Alpha_Gamma_V_gene_Dominant", "TCR_Beta_Delta_V_gene_Dominant", "TRAV10_TRAJ18", "TRAV1_TRAJ33", "TRAV1",
  "Egress_score1", "Effectorness1", "Naiveness1"
  )
seur.human@meta.data <- seur.human@meta.data[,seur_human_integrated_columns_to_keep]

# Rename columns for clarity
colnames(seur.human@meta.data) <- c(
  "nCount_RNA", "nFeature_RNA",
  "percent_mitochondrial", "tcell_lineage",
  "tcell_lineage_tissue", "donor_sex",
  "donor_age_weeks", "donor_id",
  "batch_id", "sequencing_method",
  "tissue", "clusters_integrated_data",
  "trav_trgv_simplified", "trbv_trdv_simplified",
  "trav10_traj18", "trav1_traj33",
  "trav1", "egress_score_old",
  "effector_score_old", "naive_score_old"
)

## /end ####

#_______________________________________
## 2.2. Update T cell lineage names ####

# Update tcell_lineage
table(seur.human@meta.data$tcell_lineage, useNA="ifany")
seur.human@meta.data$tcell_lineage <- case_when(
  seur.human@meta.data$tcell_lineage=="NKT" ~ "iNKT",
  .default=seur.human@meta.data$tcell_lineage
)
seur.human@meta.data$tcell_lineage_tissue <- case_when(
  seur.human@meta.data$tcell_lineage_tissue=="NKT_Thymus" ~ "iNKT_Thymus",
  seur.human@meta.data$tcell_lineage_tissue=="NKT_PBMC" ~ "iNKT_PBMC",
  .default=seur.human@meta.data$tcell_lineage_tissue
)


## /end ####

#___________________________________
## 2.3. Update demographics ####

# Donor sex
table(seur.human@meta.data$donor_sex, useNA="ifany")
seur.human@meta.data$donor_sex <- case_when(
  seur.human@meta.data$donor_sex=="" ~ "NA",
  .default=seur.human@meta.data$donor_sex
)

# More relatable age
table(seur.human@meta.data$donor_age_weeks, useNA="ifany")
seur.human@meta.data$donor_age <- case_when(
  seur.human@meta.data$donor_age_weeks==10   ~ "10wo",
  seur.human@meta.data$donor_age_weeks==16   ~ "16wo",
  seur.human@meta.data$donor_age_weeks==20   ~ "20wo",
  seur.human@meta.data$donor_age_weeks==1252 ~ "24yo",
  seur.human@meta.data$donor_age_weeks==1983 ~ "38yo",
  seur.human@meta.data$donor_age_weeks==2139 ~ "41yo",
  seur.human@meta.data$donor_age_weeks==3548 ~ "68yo",
  .default="NA"
)

table(seur.human@meta.data$donor_id, useNA="ifany")
table(seur.human@meta.data$batch_id, useNA="ifany")
table(seur.human@meta.data$sequencing_method, useNA="ifany")
table(seur.human@meta.data$tissue, useNA="ifany")
## /end ####

#___________________________________
## 2.4. Update clusters ####

# Update integrated clusters
table(seur.human@meta.data$clusters_integrated_data, useNA="ifany")
seur.human@meta.data$clusters_integrated_data <- factor(seur.human@meta.data$clusters_integrated_data, levels=0:17)

# RENAME THYMUS CLUSTERS
DimPlot(seur.thym.cd4, group.by="cell_annot")
DimPlot(seur.thym.cd8, group.by="cell_annot")
DimPlot(seur.thym.nkt, group.by="cell_annot")+scale_color_manual(values=cols_thym_nkt)
DimPlot(seur.thym.mait, group.by="cell_annot")+scale_color_manual(values=cols_thym_mait)
DimPlot(seur.thym.gdt, group.by="cell_annot")+scale_color_manual(values=cols_thym_gdt)
# rename CD4 clusters
seur.thym.cd4@meta.data$clusters_individual_objects <- case_when(
  seur.thym.cd4@meta.data$cell_annot=="thyCD4_ISP"      ~ "CD4_thymus_c0",
  seur.thym.cd4@meta.data$cell_annot=="thyCD4_DPp"      ~ "CD4_thymus_c1",
  seur.thym.cd4@meta.data$cell_annot=="thyCD4_DPq"      ~ "CD4_thymus_c2",
  seur.thym.cd4@meta.data$cell_annot=="thyCD4_ccr9"     ~ "CD4_thymus_c3",
  seur.thym.cd4@meta.data$cell_annot=="thyCD4_ccr7"     ~ "CD4_thymus_c4",
  seur.thym.cd4@meta.data$cell_annot=="thyCD4_Tagonist" ~ "CD4_thymus_c5",
  seur.thym.cd4@meta.data$cell_annot=="thyCD4_Treg"     ~ "CD4_thymus_c6",
  .default="NA"
)
# rename CD8 clusters
seur.thym.cd8@meta.data$clusters_individual_objects <- case_when(
  seur.thym.cd8@meta.data$cell_annot=="thyCD8_DP"     ~ "CD8_thymus_c0",
  seur.thym.cd8@meta.data$cell_annot=="thyCD8_cd8aa1" ~ "CD8_thymus_c1",
  seur.thym.cd8@meta.data$cell_annot=="thyCD8_cd8aa2" ~ "CD8_thymus_c2",
  seur.thym.cd8@meta.data$cell_annot=="thyCD8_ccr9"   ~ "CD8_thymus_c3",
  seur.thym.cd8@meta.data$cell_annot=="thyCD8_ccr7"   ~ "CD8_thymus_c4",
  seur.thym.cd8@meta.data$cell_annot=="thyCD8_idk"    ~ "CD8_thymus_c5",
  .default="NA"
)
# rename iNKT, MAIT, GDT clusters
seur.thym.nkt@meta.data$clusters_individual_objects <- gsub("NKT_", "iNKT_thymus_", seur.thym.nkt@meta.data$cell_annot)
seur.thym.mait@meta.data$clusters_individual_objects <- gsub("MAIT_", "MAIT_thymus_", seur.thym.mait@meta.data$cell_annot)
seur.thym.gdt@meta.data$clusters_individual_objects <- gsub("GDT_", "GDT_thymus_", seur.thym.gdt@meta.data$cell_annot)
# sanity checks
DimPlot(seur.thym.cd4, group.by="clusters_individual_objects")
DimPlot(seur.thym.cd8, group.by="clusters_individual_objects")
DimPlot(seur.thym.nkt, group.by="clusters_individual_objects")
DimPlot(seur.thym.mait, group.by="clusters_individual_objects")
DimPlot(seur.thym.gdt, group.by="clusters_individual_objects")

# RENAME PBMC CLUSTERS
DimPlot(seur.pbmc.cd4, group.by="new_clusters_CD4")+ scale_color_manual(values=cols_pbmc_cd4)
DimPlot(seur.pbmc.cd8, group.by="new_clusters_CD8")+ scale_color_manual(values=cols_pbmc_cd8)
DimPlot(seur.pbmc.nkt, group.by="new_clusters_NKT")+ scale_color_manual(values=cols_pbmc_nkt)
DimPlot(seur.pbmc.mait, group.by="new_clusters_MAIT")+ scale_color_manual(values=cols_pbmc_mait)
DimPlot(seur.pbmc.gdt, group.by="new_clusters_GD")+ scale_color_manual(values=cols_pbmc_gdt)
seur.pbmc.cd4@meta.data$clusters_individual_objects <- paste0("CD4_pbmc_c", seur.pbmc.cd4@meta.data$new_clusters_CD4)
seur.pbmc.cd8@meta.data$clusters_individual_objects <- paste0("CD8_pbmc_c", seur.pbmc.cd8@meta.data$new_clusters_CD8)
seur.pbmc.nkt@meta.data$clusters_individual_objects <- paste0("iNKT_pbmc_c", seur.pbmc.nkt@meta.data$new_clusters_NKT)
seur.pbmc.mait@meta.data$clusters_individual_objects <- paste0("MAIT_pbmc_c", seur.pbmc.mait@meta.data$new_clusters_MAIT)
seur.pbmc.gdt@meta.data$clusters_individual_objects <- paste0("GDT_pbmc_c", seur.pbmc.gdt@meta.data$new_clusters_GD)
DimPlot(seur.pbmc.cd4, group.by="clusters_individual_objects")
DimPlot(seur.pbmc.cd8, group.by="clusters_individual_objects")
DimPlot(seur.pbmc.nkt, group.by="clusters_individual_objects")
DimPlot(seur.pbmc.mait, group.by="clusters_individual_objects")
DimPlot(seur.pbmc.gdt, group.by="clusters_individual_objects")

# EXTRACT CLUSTERS AND ADD THEM TO INTEGRATED OBJECT
df_clusters <- rbind(seur.thym.cd4@meta.data[,c("group.ident", "clusters_individual_objects")],
                     seur.thym.cd8@meta.data[,c("group.ident", "clusters_individual_objects")],
                     seur.thym.nkt@meta.data[,c("group.ident", "clusters_individual_objects")],
                     seur.thym.mait@meta.data[,c("group.ident", "clusters_individual_objects")],
                     seur.thym.gdt@meta.data[,c("group.ident", "clusters_individual_objects")],
                     seur.pbmc.cd4@meta.data[,c("group.ident", "clusters_individual_objects")],
                     seur.pbmc.cd8@meta.data[,c("group.ident", "clusters_individual_objects")],
                     seur.pbmc.nkt@meta.data[,c("group.ident", "clusters_individual_objects")],
                     seur.pbmc.mait@meta.data[,c("group.ident", "clusters_individual_objects")],
                     seur.pbmc.gdt@meta.data[,c("group.ident", "clusters_individual_objects")])
head(df_clusters)
colnames(df_clusters) <- c("tcell_lineage_tissue", "clusters_per_lineage")
df_clusters <- df_clusters[rownames(seur.human@meta.data),]
table(seur.human$tcell_lineage_tissue, useNA="ifany")
table(df_clusters$tcell_lineage_tissue, useNA="ifany")
table(df_clusters$clusters_per_lineage)

## /end ####

#___________________________________
## 2.5. Update TCR usage ####

# Update TCR usage
table(seur.human@meta.data[,c("trav_trgv_simplified", "trav10_traj18")], useNA="ifany")
table(seur.human@meta.data[,c("trav_trgv_simplified", "trav1_traj33")], useNA="ifany")
table(seur.human@meta.data[,c("trav_trgv_simplified", "trav1")], useNA="ifany")
seur.human@meta.data$tcr_usage_simplified <- case_when(
  is.na(seur.human@meta.data$trbv_trdv_simplified) & is.na(seur.human@meta.data$trav_trgv_simplified) ~ "NA",
  seur.human@meta.data$trav10_traj18==1 ~ "TRAV10_TRAJ18",
  seur.human@meta.data$trav1_traj33==1  ~ "TRAV1-2_TRAJ33",
  seur.human@meta.data$trbv_trdv_simplified %in% c("TRDV2*01", "TRDV2*02", "TRDV2*03") & seur.human@meta.data$trav_trgv_simplified %in% c("TRGV9*01", "TRGV9*02") ~ "TRDV2_TRGV9",
  seur.human@meta.data$trbv_trdv_simplified == "TRDV1*01" ~ "TRDV1",
  seur.human@meta.data$trbv_trdv_simplified == "TRDV3*01" ~ "TRDV3",
  .default="tcr_other"
)
table(seur.human@meta.data[,c("trav_trgv_simplified", "tcr_usage_simplified")], useNA="ifany")
table(seur.human@meta.data[,c("trbv_trdv_simplified", "tcr_usage_simplified")], useNA="ifany")
table(seur.human@meta.data[,c("trav10_traj18", "tcr_usage_simplified")], useNA="ifany")
table(seur.human@meta.data[,c("trav1_traj33", "tcr_usage_simplified")], useNA="ifany")
table(seur.human@meta.data[,c("trav1", "tcr_usage_simplified")], useNA="ifany")

seur.human@meta.data[,c("trav_trgv_simplified", "trbv_trdv_simplified", "trav10_traj18", "trav1_traj33", "trav1")] <- NULL
## /end ####

#___________________________________
## 2.6. Update gene signatures ####

# Add egress, naive, effector scores
gene_signatures <- list("effector_new"=c("HOPX", "GZMB", "GZMK", "GNLY", "GZMA", "PRF1", "NKG7", "TBX21",  "KLRD1", "EOMES",
                                         "CCR6", "KLRB1", "RORC", "JUN", "JUNB", "FOS", "IL7R", "ID2", "RORA", "FOSB"),
                        "effector"= c("HOPX", "GZMB", "GZMK", "ZEB2", "NKG7", "GNLY", "TBX21", "EOMES", "TYROBP", "PRF1",
                                      "CCL4", "CCL5", "KLRB1", "GZMH", "GZMA", "KLRD1", "CST7", "KLF6", "CXCR4"),
                        "naive_new"=c("SATB1", "TCF7", "LEF1", "CCR7", "SELL", #"MYC", "EIF3E",
                                  "FOXP1", "KLF2", "SOX4", "ID3", "BACH2"),
                        "naive"=c("SATB1", "TCF7", "LEF1", "CCR7", "SELL", "MYC", "EIF3E",
                                  "SOX4", "ID3", "BACH2"),
                        "egress"=c("KLF2", "CORO1A", "CCR7", "CXCR4", "CXCR6", "FOXO1", "CXCR3", "S1PR1", "S1PR4",
                                   "S100A4", "S100A6", "EMP3"))
seur.human   <- AddModuleScore(seur.human,  name = names(gene_signatures), features=gene_signatures, seed=1)
colnames(seur.human@meta.data)[18:22]  <- names(gene_signatures)

# Offer new naive/effector gene signature to Laurent
max(seur.human@meta.data$naive_score_old) # 4.84
max(seur.human@meta.data$effector_score_old) # 3.88
ggrastr::rasterise(
  SCpubr::do_FeaturePlot(seur.human, features="naive_score_old", border.size=1, pt.size=2, order=T, min.cutoff=0)+
    scale_color_viridis_c(limits=c(0,4.9), option="B"),
  layers="Point", dpi=300) |
  ggrastr::rasterise(
    SCpubr::do_FeaturePlot(seur.human, features="effector_score_old", border.size=1, pt.size=2, order=T, min.cutoff=0)+
      scale_color_viridis_c(limits=c(0,4.9), option="B"),
    layers="Point", dpi=300)
ggsave("./data/human-PBMC/HumanData_34_NaiveEffectorScores/naiveeffectorscores_inseuratobject_082823.pdf", width=14, height=8)

max(seur.human@meta.data$naive) # 2.5
max(seur.human@meta.data$effector) # 1.9
ggrastr::rasterise(
  SCpubr::do_FeaturePlot(seur.human, features="naive", border.size=1, pt.size=2, order=T, min.cutoff=0)+
    scale_color_viridis_c(limits=c(0,2.6), option="B"),
  layers="Point", dpi=300) |
ggrastr::rasterise(
  SCpubr::do_FeaturePlot(seur.human, features="effector", border.size=1, pt.size=2, order=T, min.cutoff=0)+
    scale_color_viridis_c(limits=c(0,2.6), option="B"),
  layers="Point", dpi=300)
ggsave("./data/human-PBMC/HumanData_34_NaiveEffectorScores/naiveeffectorscores_supptable3.pdf", width=14, height=8)

max(seur.human@meta.data$naive_new) # 2.04
max(seur.human@meta.data$effector_new) # 2.14
ggrastr::rasterise(
  SCpubr::do_FeaturePlot(seur.human, features="naive_new", border.size=1, pt.size=2, order=T, min.cutoff=0)+
    scale_color_viridis_c(limits=c(0,2.2), option="B"),
  layers="Point", dpi=300) |
  ggrastr::rasterise(
    SCpubr::do_FeaturePlot(seur.human, features="effector_new", border.size=1, pt.size=2, order=T, min.cutoff=0)+
      scale_color_viridis_c(limits=c(0,2.2), option="B"),
    layers="Point", dpi=300)
ggsave("./data/human-PBMC/HumanData_34_NaiveEffectorScores/naiveeffectorscores_proposition.pdf", width=14, height=8)


# Egress score
max(seur.human@meta.data$egress_score_old) # 1.36
max(seur.human@meta.data$egress) # 1.35
ggrastr::rasterise(
  SCpubr::do_FeaturePlot(seur.human, features="egress_score_old", split.by="tissue", border.size=1, pt.size=2, order=T, min.cutoff=0,
                         use_viridis = T, viridis.palette="B"),
  layers="Point", dpi=300)
ggsave("./data/human-PBMC/HumanData_34_NaiveEffectorScores/egressscore_inseuratobject_082823_2.pdf", width=14, height=8)
ggrastr::rasterise(
  SCpubr::do_FeaturePlot(seur.human, features="egress", split.by="tissue", border.size=1, pt.size=2, order=T, min.cutoff=0,
                         use_viridis = T, viridis.palette="B"),
  layers="Point", dpi=100)
ggsave("./data/human-PBMC/HumanData_34_NaiveEffectorScores/egressscore_supptable3_2.pdf", width=14, height=8)

# remove naive_score_old, egress_score_old, 



## /end ####



# *******************************************
# 3. PREPARE HUMAN THYMIC SEURAT OBJECTS ####
# *******************************************

#___________________________
## 3.1. P ####

## /end ####


#___________________________
## 3.2. Second analysis ####

## /end ####


#___________________________
## 3.3. Third analysis ####

## /end ####


#___________________________
## 3.4. Fourth analysis ####

## /end ####




# *******************************************
# 4. PREPARE MOUSE THYMIC SEURAT OBJECTS ####
# *******************************************






