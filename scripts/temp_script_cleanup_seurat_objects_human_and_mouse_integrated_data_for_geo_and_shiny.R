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
seur.pbmc.cd8  <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.CD8_noMAIT_08_17_23.RDS")
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
table(seur.human@meta.data$tcell_lineage_tissue, useNA="ifany")

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
seur.human@meta.data <- seur.human@meta.data %>% relocate(donor_age, .after=donor_age_weeks)
table(seur.human@meta.data$donor_age, useNA="ifany")

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
table(rownames(df_clusters) == rownames(seur.human@meta.data), useNA="ifany") # all TRUE
table(seur.human$tcell_lineage_tissue, useNA="ifany")
table(df_clusters$tcell_lineage_tissue, useNA="ifany")

seur.human@meta.data$clusters_per_lineage <- df_clusters$clusters_per_lineage
seur.human@meta.data <- seur.human@meta.data %>% relocate(clusters_per_lineage, .after=clusters_integrated_data)
colnames(seur.human@meta.data)
table(seur.human@meta.data$clusters_per_lineage, useNA="ifany")
## /end ####

#___________________________________
## 2.5. Update TCR usage ####

# Update TCR usage
# table(seur.human@meta.data[,c("trav_trgv_simplified", "trav10_traj18")], useNA="ifany")
# table(seur.human@meta.data[,c("trav_trgv_simplified", "trav1_traj33")], useNA="ifany")
# table(seur.human@meta.data[,c("trav_trgv_simplified", "trav1")], useNA="ifany")
seur.human@meta.data$tcr_usage_simplified <- case_when(
  is.na(seur.human@meta.data$trbv_trdv_simplified) & is.na(seur.human@meta.data$trav_trgv_simplified) ~ "NA",
  seur.human@meta.data$trav10_traj18==1 ~ "TRAV10_TRAJ18",
  seur.human@meta.data$trav1_traj33==1  ~ "TRAV1-2_TRAJ33",
  seur.human@meta.data$trbv_trdv_simplified %in% c("TRDV2*01", "TRDV2*02", "TRDV2*03") & seur.human@meta.data$trav_trgv_simplified %in% c("TRGV9*01", "TRGV9*02") ~ "TRDV2_TRGV9",
  seur.human@meta.data$trbv_trdv_simplified == "TRDV1*01" ~ "TRDV1",
  seur.human@meta.data$trbv_trdv_simplified == "TRDV3*01" ~ "TRDV3",
  .default="tcr_other"
)
# table(seur.human@meta.data[,c("trav_trgv_simplified", "tcr_usage_simplified")], useNA="ifany")
# table(seur.human@meta.data[,c("trbv_trdv_simplified", "tcr_usage_simplified")], useNA="ifany")
# table(seur.human@meta.data[,c("trav10_traj18", "tcr_usage_simplified")], useNA="ifany")
# table(seur.human@meta.data[,c("trav1_traj33", "tcr_usage_simplified")], useNA="ifany")
# table(seur.human@meta.data[,c("trav1", "tcr_usage_simplified")], useNA="ifany")

seur.human@meta.data[,c("trav_trgv_simplified", "trbv_trdv_simplified", "trav10_traj18", "trav1_traj33", "trav1")] <- NULL

colnames(seur.human@meta.data)
## /end ####

#___________________________________
## 2.6. Update gene signatures ####

# Add egress, naive, effector scores
gene_signatures <- list("score_effector"=c("HOPX", "GZMB", "GZMK", "GNLY", "GZMA", "PRF1", "NKG7", "TBX21",  "KLRD1", "EOMES",
                                           "CCR6", "KLRB1", "RORC", "JUN", "JUNB", "FOS", "IL7R", "ID2", "RORA", "FOSB"),
                        # "effector"= c("HOPX", "GZMB", "GZMK", "ZEB2", "NKG7", "GNLY", "TBX21", "EOMES", "TYROBP", "PRF1",
                        #               "CCL4", "CCL5", "KLRB1", "GZMH", "GZMA", "KLRD1", "CST7", "KLF6", "CXCR4"),
                        "score_naive"=c("SATB1", "TCF7", "LEF1", "CCR7", "SELL", "FOXP1", "KLF2", "SOX4", "ID3", "BACH2"),
                        # "naive"=c("SATB1", "TCF7", "LEF1", "CCR7", "SELL", "MYC", "EIF3E",
                        #           "SOX4", "ID3", "BACH2"),
                        "score_egress"=c("KLF2", "CORO1A", "CCR7", "CXCR4", "CXCR6", "FOXO1", "CXCR3", "S1PR1", "S1PR4",
                                         "S100A4", "S100A6", "EMP3"),
                        "score_cd8aa"=c("GNG4", "CD8A", "NUCB2", "LEF1","PRKCH","PTPN7","SH2D1A","ELOVL5","TOX","AQP3",
                                        "BCL6","ITGA4","MYB","RTKN2","IKZF2","DUSP2","MINDY2","HIVEP3"))
seur.human   <- AddModuleScore(seur.human,  name = names(gene_signatures), features=gene_signatures, seed=1)
colnames(seur.human@meta.data)[19:22]  <- names(gene_signatures)

# Offer new naive/effector gene signature to Laurent
# max(seur.human@meta.data$naive_score_old) # 4.84
# max(seur.human@meta.data$effector_score_old) # 3.88
# ggrastr::rasterise(
#   SCpubr::do_FeaturePlot(seur.human, features="naive_score_old", border.size=1, pt.size=2, order=T, min.cutoff=0)+
#     scale_color_viridis_c(limits=c(0,4.9), option="B"),
#   layers="Point", dpi=300) |
#   ggrastr::rasterise(
#     SCpubr::do_FeaturePlot(seur.human, features="effector_score_old", border.size=1, pt.size=2, order=T, min.cutoff=0)+
#       scale_color_viridis_c(limits=c(0,4.9), option="B"),
#     layers="Point", dpi=300)
# ggsave("./data/human-PBMC/HumanData_34_NaiveEffectorScores/naiveeffectorscores_inseuratobject_082823.pdf", width=14, height=8)
# 
# max(seur.human@meta.data$naive) # 2.5
# max(seur.human@meta.data$effector) # 1.9
# ggrastr::rasterise(
#   SCpubr::do_FeaturePlot(seur.human, features="naive", border.size=1, pt.size=2, order=T, min.cutoff=0)+
#     scale_color_viridis_c(limits=c(0,2.6), option="B"),
#   layers="Point", dpi=300) |
# ggrastr::rasterise(
#   SCpubr::do_FeaturePlot(seur.human, features="effector", border.size=1, pt.size=2, order=T, min.cutoff=0)+
#     scale_color_viridis_c(limits=c(0,2.6), option="B"),
#   layers="Point", dpi=300)
# ggsave("./data/human-PBMC/HumanData_34_NaiveEffectorScores/naiveeffectorscores_supptable3.pdf", width=14, height=8)
# 
# max(seur.human@meta.data$naive_new) # 2.04
# max(seur.human@meta.data$effector_new) # 2.14
# ggrastr::rasterise(
#   SCpubr::do_FeaturePlot(seur.human, features="naive_new", border.size=1, pt.size=2, order=T, min.cutoff=0)+
#     scale_color_viridis_c(limits=c(0,2.2), option="B"),
#   layers="Point", dpi=300) |
#   ggrastr::rasterise(
#     SCpubr::do_FeaturePlot(seur.human, features="effector_new", border.size=1, pt.size=2, order=T, min.cutoff=0)+
#       scale_color_viridis_c(limits=c(0,2.2), option="B"),
#     layers="Point", dpi=300)
# ggsave("./data/human-PBMC/HumanData_34_NaiveEffectorScores/naiveeffectorscores_proposition.pdf", width=14, height=8)


# Egress score
# max(seur.human@meta.data$egress_score_old) # 1.36
# max(seur.human@meta.data$egress) # 1.35
# ggrastr::rasterise(
#   SCpubr::do_FeaturePlot(seur.human, features="egress_score_old", split.by="tissue", border.size=1, pt.size=2, order=T, min.cutoff=0,
#                          use_viridis = T, viridis.palette="B"),
#   layers="Point", dpi=300)
# ggsave("./data/human-PBMC/HumanData_34_NaiveEffectorScores/egressscore_inseuratobject_082823_2.pdf", width=14, height=8)
# ggrastr::rasterise(
#   SCpubr::do_FeaturePlot(seur.human, features="score_egress", split.by="tissue", border.size=1, pt.size=2, order=T, min.cutoff=0,
#                          use_viridis = T, viridis.palette="B"),
#   layers="Point", dpi=100)
# ggsave("./data/human-PBMC/HumanData_34_NaiveEffectorScores/egressscore_supptable3_2.pdf", width=14, height=8, dpi=100)

# CD8aa score
# ggrastr::rasterise(
#   SCpubr::do_FeaturePlot(seur.human, features="score_cd8aa", border.size=1, pt.size=2, order=T, min.cutoff=0,
#                          use_viridis = T, viridis.palette="B"),
#   layers="Point", dpi=300)
# ggsave("./data/human-PBMC/HumanData_34_NaiveEffectorScores/cd8aascore_supptable3.pdf", width=7, height=8)

# remove naive_score_old, egress_score_old, effector_score_old
seur.human@meta.data[,c("naive_score_old", "egress_score_old", "effector_score_old")] <- NULL
colnames(seur.human@meta.data)


# Plot for fig1
# SCpubr::do_FeaturePlot(seur.human, features="score_naive", border.size=1, pt.size=3, order=T, min.cutoff=0, max.cutoff=2,
#                        raster=T, raster.dpi=2048, use_viridis = T, viridis.palette = "B") +
#   scale_color_viridis_c(limits=c(0,2), option="B") |
# SCpubr::do_FeaturePlot(seur.human, features="score_effector", border.size=1, pt.size=3, order=T, min.cutoff=0, #max.cutoff=2,
#                        raster=T, raster.dpi=2048, use_viridis = T, viridis.palette = "B") +
#   scale_color_viridis_c(limits=c(0,2), option="B")
# ggsave("./data/human-PBMC/HumanData_34_NaiveEffectorScores/naiveeffectorscores_proposition_laurentparameters_ptsize3_maxcutoff2.pdf", width=14, height=8)
# 
# SCpubr::do_FeaturePlot(seur.human, features="score_egress", split.by="tissue", border.size=1, pt.size=5, order=T, min.cutoff=0, max.cutoff=1,
#                        raster=T, raster.dpi=2048, use_viridis = T, viridis.palette = "B")
# ggsave("./data/human-PBMC/HumanData_34_NaiveEffectorScores/egresscore_supptable3_laurentparameters_ptsize5_maxcutoff1.pdf", width=14, height=8)


## /end ####


#____________________________
## 2.7. Update GEP usage ####

# Import GEP usage matrix
gep_usage <- read.table("./data/human-PBMC/HumanData_20_RibbonPlotCellStateToID/cNMF_output/non_imputed_cNMF_allcells.usages.k_12.dt_0_02.consensus.txt", header=T)
colnames(gep_usage) <- paste0("GEP", c(2,5,3,1,4,12,6,7,8,10,9,11), "_usage")
gep_usage$GEP12_usage <- NULL # remove GEP driven by batch effect
gep_usage <- gep_usage[,paste0("GEP", 1:11, "_usage")] # reorder columns
# head(gep_usage)

# Add gep assign
gep_usage$GEP_with_max_usage <- gsub("_.*", "", colnames(gep_usage)[apply(gep_usage, 1, which.max)])
table(gep_usage$GEP_with_max_usage, useNA="ifany")

# Add to seur object
table(rownames(gep_usage)==rownames(seur.human@meta.data), useNA="ifany")
seur.human@meta.data <- cbind(seur.human@meta.data, gep_usage)

# Sanity check
# SCpubr::do_FeaturePlot(seur.human, features=paste0("GEP", 1:11, "_usage"), border.size=1, pt.size=1, order=T, use_viridis = T, viridis.palette="B",
#                        ncol=2)

colnames(seur.human@meta.data)

## /end ####


#____________________________
## 2.8. Rename dimension reduction of interest and save ####
names(seur.human@reductions)[4] <- "umap_integrated"

head(as.data.frame(seur.human@meta.data))

# Save seurat
saveRDS(seur.human, "./data/clean_data/seurat_human_integrated_object_23_12_01.rds")

# Save metadata
mtdata <- as.data.frame(seur.human@meta.data) %>% rownames_to_column("cellid")
write_csv(mtdata, "./data/clean_data/metadata_human_integrated_object_23_12_01.csv")

# SAVE COUNTS
seur.human@assays$RNA@counts[1:5,1:5]
# mart.hu <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
# mart.df <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'ensembl_gene_id_version', 'external_gene_name', 'external_synonym',
#                                          'hgnc_symbol'#,
#                                          # 'uniprot_gn_symbol'#, 'uniprot_gn_id',
#                                          ),
#                           filters = 'external_gene_name',
#                           values = rownames(seur.human),
#                           mart = mart.hu)
# dim(mart.df)
# table(rownames(seur.human) %in% mart.df$external_gene_name, useNA="ifany")

# external_gene_name: 3192 missing (genes_missing_externalgenename)
# external_synonym: 16,429 missing, only 775 present (genes_found_synonyms)
# hgnc_symbol: 3199 missing (genes_missing_hgnc)
# uniprot_gn_symbol: 3627 missing (genes_missing_uniprot)

mat <- seur.human@assays$RNA@counts
h5_target <- "./data/clean_data/rawcounts_human_integrated_object_23_12_01.h5"
ref_name <- "hu_rawcounts_postqc"

rhdf5::h5createFile(h5_target) # create target file
rhdf5::h5createGroup(h5_target, ref_name) # create data group
rhdf5::h5write(colnames(mat), h5_target, paste0("/",ref_name,"/barcodes")) # store cell IDs
rhdf5::h5write(rownames(mat), h5_target, paste0("/",ref_name,"/gene_names")) # store gene names
rhdf5::h5write(rownames(mat), h5_target, paste0("/",ref_name,"/genes")) # store gene names (need to add ENSEMBL IDs later)

rhdf5::h5write(dim(mat), h5_target, paste0("/",ref_name,"/shape")) # store dimensions as shape
# store values from mat@x as data
rhdf5::h5createDataset(h5_target, paste0("/",ref_name,"/data"), dims = length(mat@x), storage.mode = "integer", chunk = 1000, level = 4)
rhdf5::h5write(mat@x, h5_target,  paste0("/",ref_name,"/data"))
# store row indices from mat@i as indices
rhdf5::h5createDataset(h5_target, paste0("/",ref_name,"/indices"), dims = length(mat@i), storage.mode = "integer", chunk = 1000, level = 4)
rhdf5::h5write(mat@i, h5_target,  paste0("/",ref_name,"/indices"))
# store column pointers from mat@p as indptr
rhdf5::h5write(mat@p, h5_target, paste0("/",ref_name,"/indptr"))


## /end ####




# *******************************************
# 3. PREPARE MOUSE THYMIC SEURAT OBJECTS ####
# *******************************************

#___________________________
## 3.1. Keep columns of interest ####

seur_mouse_integrated_columns_to_keep <- c(
  "orig.ident", "nCount_RNA", "nFeature_RNA",
  "percent.mt", "TRAV",
  "TRBV", "TRAJ", "TRBJ",
  "NKT_TCR", "MAIT_TCR",
  "New.ident",
  "Cell.type",
  "New_clusters"
)
seur.ms@meta.data <- seur.ms@meta.data[,seur_mouse_integrated_columns_to_keep]

# Rename columns for clarity
colnames(seur.ms@meta.data) <- c(
  "ms_replicate", "nCount_RNA", "nFeature_RNA", "percent_mitochondrial", 
  "TRAV", "TRBV", "TRAJ", "TRBJ",
  "iNKT_TCR", "MAIT_TCR", "study", "tcell_lineage",
  "clusters_integrated_data"
)
colnames(seur.ms@meta.data)

## /end ####


#___________________________
## 3.2. Update study & ms_replicate ####
seur.ms@meta.data$study <- case_when(
  seur.ms@meta.data$study=="Chandra_MAIT" ~ "Chandra_MAIT",
  seur.ms@meta.data$study=="KMB_NKT"      ~ "MaasBauer_iNKT",
  seur.ms@meta.data$study=="Koay_MAIT"    ~ "Koay_MAIT",
  seur.ms@meta.data$study=="Krovi_NKT"    ~ "Krovi_iNKT",
  seur.ms@meta.data$study=="Lee_Gd"       ~ "Lee_GD",
  seur.ms@meta.data$study=="Lee_MAIT"     ~ "Lee_MAIT",
  seur.ms@meta.data$study=="Lee_NKT"      ~ "Lee_iNKT",
  seur.ms@meta.data$study=="Legoux_MAIT"  ~ "Legoux_MAIT",
  seur.ms@meta.data$study=="Li_Gd"        ~ "Li_GD",
  seur.ms@meta.data$study=="Paget_NKT"    ~ "Baranek_iNKT",
  seur.ms@meta.data$study=="SAP_MAIT"     ~ "Legoux_MAIT_SAP_KO",
  seur.ms@meta.data$study=="Stage0_NKT"   ~ "Wang_iNKT_stage0",
  seur.ms@meta.data$study=="Stage1_NKT"   ~ "Wang_iNKT_stage1",
  seur.ms@meta.data$study=="Stage2_NKT"   ~ "Wang_iNKT_stage2",
  seur.ms@meta.data$study=="Stage3_NKT"   ~ "Wang_iNKT_stage3"
)

seur.ms@meta.data$ms_replicate <- case_when(
  seur.ms@meta.data$ms_replicate=="shCh01x6"    ~ "MAIT_C56BL6_Chandra",
  seur.ms@meta.data$ms_replicate=="FVB_1NKT"    ~ "iNKT_FVBN_MaasBauer_rep1",
  seur.ms@meta.data$ms_replicate=="FVB_2NKT"    ~ "iNKT_FVBN_MaasBauer_rep2",
  seur.ms@meta.data$ms_replicate=="Thymus"      ~ "MAIT_C57BL6_Koay",
  seur.ms@meta.data$ms_replicate=="B6NKT"       ~ "iNKT_C57BL6_Krovi_rep1",
  seur.ms@meta.data$ms_replicate=="B6_2NKT"     ~ "iNKT_C57BL6_Krovi_rep2",
  seur.ms@meta.data$ms_replicate=="Lee"         ~ "GD_BALBc_Lee",
  seur.ms@meta.data$ms_replicate=="Lee_MAIT"    ~ "MAIT_BALBc_Lee",
  seur.ms@meta.data$ms_replicate=="Lee_NKT"     ~ "iNKT_BALBc_Lee",
  seur.ms@meta.data$ms_replicate=="Legoux_MAIT" ~ "MAIT_B6MAITCast_Legoux",
  seur.ms@meta.data$ms_replicate=="Li_Gd"       ~ "GD_C57BL6_Li",
  seur.ms@meta.data$ms_replicate=="Paget_NKT"   ~ "iNKT_C57BL6_Baranek",
  seur.ms@meta.data$ms_replicate=="SAP_MAIT"    ~ "MAIT_Sh2d1aKO_Legoux",
  seur.ms@meta.data$ms_replicate=="Stage0_1NKT" ~ "iNKT_stage0_Va14Tg_Wang_rep1",
  seur.ms@meta.data$ms_replicate=="Stage0_2NKT" ~ "iNKT_stage0_Va14Tg_Wang_rep2",
  seur.ms@meta.data$ms_replicate=="Stage1_1NKT" ~ "iNKT_stage1_C57BL6_Wang_rep1",
  seur.ms@meta.data$ms_replicate=="Stage1_2NKT" ~ "iNKT_stage1_C57BL6_Wang_rep2",
  seur.ms@meta.data$ms_replicate=="Stage2_1NKT" ~ "iNKT_stage2_C57BL6_Wang_rep1",
  seur.ms@meta.data$ms_replicate=="Stage2_2NKT" ~ "iNKT_stage2_C57BL6_Wang_rep2",
  seur.ms@meta.data$ms_replicate=="Stage3_1NKT" ~ "iNKT_stage3_C57BL6_Wang_rep1",
  seur.ms@meta.data$ms_replicate=="Stage3_2NKT" ~ "iNKT_stage3_C57BL6_Wang_rep2"
)

seur.ms@meta.data <- seur.ms@meta.data %>% relocate(ms_replicate, .after=study)
seur.ms@meta.data <- seur.ms@meta.data %>% relocate(tcell_lineage, .after=percent_mitochondrial)

table(seur.ms@meta.data[,c("ms_replicate", "tcell_lineage")], useNA="ifany")
table(seur.ms@meta.data[,c("study", "tcell_lineage")], useNA="ifany")
## /end ####


#___________________________
## 3.3. Add annotation ####
seur.ms@meta.data$clusters_annotation <- case_when(
  seur.ms@meta.data$clusters_integrated_data==0  ~ "0_immature",
  seur.ms@meta.data$clusters_integrated_data==1  ~ "1_GD_Gzma",
  seur.ms@meta.data$clusters_integrated_data==2  ~ "2_GD_Gzma",
  seur.ms@meta.data$clusters_integrated_data==3  ~ "signaling",
  seur.ms@meta.data$clusters_integrated_data==4  ~ "signaling",
  seur.ms@meta.data$clusters_integrated_data==5  ~ "cycling",
  seur.ms@meta.data$clusters_integrated_data==6  ~ "transition",
  seur.ms@meta.data$clusters_integrated_data==7  ~ "typeII",
  seur.ms@meta.data$clusters_integrated_data==8  ~ "typeI",
  seur.ms@meta.data$clusters_integrated_data==9  ~ "typeI",
  seur.ms@meta.data$clusters_integrated_data==10 ~ "typeIII",
  seur.ms@meta.data$clusters_integrated_data==11 ~ "typeIII",
  seur.ms@meta.data$clusters_integrated_data==12 ~ "12"
)
seur.ms@meta.data$clusters_annotation <- factor(seur.ms@meta.data$clusters_annotation, levels=c(
  "0_immature",
  "1_GD_Gzma",
  "2_GD_Gzma",
  "signaling",
  "cycling",
  "transition",
  "typeI",
  "typeII",
  "typeIII",
  "12"
))
table(seur.ms@meta.data[,c("clusters_integrated_data", "clusters_annotation")], useNA="ifany")
## /end ####


#___________________________
## 3.4. Update TCR info & save ####

# rename lineage
seur.ms@meta.data$tcell_lineage <- case_when(
  seur.ms@meta.data$tcell_lineage=="NKT"  ~ "iNKT",
  .default=seur.ms@meta.data$tcell_lineage
)
table(seur.ms@meta.data$tcell_lineage, useNA="ifany")

# rename TCR info
seur.ms@meta.data$inkt_mait_tcr <- case_when(
  seur.ms@meta.data$iNKT_TCR==1 & seur.ms@meta.data$TRAV %in% c("TRAV11", "TRAV11D") & seur.ms@meta.data$TRAJ=="TRAJ18"  ~ "TRAV11-TRAJ18 (iNKT)",
  seur.ms@meta.data$MAIT_TCR==1 & seur.ms@meta.data$TRAV=="TRAV1" & seur.ms@meta.data$TRAJ=="TRAJ33" ~ "TRAV1-TRAJ33 (MAIT)",
  is.na(seur.ms@meta.data$iNKT_TCR) | is.na(seur.ms@meta.data$MAIT_TCR) ~ "NA",
  .default=NA
)
table(seur.ms@meta.data[,c("iNKT_TCR", "inkt_mait_tcr")], useNA="ifany")
table(seur.ms@meta.data[,c("MAIT_TCR", "inkt_mait_tcr")], useNA="ifany")
table(seur.ms@meta.data[,c("TRAV", "inkt_mait_tcr")], useNA="ifany")
table(seur.ms@meta.data[,c("TRAJ", "inkt_mait_tcr")], useNA="ifany")

seur.ms@meta.data$iNKT_TCR <- NULL
seur.ms@meta.data$MAIT_TCR <- NULL
seur.ms@meta.data <- seur.ms@meta.data %>% relocate(inkt_mait_tcr, .after="tcell_lineage")

# reorder columns a bit
seur.ms@meta.data <- seur.ms@meta.data %>% relocate(study, .after="percent_mitochondrial")
seur.ms@meta.data <- seur.ms@meta.data %>% relocate(ms_replicate, .after="study")

# Save seurat
saveRDS(seur.ms, "./data/clean_data/seurat_mouse_integrated_object_23_12_02.rds")

## /end ####