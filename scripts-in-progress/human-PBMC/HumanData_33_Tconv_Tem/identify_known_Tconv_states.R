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
library(ggalluvial)
source("./scripts-final/colors_universal.R")

# Import integrated data
seur.human <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_08_28_23.RDS")
DimPlot(seur, group.by = "new_clusters", label=T, repel=T, reduction="UMAP_50")+
  scale_color_manual(values=cols_integrated)

# Import PBMC data
seur.cd4  <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.CD4_03_16_23.RDS")
seur.cd8  <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.CD8_03_16_23.RDS")
# seur.gdt  <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.GD_03_16_23.RDS")
# seur.mait <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.MAIT_03_16_23.RDS")
# seur.nkt  <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.NKT_03_16_23.RDS")

# Import GEP usage
gep_usage <- read.table("./data/human-PBMC/HumanData_20_RibbonPlotCellStateToID/cNMF_output/non_imputed_cNMF_allcells.usages.k_12.dt_0_02.consensus.txt", header=T)
colnames(gep_usage) <- paste0("gep", c(2,5,3,1,4,12,6,7,8,10,9,11), "_usage")
head(gep_usage)

# Import gene lists
genes.gep     <- read.csv("./data/human-PBMC/HumanData_27_ComparisonGEPsGarner/limited_nonimputed_genes_per_gep_post_rank_threshold_k12.csv", row.names=1)
genes.cano    <- readxl::read_excel("./data/human-PBMC/HumanData_17_GEPsOnCanogamezData/canogamez_supp.xlsx", skip=3, sheet=1)
genesCD8.rose <- readxl::read_excel("./data/human-PBMC/HumanData_22_CompareGeneLists/rose_supp_clusmodulesCD8.xlsx", sheet=1)
genesCD4.rose <- readxl::read_excel("./data/human-PBMC/HumanData_22_CompareGeneLists/rose_supp_clusmodulesCD4.xlsx", sheet=1)
genes.poon    <- readxl::read_excel("./data/human-PBMC/HumanData_22_CompareGeneLists/poon_supp.xlsx", sheet=6)[,-1]
genes.terekhova <- readxl::read_excel("./data/human-PBMC/HumanData_33_Tconv_Tem/terekhova_supptable.xlsx", sheet=5)




# ************************************
# 2. PREPARE AND SCORE GENE LISTS ####
# ************************************

genes.long <- list()

# ***************
## 2.1. GEPs ####
head(genes.gep)
gapindf <- gather(genes.gep, key=geneprogram, value=gene, colnames(genes.gep)) %>%
  filter(!is.na(gene)) %>%
  mutate(geneprogram=sub("_", "", geneprogram))
# head(gapindf)
# table(gapindf$geneprogram, useNA="ifany")
# table(is.na(gapindf$gene))
genes.long[["gapin"]] <- gapindf
## end geps ####


# **************************
## 2.2. Cano-Gamez data ####
head(genes.cano)
# take a look at the different clusters/gene programs available
# table(genes.cano$cluster, useNA="ifany")
# rename the clusters a bit
canodf <- genes.cano %>%
  mutate(cluster=gsub(" ", "", cluster),
         cluster=paste0("CanoGamez_CD4_", cluster)) %>%
  # keep only significant genes & with minimum log2FC
  filter(p_val_adj < 0.05 & LFC > 0.25) %>% # lowest log2FC is 0.25
  # convert to long format
  select(cluster, gene) %>%
  rename(geneprogram=cluster)
head(canodf)
genes.long[["canogamez"]] <- canodf
# table(is.na(genes.long$canogamez$gene)) # sanity check
# table(genes.long$canogamez$gene %in% rownames(seur.human))
## end cano-gamez ####


# **************************
## 2.3. Rose data ####
head(genesCD4.rose)
head(genesCD8.rose)

roseCD4 <- genesCD4.rose %>%
  select(SYMBOL, k.5.cluster) %>%
  filter(!is.na(SYMBOL)) %>% # the few NA symbols are miRNA or uncharacterized loci
  rename(gene=SYMBOL, geneprogram=k.5.cluster) %>%
  mutate(geneprogram=ifelse(geneprogram==1, "Rose_CD4_modul1_Tem_cm",
                            ifelse(geneprogram==2, "Rose_CD4_modul2_Tcm_em",
                                   ifelse(geneprogram==3, "Rose_CD4_modul3_Tem_cm",
                                          ifelse(geneprogram==4, "Rose_CD4_modul4_Tem",
                                                 ifelse(geneprogram==5, "Rose_CD4_modul5_Tnaive", "?"))))))

roseCD8 <- genesCD8.rose %>%
  select(SYMBOL, k.5.cluster) %>%
  filter(!is.na(SYMBOL)) %>% # the few NA symbols are miRNA or uncharacterized loci
  rename(gene=SYMBOL, geneprogram=k.5.cluster) %>%
  mutate(geneprogram=ifelse(geneprogram==1, "Rose_CD8_modul1_Temra",
                            ifelse(geneprogram==2, "Rose_CD8_modul2_Tcm_em",
                                   ifelse(geneprogram==3, "Rose_CD8_modul3_Tem_emra",
                                          ifelse(geneprogram==4, "Rose_CD8_modul4_Tnaive",
                                                 ifelse(geneprogram==5, "Rose_CD8_modul5_Tnaive_cm", "?"))))))

head(rbind(roseCD4, roseCD8))
genes.long[["rose"]] <- rbind(roseCD4, roseCD8)
# table(is.na(genes.long$rose$gene)) # sanity check
# table(genes.long$rose$gene %in% rownames(seur.human))
# table(genes.long$rose$geneprogram, useNA="ifany")
## end rose ####


# ********************
## 2.4. Poon data ####
head(genes.poon)
colnames(genes.poon) <- gsub(" ", "", colnames(genes.poon))

# Transform to long format gene program by gene program (each gene program has 3 columns: genes symbols, padj, logFC)
poondf <- data.frame()
for (i in seq(1,ncol(genes.poon), 3)){
  subdf <- genes.poon[,i:(i+2)]
  geneprogram <- sub("_.*", "",colnames(subdf)[1])
  colnames(subdf) <- c("gene", "padj", "logFC")
  subdf <- subdf %>%
    # OPTION 1: full list of Poon genes
    # filter(padj<0.01 & logFC>0.25)
    # filter(padj<0.001 & logFC>0.25)
    # OPTION 2: short list of Poon genes
    slice_head(n=500) # take 500 top genes
  cat(geneprogram, "\n")
  cat("# genes:", nrow(subdf), "\n")
  cat("min log2FC:", min(subdf$logFC), "\n\n")
  subdf$geneprogram <- paste0("Poon_", geneprogram) # add column with the name of the geneprogram
  poondf <- rbind(poondf, subdf)
}
# head(poondf)
# table(poondf$geneprogram, useNA="ifany")
# table(is.na(poondf$gene))

genes.long[["poon"]] <- poondf[,c("geneprogram", "gene")]
## end poon ####

# ********************
## 2.5. Terekhova data ####
head(genes.terekhova)
terekhovadf <- genes.terekhova %>%
  pivot_longer(cols=everything(), names_to="geneprogram", values_to = "gene") %>%
  mutate(geneprogram=paste0("Terekhova_", geneprogram))
genes.long[["terekhova"]] <- terekhovadf
## end terekhova ####


# *********************************
## 2.5. Get it all in one list ####
longdf <- bind_rows(genes.long, .id="dataset") %>%
  filter(gene %in% rownames(seur.human)) %>%
  filter(geneprogram != "GEP12")
# table(longdf$dataset, useNA="ifany")
# table(longdf$geneprogram, useNA="ifany")
# canogamez     gapin      poon      rose  terekhova
#   341         3,787     3,885     1,542      1,110

# make a list of CD4 gene signatures
cd4_genesig <- c("CanoGamez_CD4_nTreg", "CanoGamez_CD4_TCM", "CanoGamez_CD4_TEM", "CanoGamez_CD4_TEMRA", "CanoGamez_CD4_Tnaive",
                 "GEP3", "GEP4", "GEP5", "GEP6",
                 "Poon_CD4Naive", "Poon_CD4TCM/TFH", "Poon_CD4Treg", "Poon_CD4TRM",
                 "Rose_CD4_modul1_Tem_cm", "Rose_CD4_modul2_Tcm_em", "Rose_CD4_modul3_Tem_cm", "Rose_CD4_modul4_Tem", "Rose_CD4_modul5_Tnaive")
genesiglist_cd4 <- list()
for(geneprog in cd4_genesig){
  elementname <- gsub("/", "_", geneprog)
  genesiglist_cd4[[elementname]] <- longdf[longdf$geneprogram==geneprog, "gene"]
}
lapply(genesiglist_cd4, function(x) length(x))


# make a list of CD8 gene signatures
cd8_genesig <- c("GEP3", "GEP4", "GEP5", "GEP6",
                 "Poon_CD8MAIT", "Poon_CD8Naive", "Poon_CD8TEM/TEMRA", "Poon_CD8TRM",
                 "Rose_CD8_modul1_Temra", "Rose_CD8_modul2_Tcm_em", "Rose_CD8_modul3_Tem_emra", "Rose_CD8_modul4_Tnaive", "Rose_CD8_modul5_Tnaive_cm",
                 grep("Terekhova", unique(longdf$geneprogram), value=T))
genesiglist_cd8 <- list()
for(geneprog in cd8_genesig){
  elementname <- gsub("/", "_", geneprog)
  genesiglist_cd8[[elementname]] <- longdf[longdf$geneprogram==geneprog, "gene"]
}
lapply(genesiglist_cd8, function(x) length(x))
## end 2.5 ####


# *********************************************
## 2.6. Score gene lists on seurat objects ####
seur.cd4 <- AddModuleScore(seur.cd4, features=genesiglist_cd4, name=names(genesiglist_cd4), seed=1)
seur.cd8 <- AddModuleScore(seur.cd8, features=genesiglist_cd8, name=names(genesiglist_cd8), seed=1)
colnames(seur.cd4@meta.data)[42:59] <- names(genesiglist_cd4)
colnames(seur.cd8@meta.data)[42:66] <- names(genesiglist_cd8)
# FeaturePlot(seur.cd4, features=c("GEP3", "GEP4", "GEP5", "GEP6"), order=T)
# FeaturePlot(seur.cd8, features=c("GEP3", "GEP4", "GEP5", "GEP6"), order=T)

genesiglist_both <- c(genesiglist_cd4, genesiglist_cd8[5:25])
seur.human <- AddModuleScore(seur.human, features=genesiglist_both, name=names(genesiglist_both), seed=1)
colnames(seur.human@meta.data)[50:88] <- names(genesiglist_both)
# FeaturePlot(seur.human, features=c("GEP_Scores_3", "GEP3", "GEP_Scores_5", "GEP5"), order=T, reduction="UMAP_50")
## end 2.6 ####

# *********************************************
## 2.7. Tem gene lists ####

# CD8
FeaturePlot(seur.cd8, features=c("Poon_CD8Naive", "Rose_CD8_modul4_Tnaive", "Terekhova_CD8_Tnaive",
                                 "Poon_CD8TEM_TEMRA", "Rose_CD8_modul2_Tcm_em",
                                 "Rose_CD8_modul3_Tem_emra", "Rose_CD8_modul1_Temra",
                                 "Terekhova_CD8_Tem_GZMKpos", "Terekhova_CD8_Tem_GZMBpos", "Terekhova_CD8_Temra"), order=T)
# FeaturePlot(subset(seur.human, subset=group.ident=="CD8_PBMC"),
#             features=c("Poon_CD8Naive", "Rose_CD8_modul4_Tnaive", "Terekhova_CD8_Tnaive",
#                        "Poon_CD8TEM_TEMRA",
#                        "Rose_CD8_modul2_Tcm_em", "Rose_CD8_modul3_Tem_emra", "Rose_CD8_modul1_Temra",
#                        "Terekhova_CD8_Tem_GZMKpos", "Terekhova_CD8_Tem_GZMBpos", "Terekhova_CD8_Temra"),
#             order=T, reduction="UMAP_50")
VlnPlot(subset(seur.human, subset=group.ident=="CD8_PBMC"),
            features=c("Poon_CD8Naive", "Rose_CD8_modul4_Tnaive", "Terekhova_CD8_Tnaive", "Poon_CD8TEM_TEMRA",
                       "Rose_CD8_modul2_Tcm_em", "Rose_CD8_modul3_Tem_emra", "Rose_CD8_modul1_Temra",
                       "Terekhova_CD8_Tem_GZMKpos", "Terekhova_CD8_Tem_GZMBpos", "Terekhova_CD8_Temra"),
            group.by = "new_clusters")+
  scale_fill_manual(values=cols_integrated)

VlnPlot(subset(seur.human, subset=group.ident=="CD8_PBMC"),
        features=c("Terekhova_CD8_Temra"),
        group.by = "new_clusters")+
  scale_fill_manual(values=cols_integrated)+
  labs(x="")+
  theme(legend.position="none")


## end 2.7 ####




# ****************
# 3. ANALYSIS ####
# ****************

#___________________________
## 3.1. Identify max GEP and 2nd max GEP per cell ####
# remove gep12 (because it's batch driven)
df <- gep_usage[, colnames(gep_usage)[!colnames(gep_usage) %in% "gep12_usage"] ]

# Get the max
df$score_max <- apply(df, 1, max)
df$gep_assign <- gsub("_.*", "", colnames(df)[apply(df, 1, which.max)])
head(df)
table(df$gep_assign, useNA="ifany")

# Get the 2nd max
df$score_2ndmax <- apply(df[,1:11], 1, function(x) sort(x, decreasing = T)[2])
df$gep_assign_2ndmax <- gsub("_.*", "", apply(df[,1:11], 1, function(x) names(x)[order(x, decreasing=T)[2]]))
head(df)

# add cluster and group.ident info
table(rownames(df)==rownames(seur.human@meta.data), useNA="ifany")
df <- cbind(df, seur.human@meta.data[,c("group.ident", "new_clusters")])


# Look at GEP usage in CD8 PBMC in clusters 13-17
df %>%
  as_tibble() %>%
  filter(group.ident=="CD8_PBMC" & new_clusters %in% 13:17) %>%
  mutate(scoremax_minus_2ndmax=score_max-score_2ndmax,
         new_clusters = paste0("CD8 in cluster ", new_clusters)) %>%
  # group_by(new_clusters, gep_assign, gep_assign_2ndmax) %>%
  # count() %>%
  ggplot(aes(x=gep_assign, y=scoremax_minus_2ndmax))+
  facet_wrap(~new_clusters)+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.1, size=0.1)+
  ylim(c(0,1))+
  labs(x="GEP with max usage", y="max GEP usage - 2nd max GEP usage")

df %>%
  as_tibble() %>%
  filter(group.ident=="CD8_PBMC" & new_clusters %in% 13:17) %>%
  # get nb of cells per gep assignment
  group_by(new_clusters, gep_assign, gep_assign_2ndmax) %>%
  summarise(ncells=n()) %>%
  # get %cells in each lineage
  ungroup() %>%
  group_by(new_clusters) %>%
  mutate(totalcells=sum(ncells),
         freq = ncells*100/totalcells) %>%
  ungroup() %>%
  # rename a few variables
  mutate(gep_assign=toupper(gep_assign)) %>%
  mutate(gep_assign=replace(gep_assign, !gep_assign%in%c("GEP3", "GEP4", "GEP5", "GEP6"), "other"),
         gep_assign = factor(gep_assign, levels=c("GEP3", "GEP4", "GEP5", "GEP6", "other"))) %>%
  # PLOT
  ggplot(aes(axis1=new_clusters, axis2=gep_assign, axis3=gep_assign_2ndmax, y=ncells)) +
  # geom_alluvium(aes(fill=new_clusters))+
  # scale_fill_manual(values=cols_integrated)+
  geom_alluvium(aes(fill=new_clusters))+
  scale_fill_manual(values=cols_integrated)+
  geom_stratum()+
  geom_text(stat="stratum", aes(label=after_stat(stratum)), size=8)+
  # theme_void()+
  scale_x_discrete(limits=c("cluster", "GEP max", "GEP 2nd max")) + theme_classic() +
  # scale_y_reverse()+
  # coord_flip()+
  theme(legend.position="none", axis.text.x = element_text(size=20))


## /end ####


#___________________________
## 3.2. Marker genes ####
DotPlot(subset(seur.human, subset=group.ident=="CD8_PBMC"),
        features=c("CCR7", "GZMK", "GZMB", "IL7R", "CD69", "KLRB1", "CCR5", "CCR2", "TNFSF13B"),
        group.by="new_clusters")+
  labs(x="", y="clusters", title="CD8 PBMCs")
## /end ####


#___________________________
## 3.3. Third analysis ####

## /end ####


#___________________________
## 3.4. Fourth analysis ####

## /end ####




# *******************************************
# 4. INVESTIGATE TEREKHOVA GD GENE LISTS ####
# *******************************************

# Import GD object and Terekhova GD gene list scoring
seur.gdt  <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.GD_03_16_23.RDS")
terekhova_gdscoring <- read.csv("./data/human-PBMC/HumanData_33_Tconv_Tem/terekhova_GDgenelists_scoring_integratedobject.csv", row.names=1)
head(terekhova_gdscoring)
colnames(terekhova_gdscoring)[2:6] <- gsub("\\.", "+", colnames(terekhova_gdscoring)[2:6])
table(rownames(seur.gdt@meta.data)==rownames(terekhova_gdscoring[terekhova_gdscoring$group.ident=="GD_PBMC",]), useNA="ifany")
seur.gdt@meta.data <- cbind(seur.gdt@meta.data,
                            terekhova_gdscoring[terekhova_gdscoring$group.ident=="GD_PBMC",colnames(terekhova_gdscoring) != "group.ident"])

# add the gep assign
table(rownames(seur.gdt@meta.data)==rownames(df[df$group.ident=="GD_PBMC",]), useNA="ifany")
seur.gdt@meta.data <- cbind(seur.gdt@meta.data,
                            df[df$group.ident=="GD_PBMC", !colnames(df) %in% c("group.ident", "new_clusters")])


# Plot gene signatures on GD object
SCpubr::do_FeaturePlot(seur.gdt,
                       features=colnames(seur.gdt@meta.data)[42:46],
                       reduction="umap",
                       viridis.palette = "B",
                       use_viridis = T,
                       order=F,
                       ncol=5
)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_22_CompareGeneLists/plots/genelists_with_terekhova/terekhova_GDscoring_GDumap_GDcellshighlight_orderF.jpeg", width=22, height=6)

# Plot gene signatures by gep assign
# SCpubr::do_FeaturePlot(seur.gdt,
#                        features=c("gep5_usage", "gep6_usage"),
#                        reduction="umap",
#                        viridis.palette = "B",
#                        use_viridis = T,
#                        order=F,
#                        ncol=2
# )

# Identify VD2 and VD1 cells
# simplify VG and VD
seur.gdt@meta.data <- seur.gdt@meta.data %>%
  mutate(tcr_gamma=case_when(
    TCR_Alpha_Gamma_V_gene_Dominant %in% grep("TRAV", TCR_Alpha_Gamma_V_gene_Dominant, value=T) ~ "Valpha",
    TCR_Alpha_Gamma_V_gene_Dominant %in% grep("TRGV9", TCR_Alpha_Gamma_V_gene_Dominant, value=T) ~ "VG9",
    TCR_Alpha_Gamma_V_gene_Dominant %in% grep("TRGV", TCR_Alpha_Gamma_V_gene_Dominant, value=T) ~ "VGother",
    .default = "unknown"
  ),
  tcr_delta=case_when(
    TCR_Beta_Delta_V_gene_Dominant %in% grep("TRBV", TCR_Beta_Delta_V_gene_Dominant, value=T) ~ "Vbeta",
    TCR_Beta_Delta_V_gene_Dominant %in% grep("TRDV1", TCR_Beta_Delta_V_gene_Dominant, value=T) ~ "VD1",
    TCR_Beta_Delta_V_gene_Dominant %in% grep("TRDV2", TCR_Beta_Delta_V_gene_Dominant, value=T) ~ "VD2",
    TCR_Beta_Delta_V_gene_Dominant %in% grep("TRDV3", TCR_Beta_Delta_V_gene_Dominant, value=T) ~ "VD3",
    .default = "unknown"
  ),
  ) %>%
    mutate(tcr_vd2vg9=case_when(
      cell.ident=="GD" & tcr_delta=="VD2" & tcr_gamma=="VG9" ~ "GD(TRDV2_TRGV9)",
      cell.ident=="GD" & tcr_delta %in% c("VD1", "VD3")      ~ "GD(TRDV2neg)",
      .default = "other"
    ))

# Violin plot
VlnPlot(seur.gdt, features=colnames(seur.gdt@meta.data)[42:46], group.by="new_clusters_GD", ncol=5, cols=cols_pbmc_gdt, same.y.lims =T)
ggsave("./scripts-in-progress/human-PBMC/HumanData_22_CompareGeneLists/plots/genelists_with_terekhova/terekhova_GDscoring_GDumap_vlnplot_bycluster.jpeg", width=20, height=5)

VlnPlot(seur.gdt, features=colnames(seur.gdt@meta.data)[42:46], group.by="gep_assign", ncol=5, same.y.lims =T)
ggsave("./scripts-in-progress/human-PBMC/HumanData_22_CompareGeneLists/plots/genelists_with_terekhova/terekhova_GDscoring_GDumap_vlnplot_byGEPmax.jpeg", width=20, height=5)

VlnPlot(subset(seur.gdt, subset=tcr_delta %in% c("VD1", "VD2", "VD3")), features=colnames(seur.gdt@meta.data)[42:46], group.by="tcr_delta", ncol=5, same.y.lims =T)
ggsave("./scripts-in-progress/human-PBMC/HumanData_22_CompareGeneLists/plots/genelists_with_terekhova/terekhova_GDscoring_GDumap_vlnplot_byVdelta_chain.jpeg", width=20, height=5)


# Ribbon plot GEP usage
seur.gdt@meta.data %>%
  as_tibble() %>%
  # get nb of cells per gep assignment
  group_by(new_clusters_GD, gep_assign) %>%
  summarise(ncells=n()) %>%
  # get %cells in each gep assignment
  ungroup() %>%
  group_by(new_clusters_GD) %>%
  mutate(totalcells=sum(ncells),
         freq = ncells*100/totalcells) %>%
  ungroup() %>%
  # rename a few variables
  mutate(gep_assign=toupper(gep_assign)) %>%
  mutate(gep_assign=replace(gep_assign, !gep_assign%in%c("GEP3", "GEP4", "GEP5", "GEP6"), "other"),
         gep_assign = factor(gep_assign, levels=c("GEP3", "GEP4", "GEP5", "GEP6", "other"))) %>%
  # filter(!tcr_delta %in% c("unknown", "Vbeta") & gep_assign != "other") %>%
  # PLOT
  ggplot(aes(axis1=new_clusters_GD, axis2=gep_assign, y=ncells)) +
  geom_alluvium(aes(fill=new_clusters_GD))+
  geom_stratum()+
  geom_text(stat="stratum", aes(label=after_stat(stratum)), size=8)+
  scale_fill_manual(values=cols_pbmc_gdt)+
  scale_x_discrete(limits=c("GD clusters", "GEP with max usage")) + theme_classic()+
  # theme_void()+
  # scale_y_reverse()+
  # coord_flip()+
  theme(legend.position="none")
ggsave("./scripts-in-progress/human-PBMC/HumanData_22_CompareGeneLists/plots/genelists_with_terekhova/GDpbmc_flowchart_clusters_to_gep.jpeg", width=7, height=5)


seur.gdt@meta.data %>%
  as_tibble() %>%
  # get nb of cells per gep assignment
  group_by(tcr_delta, gep_assign) %>%
  summarise(ncells=n()) %>%
  # get %cells in each gep assignment
  ungroup() %>%
  group_by(tcr_delta) %>%
  mutate(totalcells=sum(ncells),
         freq = ncells*100/totalcells) %>%
  ungroup() %>%
  # rename a few variables
  mutate(gep_assign=toupper(gep_assign)) %>%
  mutate(gep_assign=replace(gep_assign, !gep_assign%in%c("GEP3", "GEP4", "GEP5", "GEP6"), "other"),
         gep_assign = factor(gep_assign, levels=c("GEP3", "GEP4", "GEP5", "GEP6", "other"))) %>%
  filter(!tcr_delta %in% c("unknown", "Vbeta") & gep_assign != "other") %>%
  # PLOT
  ggplot(aes(axis1=tcr_delta, axis2=gep_assign, y=ncells)) +
  geom_alluvium(aes(fill=tcr_delta))+
  geom_stratum()+
  geom_text(stat="stratum", aes(label=after_stat(stratum)), size=8)+
  # scale_fill_manual(values=cols_pbmc_gdt)+
  scale_x_discrete(limits=c("Known Vdelta", "GEP with max usage")) + theme_classic()+
  # theme_void()+
  # scale_y_reverse()+
  # coord_flip()+
  theme(legend.position="none")
ggsave("./scripts-in-progress/human-PBMC/HumanData_22_CompareGeneLists/plots/genelists_with_terekhova/GDpbmc_flowchart_Vdelta_to_gep.jpeg", width=7, height=5)



