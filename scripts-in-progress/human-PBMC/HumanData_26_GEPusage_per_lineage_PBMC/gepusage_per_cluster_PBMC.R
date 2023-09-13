# Purpose: make ribbon plot of cell state vs identity
# Author: Salomé Carcy
# Date: August 2023


# **************
# 1. IMPORT ####
# **************

library(Seurat)
# library(GGally)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
library(RColorBrewer)
# library(ggalluvial)
source("./scripts-final/colors_universal.R")


# Import data
# seur.cd4  <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.CD4_03_16_23.RDS")
# seur.cd8  <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.CD8_03_16_23.RDS")
seur.gdt  <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.GD_03_16_23.RDS")
seur.mait <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.MAIT_03_16_23.RDS")
seur.nkt  <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.nkt_03_16_23.RDS")

# DimPlot(seur.cd4,  group.by = "new_clusters_CD4", label=T, repel=T,  reduction="umap")
# DimPlot(seur.cd8,  group.by = "new_clusters_CD8", label=T, repel=T,  reduction="umap")
DimPlot(seur.gdt,  group.by = "new_clusters_GD", label=T, repel=T,   reduction="umap")
DimPlot(seur.mait, group.by = "new_clusters_MAIT", label=T, repel=T, reduction="umap")
DimPlot(seur.nkt,  group.by = "new_clusters_NKT", label=T, repel=T,  reduction="umap")






# ******************************
# 2. GET GEP USAGE PER CELL ####
# ******************************

# Import GEP usage
gep_usage <- read.table("./data/human-PBMC/HumanData_20_RibbonPlotCellStateToID/cNMF_output/non_imputed_cNMF_allcells.usages.k_12.dt_0_02.consensus.txt", header=T)
dim(gep_usage)
colnames(gep_usage) <- paste0("GEP", c(2,5,3,1,4,12,6,7,8,10,9,11), "_usage") # rename column names
head(gep_usage)

# Check that every row sums to 1
hist(rowSums(gep_usage), breaks=100)
min(rowSums(gep_usage))
max(rowSums(gep_usage))

# Remove GEP12 that's batch-driven
gep_usage <- gep_usage[,paste0("GEP", 1:11, "_usage")]

# Get the gep with max usage per cell
gep_usage$score_max <- apply(gep_usage, 1, max)
gep_usage$gep_assign <- gsub("_.*", "", colnames(gep_usage)[apply(gep_usage, 1, which.max)])
head(gep_usage)
# table(gep_usage$score_max<0.5, useNA="ifany") # only ~24,000 cells have a max score below 0.5
# table(gep_usage$gep_assign, useNA="ifany")

# How many have a very little difference between the max score and 2nd max score?
gep_usage$score_2ndmax <- apply(gep_usage[,1:11], 1, function(x) sort(x, decreasing = T)[2])
table(gep_usage$score_max-gep_usage$score_2ndmax<0.1, useNA="ifany") # 12% of the cells have less than 0.1 difference in score between top score and 2nd top score
ggplot(gep_usage)+
  geom_histogram(aes(x=score_max-score_2ndmax), bins=1000)+
  labs(y="# cells", title="assign cells to GEPs based on cNMF usage")+
  xlim(0,1) + ylim(0,300)




# **********************************
# 3. PLOT GEP USAGE PER CLUSTER ####
# **********************************


## 3.1. NKT cells ####
# Add to NKT seurat object
gep_usage_nkt <- gep_usage[colnames(seur.nkt),]
table(rownames(gep_usage_nkt)==rownames(seur.nkt@meta.data), useNA="ifany")
seur.nkt@meta.data <- cbind(seur.nkt@meta.data, gep_usage_nkt)
# Plot on seurat object GEP usage
SCpubr::do_FeaturePlot(seur.nkt,
                       reduction="umap",
                       features=paste0("GEP", 1:11, "_usage"),
                       ncol=3, viridis_color_map = "B", order=F)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/plots/nkt_cNMFnonimput_geps_usage_orderFALSE.jpeg", width=17, height=20)

# Plot % NKT cells in each cluster, per GEP
# seur.nkt@meta.data %>%
#   rownames_to_column("cellid") %>%
#   as_tibble %>%
#   select(cellid, new_clusters_NKT, gep_assign) %>%
#   group_by(new_clusters_NKT, gep_assign) %>%
#   count() %>%
#   ggplot(aes(x=gep_assign, y=n, fill=new_clusters_NKT))+
#   geom_bar(position="stack", stat="identity")+
#   scale_fill_manual(values=cols_pbmc_nkt)+
#   labs(y="# cells")

# Plot GEP usage per cluster
p_nkt <- seur.nkt@meta.data %>%
  rownames_to_column("cellid") %>%
  as_tibble %>%
  select(cellid, new_clusters_NKT, gep3_usage, gep4_usage, gep5_usage, gep6_usage) %>%
  pivot_longer(cols=starts_with("gep"), names_to="gep", values_to="usage") %>%
  mutate(gep=gsub("_usage", "", gep),
         gep=toupper(gep)) %>%
  ggplot(aes(x=new_clusters_NKT, y=usage))+
  facet_wrap(~gep, nrow=1)+
  geom_violin(aes(fill=new_clusters_NKT), width=1)+
  geom_boxplot(outlier.shape=NA, width=0.2)+
  scale_fill_manual(values=cols_pbmc_nkt)+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1))+
  labs(x="", y="GEP usage", title="iNKT")+
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        strip.background=element_blank(),
        strip.text = element_text(size=20))

## end NKT ####


## 3.2. MAIT cells ####

# Add to MAIT seurat object
gep_usage_mait <- gep_usage[colnames(seur.mait),]
table(rownames(gep_usage_mait)==rownames(seur.mait@meta.data), useNA="ifany")
seur.mait@meta.data <- cbind(seur.mait@meta.data, gep_usage_mait)
# Plot on seurat object GEP usage
SCpubr::do_FeaturePlot(seur.mait,
                       reduction="umap",
                       features=paste0("gep", 1:11, "_usage"),
                       ncol=3, viridis_color_map = "B", order=F)

# Plot GEP usage per cluster
p_mait <- seur.mait@meta.data %>%
  rownames_to_column("cellid") %>%
  as_tibble %>%
  select(cellid, new_clusters_MAIT, gep3_usage, gep4_usage, gep5_usage, gep6_usage) %>%
  pivot_longer(cols=starts_with("gep"), names_to="gep", values_to="usage") %>%
  mutate(gep=gsub("_usage", "", gep),
         gep=toupper(gep)) %>%
  ggplot(aes(x=new_clusters_MAIT, y=usage))+
  facet_wrap(~gep, nrow=1)+
  geom_violin(aes(fill=new_clusters_MAIT), width=1)+
  geom_boxplot(outlier.shape=NA, width=0.2)+
  scale_fill_manual(values=cols_pbmc_mait)+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1))+
  labs(x="", y="GEP usage", title="MAIT")+
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        strip.background=element_blank(),
        strip.text = element_blank())

## end MAIT ####


## 3.3. GDT cells ####

# Add to GDT seurat object
gep_usage_gdt <- gep_usage[colnames(seur.gdt),]
table(rownames(gep_usage_gdt)==rownames(seur.gdt@meta.data), useNA="ifany")
seur.gdt@meta.data <- cbind(seur.gdt@meta.data, gep_usage_gdt)
# Plot on seurat object GEP usage
SCpubr::do_FeaturePlot(seur.gdt,
                       reduction="umap",
                       features=paste0("GEP", 1:11, "_usage"),
                       ncol=3, viridis_color_map = "B", order=F)

# Plot GEP usage per cluster
p_gdt <- seur.gdt@meta.data %>%
  rownames_to_column("cellid") %>%
  as_tibble %>%
  select(cellid, new_clusters_GD, gep3_usage, gep4_usage, gep5_usage, gep6_usage) %>%
  pivot_longer(cols=starts_with("gep"), names_to="gep", values_to="usage") %>%
  mutate(gep=gsub("_usage", "", gep),
         gep=toupper(gep)) %>%
  ggplot(aes(x=new_clusters_GD, y=usage))+
  facet_wrap(~gep, nrow=1)+
  geom_violin(aes(fill=new_clusters_GD), width=1)+
  geom_boxplot(outlier.shape=NA, width=0.2)+
  scale_fill_manual(values=cols_pbmc_gdt)+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1))+
  labs(x="", y="GEP usage", title="GD")+
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        strip.background=element_blank(),
        strip.text = element_blank())

## end GDT ####


## 3.4. CD4 cells ####

# Add to CD4 seurat object
gep_usage_cd4 <- gep_usage[colnames(seur.cd4),]
table(rownames(gep_usage_cd4)==rownames(seur.cd4@meta.data), useNA="ifany")
seur.cd4@meta.data <- cbind(seur.cd4@meta.data, gep_usage_cd4)
# Plot on seurat object GEP usage
SCpubr::do_FeaturePlot(seur.cd4,
                       reduction="umap",
                       features=paste0("gep", 1:11, "_usage"),
                       ncol=3, viridis_color_map = "B", order=F)

# Plot GEP usage per cluster
p_cd4 <- seur.cd4@meta.data %>%
  rownames_to_column("cellid") %>%
  as_tibble %>%
  select(cellid, new_clusters_CD4, gep3_usage, gep4_usage, gep5_usage, gep6_usage) %>%
  pivot_longer(cols=starts_with("gep"), names_to="gep", values_to="usage") %>%
  mutate(gep=gsub("_usage", "", gep),
         gep=toupper(gep)) %>%
  ggplot(aes(x=new_clusters_CD4, y=usage))+
  facet_wrap(~gep, nrow=1)+
  geom_violin(aes(fill=new_clusters_CD4), width=1)+
  geom_boxplot(outlier.shape=NA, width=0.2)+
  scale_fill_manual(values=cols_pbmc_cd4)+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1))+
  labs(x="", y="GEP usage", title="CD4")+
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        strip.background=element_blank(),
        strip.text = element_blank())

## end CD4 ####


## 3.8. CD8 cells ####

# Add to CD8 seurat object
seur.cd8 <- seur.cd8[,colnames(seur.cd8) %in% rownames(gep_usage)] # remove the MAITs contaminating (need to get seurat obj from Laurent)
gep_usage_cd8 <- gep_usage[colnames(seur.cd8),]
table(rownames(gep_usage_cd8)==rownames(seur.cd8@meta.data), useNA="ifany")
seur.cd8@meta.data <- cbind(seur.cd8@meta.data, gep_usage_cd8)
# Plot on seurat object GEP usage
SCpubr::do_FeaturePlot(seur.cd8,
                       reduction="umap",
                       features=paste0("gep", 1:11, "_usage"),
                       ncol=3, viridis_color_map = "B", order=F)

# Plot GEP usage per cluster
p_cd8 <- seur.cd8@meta.data %>%
  rownames_to_column("cellid") %>%
  as_tibble %>%
  select(cellid, new_clusters_CD8, gep3_usage, gep4_usage, gep5_usage, gep6_usage) %>%
  pivot_longer(cols=starts_with("gep"), names_to="gep", values_to="usage") %>%
  mutate(gep=gsub("_usage", "", gep),
         gep=toupper(gep)) %>%
  ggplot(aes(x=new_clusters_CD8, y=usage))+
  facet_wrap(~gep, nrow=1)+
  geom_violin(aes(fill=new_clusters_CD8), width=1)+
  geom_boxplot(outlier.shape=NA, width=0.2)+
  scale_fill_manual(values=cols_pbmc_cd8)+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1))+
  labs(x="", y="GEP usage")+
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        strip.background=element_blank(),
        strip.text = element_blank())

## end CD8 ####

library(patchwork)
p_nkt / p_mait / p_gdt / p_cd4
# ggsave("./scripts-in-progress/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/plots/cNMFnonimput_geps_usage_per_clust.jpeg", width=6, height=15)




# **********************************************
# 4. PLOT GEP ASSIGN PER CD4/CD8 EXPRESSION ####
# **********************************************

## 4.1. NKT cells ####

nkt_counts <- as.data.frame(t(as.data.frame(seur.nkt@assays$RNA@data[c("CD4", "CD8A"),])))
table(rownames(nkt_counts)==rownames(seur.nkt@meta.data), useNA="ifany")
nkt_counts$gep_assign <- seur.nkt@meta.data$gep_assign
table(nkt_counts$gep_assign, useNA="ifany")

nkt_counts %>%
  rownames_to_column("cellid") %>%
  as_tibble() %>%
  mutate(status=ifelse(CD4==0 & CD8A==0, "DN",
                       ifelse(CD4>0 & CD8A==0, "CD4+",
                              ifelse(CD4==0 & CD8A>0, "CD8+",
                                     ifelse(CD4>0 & CD8A>0, "DP", "other"))))) %>%
  group_by(gep_assign, status) %>%
  count() %>%
  # plot
  ggplot(aes(x=gep_assign, y=n, fill=gep_assign))+
  geom_bar(stat="identity")+
  facet_wrap(~factor(status, levels=c("DN", "CD4+", "CD8+", "DP")), nrow=1)+
  labs(x="", y="# cells", title="iNKT")+
  scale_fill_manual(values=c("gold", "#a40000", "#72bcd5", "#0a2e57"), name="GEP with max usage")+
  theme_cowplot()+
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1))
# ggsave("./scripts-in-progress/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/plots/nkt_cd4cd8_per_gep.jpeg", width=8, height=4)


## end NKT ####


## 4.2. MAIT cells ####

mait_counts <- as.data.frame(t(as.data.frame(seur.mait@assays$RNA@data[c("CD4", "CD8A"),])))
table(rownames(mait_counts)==rownames(seur.mait@meta.data), useNA="ifany")
mait_counts$gep_assign <- seur.mait@meta.data$gep_assign
table(mait_counts$gep_assign, useNA="ifany")

mait_counts %>%
  rownames_to_column("cellid") %>%
  as_tibble() %>%
  filter(!gep_assign %in% c("GEP7", "GEP11")) %>%
  mutate(status=ifelse(CD4==0 & CD8A==0, "DN",
                       ifelse(CD4>0 & CD8A==0, "CD4+",
                              ifelse(CD4==0 & CD8A>0, "CD8+",
                                     ifelse(CD4>0 & CD8A>0, "DP", "other"))))) %>%
  group_by(gep_assign, status) %>%
  count() %>%
  ggplot(aes(x=gep_assign, y=n, fill=gep_assign))+
  geom_bar(stat="identity")+
  facet_wrap(~factor(status, levels=c("DN", "CD4+", "CD8+", "DP")), nrow=1)+
  labs(x="", y="# cells", title="MAIT")+
  scale_fill_manual(values=c("gold", "#a40000", "#72bcd5", "#0a2e57"), name="GEP with max usage")+
  theme_cowplot()+
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1))
# ggsave("./scripts-in-progress/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/plots/mait_cd4cd8_per_gep.jpeg", width=8, height=4)


## end MAIT ####


## 4.3. GDT cells ####

gdt_counts <- as.data.frame(t(as.data.frame(seur.gdt@assays$RNA@data[c("CD4", "CD8A"),])))
table(rownames(gdt_counts)==rownames(seur.gdt@meta.data), useNA="ifany")
gdt_counts$gep_assign <- seur.gdt@meta.data$gep_assign
table(gdt_counts$gep_assign, useNA="ifany")

gdt_counts %>%
  rownames_to_column("cellid") %>%
  as_tibble() %>%
  filter(!gep_assign %in% c("GEP2", "GEP7", "GEP11")) %>%
  mutate(status=ifelse(CD4==0 & CD8A==0, "DN",
                       ifelse(CD4>0 & CD8A==0, "CD4+",
                              ifelse(CD4==0 & CD8A>0, "CD8+",
                                     ifelse(CD4>0 & CD8A>0, "DP", "other"))))) %>%
  group_by(gep_assign, status) %>%
  count() %>%
  ggplot(aes(x=gep_assign, y=n, fill=gep_assign))+
  geom_bar(stat="identity")+
  facet_wrap(~factor(status, levels=c("DN", "CD4+", "CD8+", "DP")), nrow=1)+
  labs(x="", y="# cells", title="GD")+
  scale_fill_manual(values=c("gold", "#a40000", "#72bcd5", "#0a2e57"), name="GEP with max usage")+
  theme_cowplot()+
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1))
# ggsave("./scripts-in-progress/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/plots/gdt_cd4cd8_per_gep.jpeg", width=8, height=4)

## end GDT ####
