# Purpose: Check quality of data (WTA vs WTA_VDJ)
# Author: Salomé Carcy
# Date: March 2023


# **************
# 1. IMPORT ####
# **************

# Import librairies
library(Seurat)
library(cowplot)
library(ggplot2)
library(tidyverse)


# Import data
path <- "~/Projects/HumanThymusProject/data/human-thymus/HumanData_17_GEPsOnParkData/"

seur <- readRDS("~/Projects/HumanThymusProject/data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
cols_integrated <- c("0" = "#f4c40f", "1" = "#b75347", "2" = "#d8443c", "3" = "#e09351", "4" = "#2b9b81", 
                     "5" = "#421401", "6" = "#92c051", "7" = "#9f5691", "8" = "#17154f", "9" = "#74c8c3", 
                     "10" = "#5a97c1", "11" = "gold", "12" = "#a40000", "13" = "#72bcd5", "14" = "grey50",
                     "15" = "orange", "16" = "blueviolet", "17" = "#0a2e57", "18" = "#bdbdbd")
DimPlot(seur, group.by = "new_clusters", label=T, repel=T, reduction="UMAP_50")+
  scale_color_manual(values=cols_integrated)


# Export metadata
df <- as.data.frame(seur@meta.data)



# **********************
# 2. QUALITY PLOTS ####
# **********************

# Look at count/gene relationship per method
ggplot(df, aes(x=nCount_RNA, y=nFeature_RNA, color=Method))+
  facet_wrap(~Batch)+
  geom_point(size=0.1)+
  scale_color_manual(values=c("#d01c8b", "#4dac26"))+
  labs(x="# reads", y="# genes detected")
# ggsave("~/Projects/HumanThymusProject/data/human-thymus/HumanData_08_QualityControl/qc_WTA_VDJ.jpeg", width=8, height=8)


# Check which batches have WTA_VDJ
table(df[,c("Batch", "Method")]) # D, F,G,H

# Create column specifying whether a cell has TCR information
df <- df %>%
  mutate(TCR_info=ifelse(is.na(TRAV10_TRAJ18)==T, "no", "yes"))

table(df$TCR_info) # 9540 yes

# Plot again
ggplot(df, aes(x=nCount_RNA, y=nFeature_RNA, color=TCR_info))+
  facet_wrap(~Batch)+
  geom_point(size=0.1)+
  scale_color_manual(values=c("#404040", "#ca0020"), name="TCR info?")+
  labs(x="# reads", y="# genes detected")
# ggsave("~/Projects/HumanThymusProject/data/human-thymus/HumanData_08_QualityControl/qc_WTA_VDJ_TCRinfo.jpeg", width=8, height=8)




# Try different metadata columns
ggplot(df %>% filter(Method=="WTA_VDJ"), aes(x=nCount_RNA, y=nFeature_RNA, color=Batch))+
  facet_wrap(~Donor)+
  geom_point(size=0.1)+
  # scale_color_manual(values=c("#404040", "#ca0020"), name="TCR info?")+
  labs(x="# reads", y="# genes detected")

table(df[, c("TCR_info", "Method")])


test <- GetAssayData(seur, "data")
test$gene <- rownames(test)

hist(test, breaks=100, xlab="normalized counts", main="Distribution of normalized counts")

