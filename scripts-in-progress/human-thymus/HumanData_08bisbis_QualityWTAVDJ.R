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
path.plots <- "~/Projects/HumanThymusProject/data/human-thymus/HumanData_08_QualityControl"

seur <- readRDS("~/Projects/HumanThymusProject/data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
cols_integrated <- c("0" = "#f4c40f", "1" = "#b75347", "2" = "#d8443c", "3" = "#e09351", "4" = "#2b9b81", 
                     "5" = "#421401", "6" = "#92c051", "7" = "#9f5691", "8" = "#17154f", "9" = "#74c8c3", 
                     "10" = "#5a97c1", "11" = "gold", "12" = "#a40000", "13" = "#72bcd5", "14" = "grey50",
                     "15" = "orange", "16" = "blueviolet", "17" = "#0a2e57", "18" = "#bdbdbd")
DimPlot(seur, group.by = "new_clusters", label=T, repel=T, reduction="UMAP_50")+
  scale_color_manual(values=cols_integrated)


# Export metadata
# rownames(seur)[stringr::str_detect(rownames(seur), "^RP[SL]")]
seur[["percent.rb"]] <- PercentageFeatureSet(seur, pattern = "^RP[SL]")
df <- as.data.frame(seur@meta.data)

# df %>%
#   as_tibble() %>%
#   select(Batch, Donor, Method, cell.ident) %>%
#   distinct() %>%
#   arrange(Batch) %>%
#   print(n=50)
# table(df[,c("cell.ident", "Batch")])


# **********************
# 2. QUALITY PLOTS ####
# **********************

# Sanity check
# test <- CreateSeuratObject(counts=seur@assays$RNA@counts, meta.data = seur@meta.data[,c("Batch", "Method")])
# FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Method")

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




# **************************************
# 3. TRY DIFFERENT METADATA COLUMNS ####
# **************************************

# Look at number of unique values in each metadata column (if there are more than 50, it will be hard to facet)
columns_of_interest <- c()
for(i in 1:ncol(df)){
  print(colnames(df)[i])
  print(length(unique(df[,i])))
  if(length(unique(df[,i]))<50){columns_of_interest <- c(columns_of_interest, colnames(df)[i])}
}

# Save plots for metadata columns with less than 50 unique values
for(column in columns_of_interest){
  print(column)
  
  p <- ggplot(df %>% filter(Method=="WTA_VDJ"), aes(x=nCount_RNA, y=nFeature_RNA, color=Batch))+
    facet_wrap(facet=column)+
    geom_point(size=0.1)+
    # scale_color_manual(values=c("#404040", "#ca0020"), name="TCR info?")+
    labs(x="# reads", y="# genes detected", title=paste0("Facet by ", column))
  
  nb_panels <- length(ggplot_build(p)$layout$layout$PANEL)
  filename <- file.path(path.plots, paste0("doublecurve_metadatacols/facet_by_", column, ".jpeg"))
  
  if(nb_panels <= 3){ggsave(filename, p, width=6, height=4)}
  else if(nb_panels == 4){ggsave(filename, p, width=6, height=6)}
  else if(nb_panels <= 9){ggsave(filename, p, width=9, height=9)}
  else if(nb_panels >= 10){ggsave(filename, p, width=12, height=12)}
}

# Check what are the remaining columns
columns_remaining <- colnames(df)[!colnames(df) %in% c(columns_of_interest, "nCount_RNA", "nFeature_RNA")]
df_remaining <- df[,columns_remaining]
is_discrete   <- function(vec) all(is.numeric(x)) && all(x%%1==0)
is_continuous <- function(vec) all(is.numeric(vec)) && !is_discrete(vec)

# Continuous columns
ggplot(df, aes(x=nCount_RNA, y=nFeature_RNA, color=TCR_Beta_Delta_V_gene_Dominant))+
  facet_wrap(~Batch)+
  geom_point(size=0.1)+
  # scale_color_gradient(low="#ffeda0", high="#f03b20")+
  labs(x="# reads", y="# genes detected")




# *********************************
# 4. PLOT SOME GENE EXPRESSION ####
# *********************************

# Get the gene expression of a few genes (housekeeping or not)
genes_of_interest <- c("ACTB", "PTPRC", "CD4", "CD8A")
counts <- t(as.data.frame(seur@assays$RNA@data[genes_of_interest,]))

# Merge with df
# table(rownames(counts)==rownames(df)) # all true
counts <- cbind(counts, df %>% select(nCount_RNA, nFeature_RNA, Batch))


# Plot some gene expressio
for(gene in genes_of_interest){
  print(gene)
  
  p <- ggplot(counts, aes_string(x="nCount_RNA", y="nFeature_RNA", color=gene))+
    facet_wrap(~Batch)+
    geom_point(size=0.1)+
    scale_color_gradient(low="#ffeda0", high="#f03b20")+
    labs(x="# reads", y="# genes detected", title=gene)
  
  filename <- file.path(path.plots, paste0("doublecurve_genexp/gene_", gene, ".jpeg"))
  ggsave(filename, p, width=9, height=9)
}




# ***************************
# 5. SEPARATE TWO CURVES ####
# ***************************

test <- df %>%
  filter(Method=="WTA_VDJ") %>%
  select(nCount_RNA, nFeature_RNA, Batch, Method)

test <- test %>%
  mutate(ratio = nCount_RNA/nFeature_RNA)

ggplot(test)+
  facet_wrap(~Batch)+
  geom_point(aes(x=nCount_RNA, y=nFeature_RNA), size=0.01)+
  geom_abline(intercept=100, slope=2500/7500, color="red")+
  xlim(c(0,15000))+
  ylim(c(0,3000))+
  # scale_color_manual(values=c("#d01c8b", "#4dac26"))+
  labs(x="# reads", y="# genes detected")
ggsave(file.path(path.plots, "attempt_separation.jpeg"), width=6, height=6)

ggplot(test)+
  facet_wrap(~Batch)+
  geom_point(aes(x=nCount_RNA, y=ratio), size=0.01)+
  geom_abline(intercept=1.8, slope=1.7/10000, color="red")+
  labs(x="# reads", y="#reads / #genes")

