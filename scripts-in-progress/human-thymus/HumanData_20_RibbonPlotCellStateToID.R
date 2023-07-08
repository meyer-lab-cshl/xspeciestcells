# Purpose: make ribbon plot of cell state vs identity
# Author: Salomé Carcy
# Date: May 2023


# **************
# 1. IMPORT ####
# **************

library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(ggalluvial)


# Import data
seur <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
cols_integrated <- c("0" = "#f4c40f", "1" = "#b75347", "2" = "#d8443c", "3" = "#e09351", "4" = "#2b9b81", 
                     "5" = "#421401", "6" = "#92c051", "7" = "#9f5691", "8" = "#17154f", "9" = "#74c8c3", 
                     "10" = "#5a97c1", "11" = "gold", "12" = "#a40000", "13" = "#72bcd5", "14" = "grey50",
                     "15" = "orange", "16" = "blueviolet", "17" = "#0a2e57", "18" = "#bdbdbd")
DimPlot(seur, group.by = "new_clusters", label=T, repel=T, reduction="UMAP_50")+
  scale_color_manual(values=cols_integrated)


# Import GEPs
gep_topgenes <- read.csv("./data/human-thymus/HumanData_17_GEPsOnParkData/genes_per_GEP_df_2023-04-07.csv", row.names=1)
dim(gep_topgenes) # 12 GEPs and 100 genes per GEP

# table(gep_topgenes$GEP_1[1:200] %in% gep_topgenes$GEP_5[1:200])


# ********************
# 2. COMPUTE GEPs ####
# ********************

gep_list <- list()
for (i in 1:ncol(gep_topgenes)){
  print(paste0("GEP", i))
  gep_allgenes <- gep_topgenes[1:200,i][!is.na(gep_topgenes[1:200,i])]
  # gep_allgenes <- na.omit(gep_topgenes[,i])
  print(length(gep_allgenes))
  gep_list[[colnames(gep_topgenes)[i]]] <- gep_allgenes
}

seur.geps     <- AddModuleScore(seur, name = "GEP", features=gep_list)
# seur.geps <- readRDS("./data/human-thymus/HumanData_22_CompareGeneLists/seuratobj_gepscores_top200genes.rds")

# Sanity check
SCpubr::do_FeaturePlot(seur.geps, reduction="UMAP_50", features=paste0("GEP", 1:length(gep_list)), ncol=6,
                       viridis_color_map = "B", order=T)
# ggsave("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/gapin_cNMF_geps_top200genes.jpeg", width=40, height=20)




# *************************************
# 3. THRESHOLD CELLS BASED ON GEPs ####
# *************************************

gep_pbmc <- c("GEP1", "GEP4", "GEP5", "GEP6", "GEP8", "GEP11", "GEP12")
df <- as.data.frame(seur.geps@meta.data[seur.geps@meta.data$Tissue=="PBMC",gep_pbmc])
# Check the distribution
ggplot(pivot_longer(df, cols=gep_pbmc, names_to="gep", values_to="score"))+
  geom_histogram(aes(x=score), bins = 100)+
  facet_wrap(~gep, ncol=4)+
  # _______
  # geom_violin(aes(x=gep, y=score), width=1)+
  # geom_boxplot(aes(x=gep, y=score), outlier.shape = NA, width=0.05)+
  # geom_jitter(aes(x=gep, y=score), size=0.1, width = 0.05)+
  # _______
  labs(y="GEP score", title="Raw GEP score")



## 3.1. Method 1: Directly compare score levels ####
dfraw <- df

# Get the max
dfraw$score_max <- apply(dfraw[,1:7], 1, max)
dfraw$gep_max <- colnames(dfraw)[apply(dfraw, 1, which.max)]
head(dfraw)
tabl(dfraw$gep_max)
tabl(dfraw$score_max<0) # 7 of them have negative score, which is not good

# Look at distribution of max score
ggplot(dfraw)+
  geom_density(aes(x=score_max, color=gep_max))

# How many have a very little difference between the max score and 2nd max score?
dfraw$score_2ndmax <- apply(dfraw[,1:7], 1, function(x) sort(x, decreasing = T)[2])
tabl(dfraw$score_max-dfraw$score_2ndmax<0.01) # 7% of the cells have less than 0.01 difference in score between top score and 2nd top score
ggplot(dfraw)+
  geom_histogram(aes(x=score_max-score_2ndmax), bins=1000)+
  labs(y="# cells", title="Raw GEP score")+
  xlim(0,1) + ylim(0,400)

# First look
seur.gepsraw <- seur.geps
table(rownames(seur.gepsraw@meta.data[seur.gepsraw$Tissue=="PBMC",]) == rownames(dfraw))
seur.gepsraw@meta.data[seur.gepsraw$Tissue=="PBMC","gep_max"] <- dfraw$gep_max


# # Look at distribution of counts
# hist(df$GEP1, breaks=1000) # 0.1
# hist(df$GEP4, breaks=1000) # 0.15
# hist(df$GEP5, breaks=1000) # 0.1?...
# hist(df$GEP6, breaks=1000) # 0.1?...
# 
# # Remove any score below 0.1
# table(df$score_max < 0.05)
# df[df$gep_max == "GEP1" & df$score_max < 0.05, "score_max"] <- 0.0
# df[df$gep_max == "GEP4" & df$score_max < 0.05, "score_max"] <- 0.0
# df[df$gep_max == "GEP5" & df$score_max < 0.05, "score_max"] <- 0.0
# df[df$gep_max == "GEP6" & df$score_max < 0.05, "score_max"] <- 0.0
# table(df$score_max == 0)
# df[df$score_max == 0, "gep_max"]   <- "other"

# Final look
DimPlot(seur.gepsraw, group.by = "gep_max", repel=T, reduction="UMAP_50") + scale_color_manual(values=c(brewer.pal(7, "Dark2"), "#f0f0f0"))
# ggsave("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/method_comparison/umappreview_gepscomparisononrawscores.jpeg", width=7, height=6)
VlnPlot(seur.gepsraw[,seur.gepsraw$Tissue=="PBMC"], group.by = "gep_max", features=gep_pbmc, same.y.lims=T)
# ggsave("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/method_comparison/vlnpreview_gepscomparisononrawscores.jpeg", width=10, height=10)

# How many PBMCs are unassigned?
# tabl(seur.geps@meta.data[seur.geps@meta.data$Tissue=="PBMC", "gep_max"]) # 4,621 unassigned vs ...29,661 assigned

# Take a look at unassigned cells on UMAP
# DimPlot(subset(seur.geps, subset= Tissue=="PBMC" & gep_max=="other"), group.by = "new_clusters", repel=T, reduction="UMAP_50") +
#   scale_color_manual(values=cols_integrated)
# FeaturePlot(subset(seur.geps, subset= Tissue=="PBMC" & gep_max=="other"),
#             features="GEP1",
#             reduction="UMAP_50")

# Look in 2D scatter plot
library(GGally)
ggpairs(dfraw, columns=1:4, aes(color=gep_max, alpha=0.9), legend=c(1,1),
        upper=list(continuous=wrap("points", size=0.05)),
        lower=list(continuous=wrap("points", size=0.05)))+
  scale_color_manual(values=c(brewer.pal(7, "Dark2"), "#f0f0f0")) +
  scale_fill_manual(values=c(brewer.pal(7, "Dark2"), "#f0f0f0")) +
  theme_cowplot()+
  theme(legend.position="right")
# ggsave("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/method_comparison/ggpairs_pbmc_allgeps_colorbygepcategory_onrawscores.jpeg", width=10, height=10)



## 3.2. Method 2: Normalize distribution of scores between 0 and 1 ####
dfnorm <- as.data.frame(apply(df, 2, function(x) x+abs(min(x)))) # get everything to be >0
dfnorm <- as.data.frame(apply(dfnorm, 2, function(x) x/max(x))) # scale between [0,1]

# Get the max
dfnorm$score_max <- apply(dfnorm[,1:7], 1, max)
dfnorm$gep_max <- colnames(dfnorm)[apply(dfnorm, 1, which.max)]
head(dfnorm)
tabl(dfnorm$gep_max)

# Check distribution of normalized counts
ggplot(pivot_longer(dfnorm, cols=gep_pbmc, names_to="gep", values_to="score"))+
  geom_histogram(aes(x=score), bins = 100)+
  facet_wrap(~gep, ncol=4)+
  # _______
  # geom_violin(aes(x=gep, y=score), width=1)+
  # geom_boxplot(aes(x=gep, y=score), outlier.shape = NA, width=0.05)+
  # geom_jitter(aes(x=gep, y=score), size=0.1, width = 0.05)+
  # _______
  labs(y="normalized GEP score", title="Normalized GEP score: (x-xmin)/(xmax-xmin)")

# Sanity checks
ggplot(dfnorm)+
  geom_density(aes(x=score_max, color=gep_max))+
  xlim(0,1)

# How many have a very little difference between the max score and 2nd max score?
dfnorm$score_2ndmax <- apply(dfnorm[,1:7], 1, function(x) sort(x, decreasing = T)[2])
tabl(dfnorm$score_max-dfnorm$score_2ndmax<0.01) # 5% of the cells have less than 0.01 difference in score between top score and 2nd top score
ggplot(dfnorm)+
  geom_histogram(aes(x=score_max-score_2ndmax), bins=1000)+
  labs(y="# cells", title="Normalized GEP score: (x-xmin)/(xmax-xmin)")+
  xlim(0,1) + ylim(0,400)

  
# Add to seurat object
seur.gepsnorm <- seur.geps
# table(rownames(seur.gepsnorm@meta.data[seur.gepsnorm$Tissue=="PBMC",]) == rownames(dfnorm))
seur.gepsnorm@meta.data[seur.gepsnorm$Tissue=="PBMC","gep_max"] <- dfnorm$gep_max

# Take a look
DimPlot(seur.gepsnorm, group.by = "gep_max", repel=T, reduction="UMAP_50") + scale_color_manual(values=c(brewer.pal(7, "Dark2"), "#f0f0f0"))
# ggsave("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/method_comparison/umappreview_gepscomparisononnormscores.jpeg", width=7, height=6)
VlnPlot(seur.gepsnorm[,seur.gepsnorm$Tissue=="PBMC"], group.by = "gep_max", features=gep_pbmc, same.y.lims = T)
# ggsave("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/method_comparison/vlnpreview_gepscomparisononnormscores.jpeg", width=10, height=10)


# Make ggpairs plot
ggpairs(dfnorm, columns=1:4, aes(color=gep_max, alpha=0.9), legend=c(1,1),
        upper=list(continuous=wrap("points", size=0.05)),
        lower=list(continuous=wrap("points", size=0.05)))+
  scale_color_manual(values=c(brewer.pal(7, "Dark2"), "#f0f0f0")) +
  scale_fill_manual(values=c(brewer.pal(7, "Dark2"), "#f0f0f0")) +
  theme_cowplot()+
  theme(legend.position="right")
# ggsave("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/method_comparison/ggpairs_pbmc_allgeps_colorbygepcategory_onnormscores.jpeg", width=10, height=10)



## 3.3. Method 3: Directly compare score levels, do it by nb of SD away ####
head(dfraw)
# gep_mean <- apply(dfraw[,1:7], 2, mean)
# gep_sd <- apply(dfraw[,1:7], 2, sd)
for(i in gep_pbmc){
  columname <- paste0("z_", i)
  print(columname)
  avg <- mean(dfraw[,i])
  std <- sd(dfraw[,i])
  dfraw[,columname] <- (dfraw[,i]-avg)/std
}
ggplot(pivot_longer(dfraw, cols=paste0("z_", gep_pbmc), names_to="gep", values_to="zscore"))+
  # geom_histogram(aes(x=zscore), bins = 100)+
  # facet_wrap(~gep, ncol=4)+
  # _______
  # geom_violin(aes(x=gep, y=zscore), width=1)+
  geom_boxplot(aes(x=gep, y=zscore), outlier.shape = NA, width=0.05)+
  geom_jitter(aes(x=gep, y=zscore), size=0.1, width = 0.05)+
  # _______
  labs(y="GEP zscore")

# Max zscore per cell
dfraw$zscore_max <- apply(dfraw[,11:17], 1, max)
dfraw$zgep_max <- colnames(dfraw[,11:17])[apply(dfraw[,11:17], 1, which.max)]
head(dfraw)
tabl(dfraw$zgep_max)
tabl(dfraw$zscore_max<0) # 37 of them have negative zscore, which is not good

# Look at distribution of max score
ggplot(dfraw)+
  geom_density(aes(x=zscore_max, color=zgep_max))

# How many have a very little difference between the max score and 2nd max score?
dfraw$zscore_2ndmax <- apply(dfraw[,11:17], 1, function(x) sort(x, decreasing = T)[2])
tabl(dfraw$zscore_max-dfraw$zscore_2ndmax<0.01) # <1% of the cells have less than 0.01 difference in zscore between top score and 2nd top score
ggplot(dfraw)+
  geom_density(aes(x=zscore_max-zscore_2ndmax))+
  xlim(0,5)

# First look
# seur.geps_zscore <- seur.geps
# table(rownames(seur.geps_zscore@meta.data[seur.geps_zscore$Tissue=="PBMC",]) == rownames(dfraw))
# seur.geps_zscore@meta.data[seur.geps_zscore$Tissue=="PBMC","zgep_max"] <- dfraw$zgep_max
# DimPlot(seur.geps_zscore, group.by = "zgep_max", repel=T, reduction="UMAP_50") + scale_color_manual(values=c(brewer.pal(7, "Dark2"), "#f0f0f0"))

# Probably a bad idea, because now we have very high number of cells in GEP6, GEP8, GEP11, GEP12 ####

## 3.4. Method 4: Remove scores that are below 0, then normalize distribution of positive scores between 0 and 1 ####
dfnorm2 <- df
head(dfnorm2)
dfnorm2 <- dfnorm2 |>
  rownames_to_column("cellid") |>
  as_tibble() |>
  pivot_longer(cols=gep_pbmc, names_to="gep", values_to="score") |>
  filter(score>0) |>
  # normalize gep score
  group_by(gep) |>
  mutate(normscore=score/max(score)) |>
  ungroup() |>
  # get max score per cell
  group_by(cellid) |>
  mutate(score_max=max(normscore),
         gep_max=gep[which.max(normscore)],
         score_2ndmax=sort(normscore, decreasing=T)[2],
         score_2ndmax=replace_na(score_2ndmax,0.0)) |>
  ungroup()

head(dfnorm2)
tabl(dfnorm2$gep_max)

# Check distribution of normalized counts
ggplot(dfnorm2)+
  geom_histogram(aes(x=normscore), bins = 100)+
  facet_wrap(~gep, ncol=4)+
  # _______
  # geom_violin(aes(x=gep, y=score), width=1)+
  # geom_boxplot(aes(x=gep, y=score), outlier.shape = NA, width=0.05)+
  # geom_jitter(aes(x=gep, y=score), size=0.1, width = 0.05)+
  # _______
  labs(y="normalized positive GEP score", title="(1) remove negative GEP scores (2) normalize x/xmax")
# Check distribution of difference between max and 2nd max score
ggplot(dfnorm2 %>% select(cellid, score_max, score_2ndmax) %>% distinct())+
  geom_histogram(aes(x=score_max-score_2ndmax), bins=1000)+
  labs(y="# cells", title="(1) remove negative GEP scores (2) normalize x/xmax")+
  xlim(0,1) + ylim(0,400)






# *******************
# 4. RIBBON PLOT ####
# *******************

# Format data
counts <- seur.gepsraw@meta.data %>%
  as_tibble() %>%
  filter(Tissue=="PBMC") %>%
  # get nb of cells per gep assignment
  group_by(cell.ident, gep_max) %>%
  summarise(ncells=n()) %>%
  # get %cells in each gep assignment
  ungroup() %>%
  group_by(cell.ident) %>%
  mutate(totalcells=sum(ncells),
         freq = ncells*100/totalcells) %>%
  ungroup() %>%
  # rename
  mutate(# gep_max = ifelse(gep_max=="GEP1", "Th17?\n(GEP1)",
                          # ifelse(gep_max=="GEP4", "Temra\n(GEP4)",
                          #        ifelse(gep_max=="GEP5", "Tnaive\n(GEP5)",
                          #               ifelse(gep_max=="GEP6", "Tcm\n(GEP6)", "other")))),
         gep_max = factor(gep_max, levels=gep_pbmc),
         cell.ident=replace(cell.ident, cell.ident=="NKT", "iNKT"))
  

# Plot
ggplot(data=counts, aes(axis1=cell.ident, axis2=gep_max, y=freq)) +
  geom_alluvium(aes(fill=cell.ident))+
  geom_stratum()+
  geom_text(stat="stratum", aes(label=after_stat(stratum)))+
  scale_fill_manual(values=c("#2ca25f", "#dd1c77", "#045a8d", "#8856a7", "#bdc9e1"), name="cell type")+
  # scale_x_discrete(limits=c("Cell type", "GEP")) + theme_classic()
  theme_void()
# ggsave("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/method_comparison/ribbon_cellstateproportions_rawscores.jpeg", width=6, height=6)

# Sanity check
# VlnPlot(subset(seur.geps, subset= Tissue=="PBMC"), group.by = "cell.ident", features=c("GEP1", "GEP4", "GEP5", "GEP6"),
#         same.y.lims=T, ncol=4)+
#   scale_fill_manual(values=c("#2ca25f", "#dd1c77", "#045a8d", "#8856a7", "#bdc9e1"))


