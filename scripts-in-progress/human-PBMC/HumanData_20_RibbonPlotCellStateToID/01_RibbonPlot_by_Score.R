# Purpose: make ribbon plot of cell state vs identity
# Author: Salomé Carcy
# Date: May 2023


# **************
# 1. IMPORT ####
# **************

library(Seurat)
library(GGally)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(ggalluvial)
source("./scripts-final/colors_universal.R")


# Import data
seur <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
DimPlot(seur, group.by = "new_clusters", label=T, repel=T, reduction="UMAP_50")+
  scale_color_manual(values=cols_integrated)


# Import GEPs
gep_topgenes <- read.csv("./data/human-thymus/HumanData_17_GEPsOnParkData/genes_per_GEP_df_2023-04-07.csv", row.names=1)
dim(gep_topgenes) # 12 GEPs

# table(gep_topgenes$GEP_1[1:200] %in% gep_topgenes$GEP_5[1:200])




# ********************
# 2. COMPUTE GEPs ####
# ********************

gep_list <- list()
for (i in 1:ncol(gep_topgenes)){
  print(paste0("GEP", i))
  # gep_allgenes <- gep_topgenes[1:200,i][!is.na(gep_topgenes[1:200,i])]
  gep_allgenes <- gep_topgenes[,i][!is.na(gep_topgenes[,i])]
  print(length(gep_allgenes))
  gep_list[[colnames(gep_topgenes)[i]]] <- gep_allgenes
}

seur.geps     <- AddModuleScore(seur, name = "GEP", features=gep_list)
# seur.geps <- readRDS("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/seuratobj_gepscores_allgenes.rds")

# Sanity check
SCpubr::do_FeaturePlot(seur.geps, reduction="UMAP_50", features=paste0("GEP", 1:length(gep_list)), ncol=6,
                       viridis_color_map = "B", order=T)
# ggsave("./scripts-in-progress/human-thymus/HumanData_20_RibbonPlotCellStateToID/plots/gapin_cNMF_geps_allgenes.jpeg", width=40, height=20)
SCpubr::do_FeaturePlot(subset(seur.geps, Tissue=="PBMC"), reduction="UMAP_50", features=gep_pbmc, ncol=3,
                       viridis_color_map = "B", order=T)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/gapin_cNMF_geps_allgenes_pbmc.jpeg", width=17, height=15)



# *************************************
# 3. THRESHOLD CELLS BASED ON GEPs ####
# *************************************

gep_pbmc <- c("GEP1", "GEP4", "GEP5", "GEP6", "GEP8", "GEP12")
df <- as.data.frame(seur.geps@meta.data[seur.geps@meta.data$Tissue=="PBMC",gep_pbmc])
# Check the distribution
ggplot(pivot_longer(df, cols=all_of(gep_pbmc), names_to="gep", values_to="score"))+
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
dfraw$gep_assign <- colnames(dfraw)[apply(dfraw, 1, which.max)]
head(dfraw)
tabl(dfraw$score_max<0) # 7 of them have negative score, which is not good
dfraw[dfraw$score_max<0,"gep_assign"] <- "undefined"
tabl(dfraw$gep_assign)

# Look at distribution of max score
ggplot(dfraw)+
  geom_density(aes(x=score_max, color=gep_assign))

# How many have a very little difference between the max score and 2nd max score?
dfraw$score_2ndmax <- apply(dfraw[,1:7], 1, function(x) sort(x, decreasing = T)[2])
tabl(dfraw$score_2ndmax<0)
tabl(dfraw$score_max-dfraw$score_2ndmax<0.01) # 7% of the cells have less than 0.01 difference in score between top score and 2nd top score
ggplot(dfraw)+
  geom_histogram(aes(x=score_max-score_2ndmax), bins=1000)+
  labs(y="# cells", title="Raw GEP score")+
  xlim(0,1) + ylim(0,400)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/method1_rawscore_diff_btw_max_2ndmax.jpeg", width=6, height=4)

# # Look at distribution of counts
# hist(df$GEP1, breaks=1000) # 0.1
# hist(df$GEP4, breaks=1000) # 0.15
# hist(df$GEP5, breaks=1000) # 0.1?...
# hist(df$GEP6, breaks=1000) # 0.1?...
# 
# Remove any score below 0.05
tabl(dfraw$score_max < 0.05)
dfraw[dfraw$score_max < 0.05, "score_max"] <- 0.0
table(dfraw$score_max == 0)
dfraw[dfraw$score_max == 0, "gep_assign"]   <- "undefined"

# First look
seur.gepsraw <- seur.geps
table(rownames(seur.gepsraw@meta.data[seur.gepsraw$Tissue=="PBMC",]) == rownames(dfraw))
seur.gepsraw@meta.data[seur.gepsraw$Tissue=="PBMC","gep_assign"] <- dfraw$gep_assign

# Final look
DimPlot(seur.gepsraw, group.by = "gep_assign", repel=T, reduction="UMAP_50") + scale_color_manual(values=c(brewer.pal(7, "Dark2"), "#f0f0f0"))
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/method1_rawscore_umap.jpeg", width=6, height=6)
VlnPlot(seur.gepsraw[,seur.gepsraw$Tissue=="PBMC"], group.by = "gep_assign", features=gep_pbmc, same.y.lims=T)
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
         score_2ndmax=replace_na(score_2ndmax,0.0),
         gep_2ndmax=ifelse(score_2ndmax != 0, gep[which(normscore==score_2ndmax)], "none")) |>
  # create column with GEP assignment
  mutate(gep_assign=ifelse(score_max-score_2ndmax>0.01, gep_max, "undefined")) |>
  ungroup()

# Quick visualization
head(dfnorm2)
tabl(dfnorm2[,c("gep_max", "gep_assign")])
dfnorm2 %>% select(cellid, gep_assign) %>% distinct() %>% count(gep_assign)

# Check distribution of normalized scores
ggplot(dfnorm2)+
  # geom_density(aes(x=normscore, color=factor(gep_assign, levels=c(gep_pbmc, "undefined"))))+
  # facet_wrap(~factor(gep, levels=gep_pbmc), ncol=4, scales="free_y")+
  # scale_color_manual(values=RColorBrewer::brewer.pal(8, "Paired"), name="Assigned to...")+
  # labs(y="normalized positive GEP score", title="(1) remove negative GEP scores (2) normalize x/xmax")
  # _______
  geom_violin(aes(x=factor(gep_assign, levels=c(gep_pbmc, "undefined")), y=normscore), width=1)+
  geom_jitter(aes(x=factor(gep_assign, levels=c(gep_pbmc, "undefined")), y=normscore, color=factor(gep_assign, levels=c(gep_pbmc, "undefined"))), size=0.1, width = 0.05)+
  facet_wrap(~factor(gep, levels=gep_pbmc), ncol=4)+
  scale_color_manual(values=RColorBrewer::brewer.pal(8, "Paired"), name="Assigned to...")+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  labs(y="normalized positive GEP score", x="Cells assigned to...")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/method4_%max_scoreassign.jpeg", width=10, height=8)


# Check distribution of unassigned cells
ggplot(dfnorm2 %>% filter(gep_assign=="undefined") %>% group_by(cellid) %>% filter(gep %in% c(gep_max, gep_2ndmax)) %>% ungroup(),
       aes(x=factor(gep, levels=gep_pbmc), y=normscore))+
  geom_violin(width=1)+
  geom_jitter(size=0.1, width = 0.05)+
  geom_line(aes(group=cellid), linewidth=0.1)+
  # facet_wrap(~factor(gep, levels=gep_pbmc), ncol=4)+
  # scale_color_manual(values=RColorBrewer::brewer.pal(8, "Paired"), name="Assigned to...")+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  labs(y="normalized positive GEP score", x="Cells assigned to...")

# Check distribution of difference between max and 2nd max score
ggplot(dfnorm2 %>% select(cellid, score_max, score_2ndmax) %>% distinct())+
  geom_histogram(aes(x=score_max-score_2ndmax), bins=1000)+
  geom_vline(xintercept=0.01)+
  labs(y="# cells", title="(1) remove negative GEP scores (2) normalize x/xmax")+
  xlim(0,1) + ylim(0,400)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/method4_%max_diff_btw_max_2ndmax.jpeg", width=6, height=4)


# Quick look on UMAP
seur.geps_norm2 <- seur.geps
dfnorm2.temp <- dfnorm2 |>
  select(cellid, gep_assign) |>
  distinct()
tabl(rownames(seur.geps_norm2@meta.data) %in% dfnorm2.temp$cellid)
tabl(rownames(seur.geps_norm2@meta.data[dfnorm2.temp$cellid,]) == dfnorm2.temp$cellid)
seur.geps_norm2@meta.data[dfnorm2.temp$cellid,"gep_assign"] <- dfnorm2.temp$gep_assign
DimPlot(seur.geps_norm2, group.by = "gep_assign", repel=T, reduction="UMAP_50") +
  scale_color_manual(values=c(brewer.pal(8, "Dark2"), "#f0f0f0"))
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/method4_%max_umap.jpeg", width=6, height=6)

# Plot it one at a time
Idents(seur.geps_norm2) <- "gep_assign"
DimPlot(seur.geps_norm2, cells.highlight =WhichCells(seur.geps_norm2, idents="GEP8"), repel=T, reduction="UMAP_50")

# Potentially a good method ####


# 3.4bis. Method 4b ####
ggplot(dfnorm2 %>% select(gep, normscore) %>% distinct())+
  geom_histogram(aes(x=normscore), bins = 100)+
  facet_wrap(~gep, ncol=4)+
  labs(y="GEP score", title="Normalized GEP score (%max)")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/method4bis_normscore_distribution.jpeg", width=6, height=4)
# ggplot(dfnorm2 %>% select(gep, normscore) %>% distinct() %>% filter(gep=="GEP5"))+
#   geom_histogram(aes(x=normscore), bins = 500)+
#   facet_wrap(~gep, ncol=4)+
#   labs(y="GEP score", title="Normalized GEP score (%max)")

dfnorm2b <- dfnorm2 |>
  # get back to wide
  select(cellid, gep, normscore) |>
  pivot_wider(names_from=gep, values_from = normscore, values_fill=0) |>
  # thresholds
  mutate(GEP1_threshold=ifelse(GEP1>0.2, T, F),
         GEP4_threshold=ifelse(GEP4>0.2, T, F),
         GEP5_threshold=ifelse(GEP5>0.35, T, F),
         GEP6_threshold=ifelse(GEP6>0.35, T, F),
         GEP8_threshold=ifelse(GEP8>0.4, T, F),
         GEP11_threshold=ifelse(GEP11>0.4, T, F),
         GEP12_threshold=ifelse(GEP12>0.1, T, F)) |>
  # pivot longer
  pivot_longer(cols=ends_with("threshold"), names_to="gep", values_to="pass") |>
  mutate(pass=as.numeric(pass),
         gep=gsub("_threshold", "", gep)) |>
  # keep only lines with GEPs that passed threshold (1 gep assigned)
  group_by(cellid) |>
  filter(pass==1 & sum(pass)==1) |>
  # keep only columns of interest
  select(cellid, gep) |>
  dplyr::rename(gep_assign=gep) |>
  distinct()
  
# Quick look on UMAP
seur.geps_norm2b <- seur.geps
tabl(rownames(seur.geps_norm2b@meta.data) %in% dfnorm2b$cellid)
tabl(rownames(seur.geps_norm2b@meta.data[dfnorm2b$cellid,]) == dfnorm2b$cellid)
seur.geps_norm2b@meta.data[dfnorm2b$cellid,"gep_assign"] <- dfnorm2b$gep_assign
seur.geps_norm2b@meta.data[seur.geps_norm2b$Tissue=="PBMC" & is.na(seur.geps_norm2b$gep_assign),"gep_assign"] <- "undefined"
tabl(seur.geps_norm2b@meta.data[,c("Tissue", "gep_assign")]) # sanity check
DimPlot(seur.geps_norm2b, group.by = "gep_assign", repel=T, reduction="UMAP_50") +
  scale_color_manual(values=c(brewer.pal(8, "Dark2"), "#f0f0f0"))
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/method4bis_umap.jpeg", width=6, height=6)




# *******************
# 4. RIBBON PLOT ####
# *******************

# Format data
counts <- seur.geps_norm2b@meta.data %>%
  as_tibble() %>%
  filter(Tissue=="PBMC") %>%
  # get nb of cells per gep assignment
  group_by(cell.ident, gep_assign) %>%
  summarise(ncells=n()) %>%
  # get %cells in each gep assignment
  ungroup() %>%
  group_by(cell.ident) %>%
  mutate(totalcells=sum(ncells),
         freq = ncells*100/totalcells) %>%
  ungroup() %>%
  # rename to "other" for GEP8,11,12
  mutate(gep_assign=replace(gep_assign, gep_assign%in%c("GEP8", "GEP11", "GEP12", "undefined"), "other"),
         gep_assign = factor(gep_assign, levels=c(gep_pbmc[1:4], "other"))) %>%
  # if no rename to "other"
  # mutate(gep_assign=factor(gep_assign, levels=c(gep_pbmc, "undefined"))) %>%
  # rename
  mutate(cell.ident=replace(cell.ident, cell.ident=="NKT", "iNKT"))
  

# Plot
ggplot(data=counts, aes(axis1=cell.ident, axis2=gep_assign, y=freq)) +
  geom_alluvium(aes(fill=cell.ident))+
  geom_stratum()+
  geom_text(stat="stratum", aes(label=after_stat(stratum)))+
  scale_fill_manual(values=c("#2ca25f", "#dd1c77", "#045a8d", "#8856a7", "#bdc9e1"), name="cell type")+
  # scale_x_discrete(limits=c("Cell type", "GEP")) + theme_classic()
  theme_void()+
  labs(title="method4bis (threshold on distribution %max)")
ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/method4bis_ribbon.jpeg", width=6, height=6)


# Sanity check
# VlnPlot(subset(seur.geps, subset= Tissue=="PBMC"), group.by = "cell.ident", features=c("GEP1", "GEP4", "GEP5", "GEP6"),
#         same.y.lims=T, ncol=4)+
#   scale_fill_manual(values=c("#2ca25f", "#dd1c77", "#045a8d", "#8856a7", "#bdc9e1"))

# dfplot.test <- seur.geps@meta.data[,c("Tissue", "cell.ident", gep_pbmc)] %>%
#   filter(Tissue=="PBMC") %>%
#   rownames_to_column("cellid") %>%
#   as_tibble() %>%
#   select(-Tissue) %>%
#   pivot_longer(cols=starts_with("GEP"), names_to="gep", values_to="score")
# 
# ggplot(dfplot.test, aes(x=factor(gep, levels=gep_pbmc), y=score))+
#   geom_violin()+
#   geom_jitter(width=0.01, size=0.1)+
#   facet_wrap(~cell.ident)+
#   theme(axis.text.x = element_text(angle=45, hjust=1))

