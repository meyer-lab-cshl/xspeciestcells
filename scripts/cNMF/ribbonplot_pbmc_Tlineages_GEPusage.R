# Purpose: make ribbon plot of cell state vs identity
# Author: Salomé Carcy
# Date: August 2023


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
setwd("~/Projects/HumanThymusProject/")
source("./scripts-final/colors_universal.R")


# Import data
# seur <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
seur <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_08_28_23.RDS")
DimPlot(seur, group.by = "new_clusters", label=T, repel=T, reduction="UMAP_50")+
  scale_color_manual(values=cols_integrated)




# ********************
# 2. COMPUTE GEPs ####
# ********************

# # Import GEPs
# gep_topgenes <- read.csv("./data/human-thymus/HumanData_17_GEPsOnParkData/genes_per_GEP_df_2023-04-07.csv", row.names=1)
# dim(gep_topgenes) # 12 GEPs
# test <- read.csv("~/Downloads/genes_per_GEP_df.csv", row.names=1)

# table(gep_topgenes$GEP_1[1:200] %in% gep_topgenes$GEP_5[1:200])

# gep_list <- list()
# for (i in 1:ncol(gep_topgenes)){
#   print(paste0("GEP", i))
#   # gep_allgenes <- gep_topgenes[1:200,i][!is.na(gep_topgenes[1:200,i])]
#   gep_allgenes <- gep_topgenes[,i][!is.na(gep_topgenes[,i])]
#   print(length(gep_allgenes))
#   gep_list[[colnames(gep_topgenes)[i]]] <- gep_allgenes
# }
# 
# seur.geps     <- AddModuleScore(seur, name = "GEP", features=gep_list, seed=123)
# seur.geps <- readRDS("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/seuratobj_gepscores_allgenes.rds")
# 
# # Sanity check
# gep_pbmc <- c("GEP1", "GEP4", "GEP5", "GEP6", "GEP8", "GEP12")
# SCpubr::do_FeaturePlot(subset(seur.geps, Tissue=="PBMC"), reduction="UMAP_50", features=gep_pbmc, ncol=3,
#                        viridis_color_map = "B", order=T)
# # ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/gapin_cNMF_geps_allgenes_pbmc_orderTRUE.jpeg", width=17, height=15)




# ******************************
# 3. GET GEP USAGE PER CELL ####
# ******************************

# Import GEP usage
# gep_usage <- read.table("./data/human-PBMC/HumanData_20_RibbonPlotCellStateToID/cNMF_output/imputed_cNMF.usages.k_12.dt_0_02.consensus.txt", header=T)
gep_usage <- read.table("./data/human-PBMC/HumanData_20_RibbonPlotCellStateToID/cNMF_output/non_imputed_cNMF_allcells.usages.k_12.dt_0_02.consensus.txt", header=T)
dim(gep_usage)
nrow(gep_usage)==ncol(seur) # rows are cells
head(gep_usage)

# Check that every row sums to 1
hist(rowSums(gep_usage), breaks=100)
min(rowSums(gep_usage))
max(rowSums(gep_usage))

# Add to seurat object
# colnames(gep_usage) <- paste0("gep", 1:12, "_usage")
colnames(gep_usage) <- paste0("gep", c(2,5,3,1,4,12,6,7,8,10,9,11), "_usage")
table(rownames(gep_usage)==rownames(seur@meta.data))
seur@meta.data <- cbind(seur@meta.data, gep_usage)

# Plot
gepusage_pbmc <- c("gep3_usage", "gep4_usage", "gep5_usage", "gep6_usage", "gep11_usage", "gep12_usage")
SCpubr::do_FeaturePlot(subset(seur, Tissue=="PBMC"),
                       reduction="UMAP_50",
                       # features=gepusage_pbmc,
                       features=sort(colnames(gep_usage)),
                       # features=paste0("GEP_Scores_", 1:12),
                       ncol=3, viridis_color_map = "B", order=F)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/gapin_cNMFnonimput_geps_usage_pbmc_orderFALSE_reordered.jpeg", width=17, height=20)




# ****************************
# 4. ASSIGN CELLS TO GEPs ####
# ****************************

# Keep only PBMCs and remove gep12 (because it's batch driven)
df <- gep_usage[grep("PBMC", rownames(gep_usage), value=T), colnames(gep_usage)[!colnames(gep_usage) %in% "gep12_usage"] ]
dim(df) # 41,238 cells; 40986 with MAITs removed

# Get the max
df$score_max <- apply(df, 1, max)
df$gep_assign <- gsub("_.*", "", colnames(df)[apply(df, 1, which.max)])
head(df)
table(df$score_max<0.5, useNA="ifany") # only ~11,000 cells have a max score below 0.5
table(df$gep_assign, useNA="ifany")

# Look at distribution of max score
# ggplot(df)+
#   geom_density(aes(x=score_max, color=gep_assign))

# How many have a very little difference between the max score and 2nd max score?
df$score_2ndmax <- apply(df[,1:11], 1, function(x) sort(x, decreasing = T)[2])
table(df$score_max-df$score_2ndmax<0.1) # 11% of the cells have less than 0.1 difference in score between top score and 2nd top score
ggplot(df)+
  geom_histogram(aes(x=score_max-score_2ndmax), bins=1000)+
  labs(y="# cells", title="assign cells to GEPs based on GEP with highest/max usage")+
  xlim(0,1) #+ ylim(0,300)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/method1_rawscore_diff_btw_max_2ndmax.jpeg", width=6, height=4)

# First look
table(rownames(seur@meta.data[seur$Tissue=="PBMC",]) == rownames(df))
seur@meta.data[seur$Tissue=="PBMC","gep_assign"] <- df$gep_assign

# Final look
DimPlot(seur, group.by = "gep_assign", repel=T, reduction="UMAP_50") + scale_color_manual(values=c(brewer.pal(7, "Dark2"), "#f0f0f0"))
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/method5_gepusage_umap.jpeg", width=6, height=6)
# VlnPlot(seur.geps[,seur.geps$Tissue=="PBMC"], group.by = "gep_assign", features=gep_pbmc, same.y.lims=T)
# ggsave("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/method_comparison/vlnpreview_gepscomparisononrawscores.jpeg", width=10, height=10)


# Quick look at iNKTs in GEP3
df_inkt <- df[grep("NKT", rownames(df), value=T),]
dim(df_inkt) # 1801 cells
df_inkt$gep_assign_2ndmax <- gsub("_.*", "", apply(df_inkt[,1:11], 1, function(x) names(x)[order(x, decreasing=T)[2]]))
head(df_inkt)
ggplot(df_inkt)+
  geom_histogram(aes(x=score_max-score_2ndmax, fill=gep_assign), bins=1000)+
  labs(y="# cells", title="assign cells to GEPs based on cNMF usage")+
  xlim(0,1)
# difference between max usage and 2nd max usage is lower when cells are assigned to GEP3
df_inkt %>%
  filter(gep_assign=="gep3") %>%
  group_by(gep_assign_2ndmax) %>%
  summarize(usage_diff_mean=mean(score_max-score_2ndmax),
            usage_diff_std=sd(score_max-score_2ndmax),
            ncells=n())
# 664 cells / 862 assigned to GEP3 have GEP4 as second highest usage (77%)
# the average difference in usage between GEP3 and GEP4 is 0.2 (± 0.1)
df_inkt %>%
  # group_by(gep_assign, gep_assign_2ndmax) %>%
  # summarize(usage_diff_mean=mean(score_max-score_2ndmax),
  #           usage_diff_std=sd(score_max-score_2ndmax),
  #           ncells=n())
  mutate(diff_usage = score_max-score_2ndmax,
         gep_assign_2ndmax=toupper(gep_assign_2ndmax),
         gep_assign=toupper(gep_assign),
         gep_assign=paste0("iNKTs in ", gep_assign)) %>%
  ggplot(aes(x=factor(gep_assign_2ndmax, levels=paste0("GEP", 1:11)), y=diff_usage))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.1, size=0.1)+
  ylim(0,1)+
  facet_wrap(~gep_assign)+
  theme_bw()+
  labs(x="GEP with 2nd highest usage", y="(max usage) - (2nd highest usage)", title="Usage difference between top 2 GEPs")+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        strip.text=element_text(size=15))
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/method5_inkt_usagediff_top2geps.jpeg", width=6, height=6)




# *******************
# 5. RIBBON PLOT ####
# *******************

# Format data
counts <- seur@meta.data %>%
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
  # rename a few variables
  mutate(gep_assign=toupper(gep_assign),
         cell.ident=replace(cell.ident, cell.ident=="NKT", "iNKT"))


# Main figure
counts %>%
  mutate(gep_assign=replace(gep_assign, !gep_assign%in%c("GEP3", "GEP4", "GEP5", "GEP6"), "other"),
         gep_assign = factor(gep_assign, levels=c("GEP3", "GEP4", "GEP5", "GEP6", "other"))) %>%
  filter(gep_assign != "other") %>%
ggplot(aes(axis1=cell.ident, axis2=gep_assign, y=freq)) +
  geom_alluvium(aes(fill=cell.ident))+
  geom_stratum()+
  geom_text(stat="stratum", aes(label=after_stat(stratum)), size=8, angle=90)+
  scale_fill_manual(values=c("#2ca25f", "#dd1c77", "#045a8d", "#8856a7", "#bdc9e1"), name="cell type")+
  # scale_x_discrete(limits=c("Cell type", "GEP")) + theme_classic()
  theme_void()+
  scale_y_reverse()+
  # coord_flip()+
  theme(legend.position="none")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/plots/method5_ribbon_nonimput3.jpeg", width=4, height=10)



# # Supp figure
# counts %>%
#   mutate(gep_assign = factor(gep_assign, levels=c("GEP1", "GEP3", "GEP4", "GEP5", "GEP6", "GEP8", "GEP12"))) %>%
# ggplot(aes(axis1=cell.ident, axis2=gep_assign, y=freq)) +
#   geom_alluvium(aes(fill=cell.ident))+
#   geom_stratum()+
#   geom_text(stat="stratum", aes(label=after_stat(stratum)))+
#   # ggrepel::geom_text_repel(
#   #   aes(label = gep_assign),
#   #   stat = "stratum", size = 4, direction = "y", nudge_x = .3
#   # ) +
#   scale_fill_manual(values=c("#2ca25f", "#dd1c77", "#045a8d", "#8856a7", "#bdc9e1"), name="cell type")+
#   # scale_x_discrete(limits=c("Cell type", "GEP")) + theme_classic()
#   theme_void()+
#   labs(title="")


# cluster to gep flowchart
seur@meta.data %>%
  as_tibble() %>%
  filter(Tissue=="PBMC") %>%
  # get nb of cells per gep assignment
  group_by(gep_assign, new_clusters) %>%
  summarise(ncells=n()) %>%
  # get %cells in each cluster
  ungroup() %>%
  group_by(gep_assign) %>%
  mutate(totalcells=sum(ncells),
         freq = ncells*100/totalcells) %>%
  ungroup() %>%
  # rename a few variables
  mutate(gep_assign=toupper(gep_assign)) %>%
  mutate(gep_assign=replace(gep_assign, !gep_assign%in%c("GEP3", "GEP4", "GEP5", "GEP6"), "other"),
         gep_assign = factor(gep_assign, levels=c("GEP3", "GEP4", "GEP5", "GEP6", "other"))) %>%
  filter(gep_assign != "other" & new_clusters %in% 9:17) %>%
  # PLOT
  ggplot(aes(axis1=gep_assign, axis2=new_clusters, y=freq)) +
  geom_alluvium(aes(fill=new_clusters))+
  geom_stratum()+
  geom_text(stat="stratum", aes(label=after_stat(stratum)), size=8)+
  scale_fill_manual(values=cols_integrated)+
  # scale_x_discrete(limits=c("Cell type", "GEP")) + theme_classic()
  theme_void()+
  # scale_y_reverse()+
  # coord_flip()+
  theme(legend.position="none")


# lineage to cluster to gep flowchart
seur@meta.data %>%
  as_tibble() %>%
  filter(Tissue=="PBMC") %>%
  # get nb of cells per gep assignment
  group_by(cell.ident, gep_assign, new_clusters) %>%
  summarise(ncells=n()) %>%
  # get %cells in each lineage
  ungroup() %>%
  group_by(cell.ident) %>%
  mutate(totalcells=sum(ncells),
         freq = ncells*100/totalcells) %>%
  ungroup() %>%
  # rename a few variables
  mutate(gep_assign=toupper(gep_assign)) %>%
  mutate(gep_assign=replace(gep_assign, !gep_assign%in%c("GEP3", "GEP4", "GEP5", "GEP6"), "other"),
         gep_assign = factor(gep_assign, levels=c("GEP3", "GEP4", "GEP5", "GEP6", "other"))) %>%
  filter(gep_assign != "other" & new_clusters %in% 9:17) %>%
  # PLOT
  ggplot(aes(axis1=cell.ident, axis3=gep_assign, axis2=new_clusters, y=freq)) +
  # geom_alluvium(aes(fill=new_clusters))+
  # scale_fill_manual(values=cols_integrated)+
  geom_alluvium(aes(fill=cell.ident))+
  scale_fill_manual(values=cols_lineages)+
  geom_stratum()+
  geom_text(stat="stratum", aes(label=after_stat(stratum)), size=8)+
  # scale_x_discrete(limits=c("Cell type", "GEP")) + theme_classic()
  theme_void()+
  # scale_y_reverse()+
  # coord_flip()+
  theme(legend.position="none")



# # GEP usage for each donor
# gepperdonor <- function(d, age){
#   seur.geps@meta.data %>%
#     as_tibble() %>%
#     filter(Tissue=="PBMC" & Donor==d) %>%
#     # get nb of cells per gep assignment
#     group_by(cell.ident, gep_assign) %>%
#     summarise(ncells=n()) %>%
#     # get %cells in each gep assignment
#     ungroup() %>%
#     group_by(cell.ident) %>%
#     mutate(totalcells=sum(ncells),
#            freq = ncells*100/totalcells) %>%
#     ungroup() %>%
#     # remove a lineage if less than 100 cells
#     filter(totalcells > 100) %>%
#     # rename a few variables
#     mutate(gep_assign=toupper(gep_assign),
#            cell.ident=replace(cell.ident, cell.ident=="NKT", "iNKT")) %>%
#     # keep only GEP1,4,5,6
#     mutate(gep_assign=replace(gep_assign, !gep_assign%in%c("GEP1", "GEP4", "GEP5", "GEP6"), "other"),
#            gep_assign = factor(gep_assign, levels=c("GEP1", "GEP4", "GEP5", "GEP6", "other"))) %>%
#     filter(gep_assign != "other") %>%
#     # Plot
#     ggplot(aes(axis1=cell.ident, axis2=gep_assign, y=freq)) +
#     geom_alluvium(aes(fill=cell.ident))+
#     geom_stratum()+
#     geom_text(stat="stratum", aes(label=after_stat(stratum)), size=8)+
#     scale_fill_manual(values=c("#2ca25f", "#dd1c77", "#045a8d", "#8856a7", "#bdc9e1"), name="cell type")+
#     labs(title=paste0("Donor ", d, " (", age, "y.o.)"))+
#     theme_void()+
#     theme(legend.position="none", title=element_text(size=20))
# }
# 
# gepperdonor(d=5, age="24")
# gepperdonor(d=6, age="41")
# gepperdonor(d=7, age="68")
# gepperdonor(d=11, age="38")
