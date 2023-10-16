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

# Import SCENIC regulon activity matrix
scenic_matrix <- read.csv("./data/human-PBMC/HumanData_28_SCENICpruning/binarized/prunedGenes_full_obj_new_scenic_100_runs_09202023_raw_counts_reg80-target80_regulons_pruned_final_consensus_100runs_auc_mtx.csv")
dim(scenic_matrix)
head(scenic_matrix)

scenic_matrix_bin <- read.csv("./data/human-PBMC/HumanData_28_SCENICpruning/binarized/prunedGenes_full_obj_new_scenic_100_runs_09202023_raw_counts_reg80-target80_regulons_pruned_final_consensus_100runs_binarized_matrix.csv")
head(scenic_matrix_bin)

# put two matrices in a list
table(colnames(scenic_matrix)==colnames(scenic_matrix_bin), useNA="ifany")
scenic_matrix_list <- list("mat_auc"=scenic_matrix, "mat_bin"=scenic_matrix_bin)


# Import GEP usage
gep_usage <- read.table("./data/human-PBMC/HumanData_20_RibbonPlotCellStateToID/cNMF_output/non_imputed_cNMF_allcells.usages.k_12.dt_0_02.consensus.txt", header=T)
colnames(gep_usage) <- paste0("gep", c(2,5,3,1,4,12,6,7,8,10,9,11), "_usage")
dim(gep_usage)
nrow(gep_usage)==nrow(scenic_matrix) # rows are cells
head(gep_usage)





# *****************
# 2. FUNCTIONS ####
# *****************





# ****************
# 3. ANALYSIS ####
# ****************

#___________________________
## 3.1. Get matrix with the regulons activity & GEP assignment per cell ####

# GEP usage: get max GEP usage per cell
# gep_usage$score_max <- apply(gep_usage, 1, max)
gep_usage$gep_assign <- gsub("_.*", "", colnames(gep_usage)[apply(gep_usage, 1, which.max)])
head(gep_usage)
table(gep_usage$gep_assign)

# bind gep_usage matrix with scenic matrices
lapply(scenic_matrix_list, function(x) table(rownames(gep_usage)==x$Cell, useNA="ifany"))
scenic_matrix_list_pbmc <- lapply(scenic_matrix_list, function(x){
  x$gep_assign <- gep_usage$gep_assign
  x <- x[x$Cell %in% grep("PBMC", x$Cell, value=T),]
  # remove cells belonging to gep2,7,11 (too few cells)
  x <- x[!x$gep_assign %in% c("gep2", "gep7", "gep11"),]
  })
# sanity checks
# lapply(scenic_matrix_list_pbmc, function(x) dim(x)) # 40,947 cells
# lapply(scenic_matrix_list_pbmc, function(x) table(x$gep_assign, useNA="ifany"))

## /end ####


#___________________________
## 3.2. Keep only regulons that are active in at least 50% of cells ####

# Compute % cells active for each regulon and per GEP
df_regulonsactive <- data.frame()
vector_regulons <- colnames(scenic_matrix_list_pbmc$mat_auc)[!colnames(scenic_matrix_list_pbmc$mat_auc) %in% c("Cell", "gep_assign")]
length(vector_regulons)==length(unique(vector_regulons))
for(gep in unique(scenic_matrix_list_pbmc$mat_auc$gep_assign)){
  print(gep)
  df_temp <- NULL
  matrix_subgep <- NULL
  matrix_subgep <- lapply(scenic_matrix_list_pbmc, function(x) x[x$gep_assign==gep,])
  df_temp <- data.frame("regulon"=vector_regulons,
                        "gep"=gep,
                        "propcells_bin_positive"=colSums(matrix_subgep$mat_bin[,vector_regulons]>0)*100/nrow(matrix_subgep$mat_bin),
                        "propcells_auc_positive"=colSums(matrix_subgep$mat_auc[,vector_regulons]>0)*100/nrow(matrix_subgep$mat_auc),
                        "mean_auc"  = apply(matrix_subgep$mat_auc[,vector_regulons], 2, mean),
                        "median_auc"= apply(matrix_subgep$mat_auc[,vector_regulons], 2, median))
  rownames(df_temp) <- NULL
  df_regulonsactive <- rbind(df_regulonsactive, df_temp)
}

head(df_regulonsactive)
table(df_regulonsactive$gep, useNA="ifany") # should all have 184, the number of regulons
table(df_regulonsactive$regulon, useNA="ifany") # should all have 5 (number of geps)


# Plot %cells per GEP that have AUC>0 or "active regulon" (bin=1)
plot_grid(
  ggplot(df_regulonsactive, aes(x=gep, y=propcells_auc_positive))+
    geom_boxplot(outlier.shape=NA) + geom_jitter(width=0.2)+
    theme_cowplot(),
  ggplot(df_regulonsactive, aes(x=gep, y=propcells_bin_positive))+
    geom_boxplot(outlier.shape=NA) + geom_jitter(width=0.2)+
    geom_text_repel(data=subset(df_regulonsactive, propcells_bin_positive>80), aes(label=regulon), force=100)+
    theme_cowplot(),
  ncol=2)
# check correlation between the mean_AUC and propcells that have active regulon
ggplot(df_regulonsactive, aes(x=propcells_bin_positive, y=mean_auc))+
  geom_point()+
  facet_wrap(~gep)+
  geom_text_repel(data=subset(df_regulonsactive, propcells_bin_positive>80), aes(label=regulon), force=100, max.overlaps=20)


# check regulons with highest mean_auc
library(ggrepel)
ggplot(df_regulonsactive, aes(x=gep, y=mean_auc))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.2)+
  geom_text_repel(data=subset(df_regulonsactive, mean_auc>0.2), aes(label=regulon))+
  theme_cowplot()+
  labs(x="PBMCs belonging to GEP...")

# jpeg("./scripts-in-progress/human-PBMC/HumanData_28_SCENICpruning/plots/heatmap_pbmc_meanAUC_per_GEP.jpeg", width=1000, height=2000, res = 100)
htmp1 <- pheatmap::pheatmap(df_regulonsactive %>% select(regulon, mean_auc, gep) %>% pivot_wider(names_from=gep, values_from=mean_auc) %>% column_to_rownames("regulon"),
                   scale="row",
                   cluster_cols=F,
                   fontsize_col=20,
                   main="Mean AUC of PBMCs from each GEP (z-scored)")
# dev.off()
htmp1_roworder <- htmp1$tree_row$labels[htmp1$tree_row$order]


# PRUNE REGULONS
# Identify regulons that are active in >50% of cells across all GEPs
df_regulonsactive_pruned <- df_regulonsactive %>%
  as_tibble() %>%
  group_by(regulon) %>%
  filter(sum(propcells_bin_positive>50) > 0) %>%
  ungroup()
length(unique(df_regulonsactive$regulon)) # 184 regulons
length(unique(df_regulonsactive_pruned$regulon)) # 27 regulons

# Identify regulons that have the highest SD of their mean, and highest SD of propcells_bin_positive
df_regulonsactive_pruned <- df_regulonsactive_pruned %>%
  group_by(regulon) %>%
  mutate(sd_meanauc_across_geps=sd(mean_auc),
         sd_propcellsactive_across_geps=sd(propcells_bin_positive)) %>%
  ungroup()

# check which regulons have the highest variation in mean_AUC and %cells with active regulon
plot_grid(
  ggplot(df_regulonsactive_pruned, aes(x=sd_meanauc_across_geps, y=mean_auc))+
    geom_point()+
    facet_wrap(~gep, nrow=1)+
    geom_text_repel(data=subset(df_regulonsactive_pruned, sd_meanauc_across_geps>0.03), aes(label=regulon), force=100)+
    theme_cowplot()+
    labs(x="standard deviation of mean AUC", title="sd(mean_auc)"),
  ggplot(df_regulonsactive_pruned, aes(x=sd_propcellsactive_across_geps, y=mean_auc))+
    geom_point()+
    facet_wrap(~gep, nrow=1)+
    geom_text_repel(data=subset(df_regulonsactive_pruned, sd_propcellsactive_across_geps>10), aes(label=regulon), force=100)+
    theme_cowplot()+
    labs(x="standard deviation of proportion of cells with regulon active", title="sd(% cells with active regulon)"),
  nrow=2
  )


# HEATMAP WITH MEAN AUC
# jpeg("./scripts-in-progress/human-PBMC/HumanData_28_SCENICpruning/plots/heatmap_pbmc_meanAUC_per_GEP_pruned.jpeg", width=1000, height=1500, res = 150)
htmp2 <- pheatmap::pheatmap(df_regulonsactive_pruned %>%
                     select(regulon, mean_auc, gep) %>%
                     pivot_wider(names_from=gep, values_from=mean_auc) %>%
                     # arrange(factor(regulon, levels=htmp1_roworder[htmp1_roworder %in% regulon])) %>%
                     column_to_rownames("regulon"),
                   scale="row",
                   cluster_cols=F,
                   cluster_rows=T,
                   fontsize_col=20,
                   main="Pruned regulons (active in >50% cells for at least 1 GEP)\ncolor scale: mean AUC across cells")
# dev.off()
htmp2_roworder <- htmp2$tree_row$labels[htmp2$tree_row$order]


# HEATMAP WITH PROP CELLS WITH REGULON ACTIVE
jpeg("./scripts-in-progress/human-PBMC/HumanData_28_SCENICpruning/plots/heatmap_pbmc_propcellsregulonactive_per_GEP_pruned.jpeg", width=1000, height=1500, res = 150)
pheatmap::pheatmap(df_regulonsactive_pruned %>%
                              select(regulon, propcells_bin_positive, gep) %>%
                              pivot_wider(names_from=gep, values_from=propcells_bin_positive) %>%
                              arrange(factor(regulon, levels=htmp2_roworder[htmp2_roworder %in% regulon])) %>%
                              column_to_rownames("regulon"),
                            scale="none",
                            cluster_cols=F,
                            cluster_rows=F,
                            fontsize_col=20,
                            main="Pruned regulons (active in >50% cells for at least 1 GEP)\ncolor scale: %cells with regulon active")
dev.off()



## /end ####


#___________________________
## 3.3.  Keep only regulons that have the most variability in their mean_AUC and/or %cells active ####

# Identify regulons that have the highest SD of their mean, and highest SD of propcells_bin_positive
df_regulonsactive <- df_regulonsactive %>%
  group_by(regulon) %>%
  mutate(sd_meanauc_across_geps=sd(mean_auc),
         sd_propcellsactive_across_geps=sd(propcells_bin_positive)) %>%
  ungroup()

# check which regulons have the highest variation in mean_AUC and %cells with active regulon
plot_grid(
  ggplot(df_regulonsactive, aes(x=sd_meanauc_across_geps, y=mean_auc))+
    geom_point()+
    facet_wrap(~gep, nrow=1)+
    geom_text_repel(data=subset(df_regulonsactive, regulon%in%c("CEBPD", "FOSL2")), aes(label=regulon), force=100)+
    theme_cowplot()+
    labs(x="standard deviation of mean AUC", title="sd(mean_auc)"),
  ggplot(df_regulonsactive, aes(x=sd_propcellsactive_across_geps, y=mean_auc))+
    geom_point()+
    facet_wrap(~gep, nrow=1)+
    geom_text_repel(data=subset(df_regulonsactive, regulon%in%c("CEBPD", "FOSL2")), aes(label=regulon), force=200, max.overlaps=20)+
    theme_cowplot()+
    labs(x="standard deviation of proportion of cells with regulon active", title="sd(% cells with active regulon)"),
  nrow=2
)

ggplot(df_regulonsactive, aes(x=sd_propcellsactive_across_geps, y=sd_meanauc_across_geps))+
  geom_point()+
  geom_text_repel(data=subset(df_regulonsactive %>% filter(gep=="gep3"), sd_propcellsactive_across_geps>10 | sd_meanauc_across_geps>0.02),
                  aes(label=regulon), force=50, max.overlaps=20)+
  labs(x="sd(% cells with active regulon)", y="sd(mean_auc)", title="standard deviations across GEPs")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_28_SCENICpruning/plots/ggplot_sd_propcellsregulonactive_meanAUC.jpeg", width=6, height=6)


# Prune regulons
df_regulonsactive_pruned2 <- df_regulonsactive %>%
  filter(sd_propcellsactive_across_geps>10 | sd_meanauc_across_geps>0.02)
length(unique(df_regulonsactive_pruned2$regulon)) # 27 regulons
length(intersect(unique(df_regulonsactive_pruned2$regulon), unique(df_regulonsactive_pruned$regulon))) # 10 in common with previous pruning method


# HEATMAP WITH MEAN AUC
# jpeg("./scripts-in-progress/human-PBMC/HumanData_28_SCENICpruning/plots/heatmap_pbmc_meanAUC_per_GEP_pruned.jpeg", width=1000, height=1500, res = 150)
htmp3 <- pheatmap::pheatmap(df_regulonsactive_pruned2 %>%
                              select(regulon, mean_auc, gep) %>%
                              pivot_wider(names_from=gep, values_from=mean_auc) %>%
                              # arrange(factor(regulon, levels=htmp1_roworder[htmp1_roworder %in% regulon])) %>%
                              column_to_rownames("regulon"),
                            scale="row",
                            cluster_cols=F,
                            cluster_rows=T,
                            fontsize_col=20,
                            main="Pruned regulons (based on SD of mean_AUC and %cells active)\ncolor scale: mean AUC across cells")
# dev.off()
htmp3_roworder <- htmp3$tree_row$labels[htmp3$tree_row$order]


# HEATMAP WITH PROP CELLS WITH REGULON ACTIVE
# jpeg("./scripts-in-progress/human-PBMC/HumanData_28_SCENICpruning/plots/heatmap2_pbmc_propcellsregulonactive_per_GEP_pruned.jpeg", width=1000, height=1500, res = 150)
pheatmap::pheatmap(df_regulonsactive_pruned2 %>%
                     select(regulon, propcells_bin_positive, gep) %>%
                     pivot_wider(names_from=gep, values_from=propcells_bin_positive) %>%
                     arrange(factor(regulon, levels=htmp3_roworder[htmp3_roworder %in% regulon])) %>%
                     column_to_rownames("regulon"),
                   scale="none",
                   cluster_cols=F,
                   cluster_rows=F,
                   fontsize_col=20,
                   main="Pruned regulons (based on SD of mean_AUC and %cells active)\ncolor scale: %cells with regulon active")
# dev.off()



# Regulons that are different between two methods
reg_diff <- unique(df_regulonsactive_pruned$regulon)[!unique(df_regulonsactive_pruned$regulon) %in% unique(df_regulonsactive_pruned2$regulon)]
pheatmap::pheatmap(df_regulonsactive_pruned %>%
                     filter(regulon %in% reg_diff) %>%
                     select(regulon, mean_auc, gep) %>%
                     pivot_wider(names_from=gep, values_from=mean_auc) %>%
                     arrange(factor(regulon, levels=htmp2_roworder[htmp2_roworder %in% regulon])) %>%
                     column_to_rownames("regulon"),
                   scale="row",
                   cluster_cols=F,
                   cluster_rows=F,
                   fontsize_col=20,
                   main="Pruned regulons  (active in >50% cells for at least 1 GEP)\ncolor scale: mean AUC across cells")

pheatmap::pheatmap(df_regulonsactive_pruned %>%
                     filter(regulon %in% reg_diff) %>%
                     select(regulon, propcells_bin_positive, gep) %>%
                     pivot_wider(names_from=gep, values_from=propcells_bin_positive) %>%
                     arrange(factor(regulon, levels=htmp2_roworder[htmp2_roworder %in% regulon])) %>%
                     column_to_rownames("regulon"),
                   scale="none",
                   cluster_cols=F,
                   cluster_rows=F,
                   fontsize_col=20,
                   main="Pruned regulons  (active in >50% cells for at least 1 GEP)\ncolor scale: %cells with regulon active")



## /end ####


#___________________________
## 3.4. Fourth analysis ####

## /end ####

