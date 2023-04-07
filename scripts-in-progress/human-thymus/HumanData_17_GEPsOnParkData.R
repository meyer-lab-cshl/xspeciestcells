# Purpose: run GEPs on Park dataset
# Author: Salomé Carcy
# Date: March 2023


# **************
# 1. IMPORT ####
# **************

# Import librairies
library(Seurat)
library(cowplot)
library(ggplot)

path <- "~/Projects/HumanThymusProject/data/human-thymus/HumanData_17_GEPsOnParkData/"

# Import datasets
seur.park.thymocytes <- readRDS("~/Projects/HumanThymusProject/data/human-thymus/HumanData_16_ParkIntegration/park_thymus_seu_gene_names.rds")
cols_park <- c("DN(early)" = "#78c679",
               "DN(P)" = "#41ab5d",
               "DN(Q)" = "#238443",
               "γδT" = "#92c051",
               "DP(P)" = "#b75347",
               "DP(Q)" = "#d8443c",
               "αβT(entry)" = "#e09351",
               "CD8+T"= "#5a97c1",
               "CD8αα(I)" = "#421401",
               "CD8αα(II)" = "#0a2e57",
               "CD4+T"= "gold",
               "T(agonist)" = "#9f5691",
               "Treg(diff)" = "#9f5691",
               "Treg" = "blueviolet",
               "Th17" = "#a40000",
               "NKT" = "#72bcd5")
DimPlot(seur.park.thymocytes, group.by="cell_type", label=T)+
  scale_color_manual(values=cols_park)

seur <- readRDS("~/Projects/HumanThymusProject/data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
cols_integrated <- c("0" = "#f4c40f", "1" = "#b75347", "2" = "#d8443c", "3" = "#e09351", "4" = "#2b9b81", 
                     "5" = "#421401", "6" = "#92c051", "7" = "#9f5691", "8" = "#17154f", "9" = "#74c8c3", 
                     "10" = "#5a97c1", "11" = "gold", "12" = "#a40000", "13" = "#72bcd5", "14" = "grey50",
                     "15" = "orange", "16" = "blueviolet", "17" = "#0a2e57", "18" = "#bdbdbd")
DimPlot(seur, group.by = "new_clusters", label=T, repel=T, reduction="UMAP_50")+
  scale_color_manual(values=cols_integrated)


# Import GEPs
gep_topgenes <- read.csv(file.path(path, "TopGenes_perGEP_ALLcells_2023-03-13.csv"), row.names=1)
dim(gep_topgenes) # 19 GEPs and 50 genes per GEP




# ************************
# 2. RUN MODULE SCORE ####
# ************************

# Make list of GEPs
gep_list <- list()
for (i in 1:ncol(gep_topgenes)){
  print(paste0("GEP", i))
  gep_allgenes <- gep_topgenes[1:100,i]
  gep_genesinpark <- gep_allgenes[gep_allgenes %in% rownames(seur.park.thymocytes)]
  print(length(gep_genesinpark))
  gep_list[[colnames(gep_topgenes)[i]]] <- gep_genesinpark
}
# names(gep_list) # sanity check

# Check for each GEP how many genes are not present in park seurat object
# gep_list_genes_notinpark <- lapply(gep_list, function(x) x[!x %in% rownames(seur.park.thymocytes)])
# lengths(gep_list_genes_notinpark)

# Run module score on park data and our data
# test1 <- AddModuleScore(seur.park.thymocytes, name = "GEP", features=list(gep_list[[1]])) # this is only for sanity check
seur.park.gep <- AddModuleScore(seur.park.thymocytes, name = "GEP", features=gep_list)
seur.geps <- AddModuleScore(seur, name = "GEP", features=gep_list)
# sanity check
# table(test1$GEP1==seur.park.gep$GEP1) # all same




# *******************************
# 3. VISUALIZE GEPs ON UMAPs ####
# *******************************

# Park
SCpubr::do_FeaturePlot(seur.park.gep, reduction="UMAP", features=paste0("GEP", 1:19), ncol=7,
                       viridis_color_map = "B")
# ggsave(file.path(path, "park_cNMF_geps.jpeg"), width=40, height=20)


# Our data
SCpubr::do_FeaturePlot(seur.geps, reduction="UMAP_50", features=paste0("GEP", 1:19), ncol=7,
                       viridis_color_map = "B")
# ggsave(file.path(path, "gapin_cNMF_geps_parkgenes.jpeg"), width=40, height=20)




# *********************************
# 3. VISUALIZE GEPs ON BATCHES ####
# *********************************

# Our data
df <- seur.geps@meta.data %>%
  rownames_to_column("cell") %>%
  as_tibble() %>%
  select(cell, Batch, Donor, Tissue, new_clusters, paste0("GEP", 1:19)) %>%
  pivot_longer(cols=paste0("GEP", 1:19), names_to="GEP", values_to="score")

ggplot(df, aes(x=Batch, y=score))+
  facet_wrap(~factor(GEP, levels=paste0("GEP", 1:19)), ncol=7)+
  geom_violin()+
  theme_cowplot()
# ggsave(file.path(path, "gapin_cNMF_geps_parkgenes_violin-per-batch.jpeg"), width=10, height=8)

ggplot(df, aes(x=factor(new_clusters, levels=0:17), y=score, fill=new_clusters))+
  facet_wrap(~factor(GEP, levels=paste0("GEP", 1:19)), ncol=7, scales="free_x")+
  # scale_fill_manual(values=c("PBMC"="#a40000", "Thymus"="#72bcd5"))+
  scale_fill_manual(values=cols_integrated)+
  geom_violin()+
  xlab("")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1, size=7),
        legend.position="none")
ggsave(file.path(path, "gapin_cNMF_geps_parkgenes_violin-per-cluster.jpeg"), width=15, height=8)


# Park data: visualize GEPs by cluster (cell_type)
df.park <- seur.park.gep@meta.data %>%
  rownames_to_column("cell") %>%
  as_tibble() %>%
  select(cell, sample_id, cell_type, paste0("GEP", 1:19)) %>%
  pivot_longer(cols=paste0("GEP", 1:19), names_to="GEP", values_to="score")

# ggplot(df.park, aes(x=sample_id, y=score))+
#   facet_wrap(~factor(GEP, levels=paste0("GEP", 1:19)), ncol=7)+
#   geom_violin()+
#   theme_cowplot()+
#   theme(axis.text.x=element_text(angle=90, hjust=1, size=5))
# ggsave(file.path(path, "park_cNMF_geps_violin-per-batch.jpeg"), width=20, height=8)
order_park <- c("DN(early)", "DN(P)", "DN(Q)", "DP(P)", "DP(Q)", "αβT(entry)", "T(agonist)", "CD8αα(I)", "CD8αα(II)",
                "γδT", "Treg(diff)", "Treg", "CD4+T", "CD8+T", "Th17", "NKT")

ggplot(df.park, aes(x=factor(cell_type, levels=order_park), y=score, fill=cell_type))+
  facet_wrap(~factor(GEP, levels=paste0("GEP", 1:19)), ncol=7)+
  scale_fill_manual(values=cols_park)+
  geom_violin()+
  xlab("")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1, size=7), legend.position="none")
# ggsave(file.path(path, "park_cNMF_geps_violin-per-cluster.jpeg"), width=15, height=8)

