# Purpose: UMAP plot of CCR9 and CCR7
# Author: Salomé
# Date: 02.27.23

# Import librairies and seurat object
library(Seurat)
library(ggplot2)
library(cowplot)
library(SCpubr)
library(ggpubr)

seur <- readRDS("~/Projects/HumanThymusProject/data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
DimPlot(seur, reduction="UMAP_50", group.by="new_clusters")


# Plot CCR9
for(group in c("CD4_Thymus", "CD8_Thymus", "NKT_Thymus", "MAIT_Thymus", "GD_Thymus")){
  SCpubr::do_FeaturePlot(seur, 
                         cells.highlight= colnames(seur[,seur$group.ident=="CD4_Thymus"]),
                         order=T,
                         viridis_color_map="A",
                         slot="data",
                         plot.title=group,
                         legend.position="right",
                         max.cutoff=5,
                         features = c("CCR9"))
  ggsave(paste0("~/Projects/HumanThymusProject/data/human-thymus/HumanData_14_CCR9-CCR7-plots/ccr9_", group, ".jpeg"), width=8, height=7)
}

SCpubr::do_FeaturePlot(seur, 
                       # cells.highlight= colnames(seur[,seur$Tissue=="Thymus"]),
                       split.by="group.ident",
                       order=T,
                       viridis_color_map="A",
                       slot="data",
                       legend.position="right",
                       ncol=5,
                       features = c("CCR9"))
ggsave("~/Projects/HumanThymusProject/data/human-thymus/HumanData_14_CCR9-CCR7-plots/ccr9_thymocytes.jpeg", width=35, height=14)


# Plot CCR7
for(group in c("CD4_Thymus", "CD8_Thymus", "NKT_Thymus", "MAIT_Thymus", "GD_Thymus")){
  SCpubr::do_FeaturePlot(seur, 
                         cells.highlight= colnames(seur[,seur$group.ident==group]),
                         order=T,
                         viridis_color_map="A",
                         slot="data",
                         plot.title=group,
                         legend.position="right",
                         features = c("CCR7"))
  ggsave(paste0("~/Projects/HumanThymusProject/data/human-thymus/HumanData_14_CCR9-CCR7-plots/ccr7_", group, ".jpeg"), width=8, height=7)
}

SCpubr::do_FeaturePlot(seur, 
                       # cells.highlight= colnames(seur[,seur$Tissue=="Thymus"]),
                       split.by="group.ident",
                       order=T,
                       viridis_color_map="A",
                       slot="data",
                       legend.position="right",
                       ncol=5,
                       features = c("CCR7"))
ggsave("~/Projects/HumanThymusProject/data/human-thymus/HumanData_14_CCR9-CCR7-plots/ccr7_thymocytes.jpeg", width=35, height=14)



# Plot group.ident
test <- c("#74c476", "#f768a1", "#084594", "#9ecae1","#9e9ac8", rep("black", 5))
names(test) <- c("CD4_Thymus", "CD8_Thymus", "GD_Thymus", "MAIT_Thymus", "NKT_Thymus",
                 "CD4_PBMC", "CD8_PBMC", "GD_PBMC", "MAIT_PBMC", "NKT_PBMC")

SCpubr::do_DimPlot(seur, 
                   split.by = "group.ident", 
                   ncol = 5,
                   idents.keep = c("CD4_Thymus", "CD8_Thymus", "GD_Thymus", "MAIT_Thymus", "NKT_Thymus"),
                   colors.use = test,
                   # legend.position = "none",
                   font.size = 24)
ggsave("~/Projects/HumanThymusProject/data/human-thymus/HumanData_14_CCR9-CCR7-plots/thymocytes_highlighted.jpeg", width=25, height=7)
