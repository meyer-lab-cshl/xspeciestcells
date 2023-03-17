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


# Plot clusters 3 and 11 CD4/NKT
seur@meta.data$new_group <- paste(seur@meta.data$group.ident, seur@meta.data$new_clusters, sep="_")
Idents(seur) <- "new_group"

DimPlot(seur, #label=T,
        reduction = "UMAP_50", #group.by = "new_clusters",
        cells.highlight= list("3"=WhichCells(seur, idents="CD4_Thymus_3"),
                              "11"=WhichCells(seur, idents="CD4_Thymus_11")),
        cols.highlight = c("3"="#e09351", "11"="gold"),
        pt.size=0.2,
        sizes.highlight = 1,
        cols= "#bdbdbd" )+
  labs(x="UMAP 1", y="UMAP 2", title="CD4+ cells")+
  theme_cowplot()+
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(hjust=0, size=30),
        plot.title = element_text(hjust=0.5, size=40),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position="none") |
DimPlot(seur, #label=T,
          reduction = "UMAP_50", #group.by = "new_clusters",
          cells.highlight= list("3"=WhichCells(seur, idents="NKT_Thymus_3"),
                                "11"=WhichCells(seur, idents="NKT_Thymus_11")),
          cols.highlight = c("3"="#e09351", "11"="gold"),
          pt.size=0.2,
          sizes.highlight = 1,
          cols= "#bdbdbd" )+
  labs(x="UMAP 1", y="UMAP 2", title="iNKT cells")+
  theme_cowplot()+
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(hjust=0, size=30),
        plot.title = element_text(hjust=0.5, size=40),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position="none")
ggsave("~/Projects/HumanThymusProject/data/human-thymus/HumanData_14_CCR9-CCR7-plots/cd4_nkt_comparison.jpeg", width=13, height=7)


# Plot clusters 9 and 10 CD8/MAIT
DimPlot(seur, #label=T,
        reduction = "UMAP_50", #group.by = "new_clusters",
        cells.highlight= list("9"=WhichCells(seur, idents="CD8_Thymus_9"),
                              "10"=WhichCells(seur, idents="CD8_Thymus_10")),
        cols.highlight = c("9"="#74c8c3", "10"="#5a97c1"),
        pt.size=0.2,
        sizes.highlight = 1,
        cols= "#bdbdbd" )+
  labs(x="UMAP 1", y="UMAP 2", title="CD8+ cells")+
  theme_cowplot()+
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(hjust=0, size=30),
        plot.title = element_text(hjust=0.5, size=40),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position="none") |
  DimPlot(seur, #label=T,
          reduction = "UMAP_50", #group.by = "new_clusters",
          cells.highlight= list("9"=WhichCells(seur, idents="MAIT_Thymus_9"),
                                "10"=WhichCells(seur, idents="MAIT_Thymus_10")),
          cols.highlight = c("9"="#74c8c3", "10"="#5a97c1"),
          pt.size=0.2,
          sizes.highlight = 1,
          cols= "#bdbdbd" )+
  labs(x="UMAP 1", y="UMAP 2", title="MAIT cells")+
  theme_cowplot()+
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(hjust=0, size=30),
        plot.title = element_text(hjust=0.5, size=40),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position="none")
ggsave("~/Projects/HumanThymusProject/data/human-thymus/HumanData_14_CCR9-CCR7-plots/cd8_mait_comparison.jpeg", width=13, height=7)

