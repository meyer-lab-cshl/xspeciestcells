# Purpose: plot CD1d and MR1 expression in Park dataset
# Author: Salomé Carcy
# Date: August 2023


# **************
# 1. IMPORT ####
# **************

library(Seurat)
library(ggplot2)
library(tidyverse)
# library(dplyr)
library(cowplot)
library(RColorBrewer)
library(scico)
source("./scripts-final/colors_universal.R")


# Import data
seur.human <- readRDS("~/Projects/HumanThymusProject/data/human-thymus/HumanData_16_ParkIntegration/park_full_seu_gene_names.rds")
seur_metadata <- read.csv("~/Projects/HumanThymusProject/data/human-thymus/HumanData_16_ParkIntegration/park_fullmetadata.csv")

seur.mouse <- readRDS("~/Projects/HumanThymusProject/data/raw_data/mouse_data/thymus_Park/park_seu_mouse.rds")

# Incorporate metadata in seurat object
table(seur_metadata$index == rownames(seur.human@meta.data), useNA="ifany")
seur.human@meta.data <- cbind(seur.human@meta.data, seur_metadata)



# ******************
# 2. PLOT UMAPs ####
# ******************

celltypes_col <- c("#7FC97F", "#ABB6BA", "#D7B5B4", "#FDC687", "#FEF295", "#9BB5A4", "#5C56A6", "#DD0C83", "#D23740", "#A45E2E")
names(celltypes_col) <- c("T_naive", "DP", "CD8αα(I)", "DN", "Mono", "B_naive", "γδT", "B_memory", "mTEC", "cTEC")


## 2.1. HUMAN ####

# Clusters
# colnames(seur.human@meta.data)
table(seur.human@meta.data$Anno_level_3, useNA="ifany")
seur.human@meta.data[seur.human@meta.data$Anno_level_3=="TEC(neuro)", "Anno_level_3"] <- "mTEC"
seur.human@meta.data[seur.human@meta.data$Anno_level_3=="TEC(myo)", "Anno_level_3"] <- "mTEC"
table(seur.human@meta.data$Anno_level_3=="mTEC", useNA="ifany") # should have 6928 mTEC


Idents(seur.human) <- "Anno_level_3"
p1 <- ggrastr::rasterise(
  # ---
  SCpubr::do_DimPlot(seur.human,
                     reduction="UMAP",
                     # group.by="Anno_level_fig1",
                     # group.by="Anno_level_3",
                     # cells.highlight = rownames(seur.human@meta.data[seur.human@meta.data$Anno_level_3 %in% c("cTEC", "DN", "DP"),]),
                     idents.keep=c("cTEC", "mTEC", "DN", "DP", "CD8αα(I)", "Mono", "T_naive", "γδT", "B_naive", "B_memory"),
                     na.value="grey90",
                     legend.position="right", legend.ncol=1, font.size=30, legend.icon.size=10)+
    # scale_color_manual(values=grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(11))
    scale_color_manual(values=celltypes_col),
  # ---
  layers="Point", dpi=300)
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/umap_clusters_level3_highlight.jpeg", width=10, height=8)


# Plot CD1d expression
p2 <- SCpubr::do_FeaturePlot(seur.human, 
                       features = "CD1D", order = T,
                       # plot.title = "CD1D",
                       border.color = "black", border.size = 2,
                       reduction = "UMAP", pt.size = 1.2, legend.position = "right", font.size=30) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0.1, end = 0.85, direction = -1)
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/umap_cd1d_orderT.jpeg", width=10, height=8)


p1 | p2
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/umap_combined.jpeg", width=20, height=8)

# Plot CD1d expression per cluster in violin plot
celltypes_col2 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(30)
celltypes_col2 <- celltypes_col2[!celltypes_col2 %in% celltypes_col]
celltypes_col2 <- celltypes_col2[1:28]
names(celltypes_col2) <- unique(seur.human$Anno_level_3)[!unique(seur.human$Anno_level_3) %in% names(celltypes_col)]
celltypes_col_human <- c(celltypes_col, celltypes_col2)

ggrastr::rasterise(
  VlnPlot(seur.human, features = "CD1D", group.by = "Anno_level_3", raster=F)+
    scale_fill_manual(values=celltypes_col_human)+
    labs(y="CD1D (normalized expression)", title="", x="")+
    theme(legend.position="none",
          axis.text.x=element_text(size=15),
          axis.title.y=element_text(size=15)),
  layers="Point", dpi=300)
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/vlnplt_cd1d_level3.jpeg", width=12, height=6)




## 2.2. MOUSE ####
seur.mouse@meta.data <- seur.mouse@meta.data %>%
  mutate(cell.types=case_when(
    cell.types %in% c("DN(Q)", "DN(P)") ~ "DN",
    cell.types %in% c("DP(Q)", "DP(P)") ~ "DP",
    TRUE ~ cell.types
  ))
seur.mouse <- subset(seur.mouse, age != "Rag1KO") # 34,073 cells

celltypes_col_mouse_umap <- c("#FDC687", "#ABB6BA", "#c7e9c0", "#7FC97F", "#D23740", "#A45E2E")
names(celltypes_col_mouse_umap) <- c("DN", "DP","αβT(entry)", "CD4+T", "mTEC", "cTEC")

p3 <- ggrastr::rasterise(
  # ---
  SCpubr::do_DimPlot(seur.mouse,
                     reduction="umap",
                     group.by="cell.types",
                     idents.keep=c("DN", "DP","αβT(entry)", "CD4+T", "mTEC", "cTEC"),
                     na.value="grey90",
                     legend.position="right", legend.ncol=1, font.size=30, legend.icon.size=10)+
    scale_color_manual(values=celltypes_col_mouse_umap),
  # ---
  layers="Point", dpi=300)

p4 <- SCpubr::do_FeaturePlot(seur.mouse, 
                       features = "Cd1d1", order = T,
                       border.color = "black", border.size = 2,
                       reduction = "umap", pt.size = 1.2, legend.position = "right", font.size=30) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0.1, end = 0.85, direction = -1)

p3 | p4
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/umapmouse_combined.jpeg", width=20, height=8)

# plot Cd1d1 expression in all clusters
celltypes_col_mouse2 <- celltypes_col_human[names(celltypes_col_human) %in% unique(seur.mouse$cell.types)]
unique(seur.mouse$cell.types)[!unique(seur.mouse$cell.types) %in% names(celltypes_col_mouse2)]
celltypes_col_mouse1 <- c("αβT(entry)"="#c7e9c0",
                          "CD4+T"="#7FC97F",
                          "IELpA"="#77479F",
                          "CD8+T"="#006d2c",
                          "IELpB/NKT"="#FDCA89",
                          "B"="#9BB5A4",
                          "Fb"="#EC0877",
                          "TEC_early"="#cc4c02",
                          "HSC"="#5380AC",
                          "Epi_unknown"="#7B6352")
celltypes_col_mouse <- c(celltypes_col_mouse1, celltypes_col_mouse2)
table(unique(seur.mouse$cell.types) %in% names(celltypes_col_mouse), useNA="ifany")

VlnPlot(seur.mouse, features = "Cd1d1", group.by = "cell.types", raster=F)+
  scale_fill_manual(values=celltypes_col_mouse)+
  labs(y="Cd1d1 (normalized expression)", title="", x="")+
  theme(legend.position="none",
        axis.text.x=element_text(size=15),
        axis.title.y=element_text(size=15))
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/vlnplt_cd1d1_mouse.jpeg", width=10, height=6)

# Plot % cTECs with Cd1d expression > 0 in mouse vs human
ctec_hu <- seur.human@meta.data[seur.human@meta.data$Anno_level_3 %in% c("cTEC", "mTEC", "DP"), c("Anno_level_3", "donor_id")] # 10,156 cells
ctec_ms <- seur.mouse@meta.data[seur.mouse@meta.data$cell.types %in% c("cTEC", "mTEC", "DP"), c("cell.types", "sample_ID")] # 2,917 cells

ctec_hu_count <- seur.human@assays$RNA@data[,rownames(ctec_hu)]
ctec_ms_count <- seur.mouse@assays$RNA@data[,rownames(ctec_ms)]

table(colnames(ctec_hu_count)==rownames(ctec_hu), useNA="ifany")
table(colnames(ctec_ms_count)==rownames(ctec_ms), useNA="ifany")

ctec_hu$countnorm <- ctec_hu_count["CD1D",]
ctec_ms$countnorm <- ctec_ms_count["Cd1d1",]

ctec_df <- rbind(
  ctec_hu %>%
    dplyr::rename(celltype=Anno_level_3) %>%
    mutate(donor_id=factor(donor_id, levels=unique(ctec_hu$donor_id))) %>%
    group_by(celltype, donor_id) %>%
    summarise(ncells=n(),
              ncells_cd1d_pos=sum(countnorm>0)) %>%
    ungroup() %>%
    filter(ncells>10) %>%
    dplyr::rename(sample_ID=donor_id) %>%
    mutate(percent_cd1d_pos=ncells_cd1d_pos*100/ncells,
           species="human"),
  ctec_ms %>%
    dplyr::rename(celltype=cell.types) %>%
    mutate(sample_ID=factor(sample_ID, levels=unique(ctec_ms$sample_ID))) %>%
    group_by(celltype, sample_ID) %>%
    summarise(ncells=n(),
              ncells_cd1d_pos=sum(countnorm>0)) %>%
    ungroup() %>%
    filter(ncells>10) %>%
    mutate(percent_cd1d_pos=ncells_cd1d_pos*100/ncells,
           species="mouse")
)

ggplot(ctec_df, aes(x=species, y=percent_cd1d_pos))+
  facet_wrap(~celltype)+
  geom_boxplot()+
  geom_jitter(width=0.1)+
  theme_cowplot()+
  labs(x="", y="% cTECs with CD1D expression > 0")


ggplot(rbind(
  ctec_hu %>% mutate(species="human") %>% dplyr::rename(celltype=Anno_level_3, sample_ID=donor_id),
  ctec_ms %>% mutate(species="mouse") %>% dplyr::rename(celltype=cell.types)
))+
  facet_grid(species~celltype)+
  geom_histogram(aes(x=countnorm, color=species))+
  scale_y_continuous(trans="log2")




# *****************************
# 3. PLOT FLOW PERCENTAGES ####
# *****************************

# import csv file
flowdf <- read.csv("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/flow_cytometry_plots/percentages.csv", header=T)

# plot thymocytes
sample_cols <- RColorBrewer::brewer.pal(7, "Set2")
names(sample_cols) <- unique(flowdf$sample)

ggplot(flowdf %>% filter(celltype=="thymocyte"),
       aes(x=factor(cellsubtype, levels=c("DN", "DP", "CD4SP", "CD8SP")), y=percent_CD1d))+
  # geom_boxplot(aes(fill=factor(cellsubtype, levels=c("DN", "DP", "CD4SP", "CD8SP"))))+
  # scale_fill_manual(values=c("#FDC687", "#ABB6BA", "#7FC97F" , "darkgreen"), name="")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=3, width=0.15)+
  # scale_color_manual(values=sample_cols)+
  labs(x="", y="% cells CD1d+")+
  theme_cowplot()+
  ylim(c(0,100))+
  theme(axis.text.x=element_text(size=20, angle=45, hjust=1),
        axis.title.y=element_text(size=20))
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/flow_boxplot_thymocyte.jpeg", width=4, height=5)


# plot TECs
ggplot(flowdf %>% filter(celltype=="tec"),
       aes(x=factor(cellsubtype, levels=c("mTEC", "cTEC")), y=percent_CD1d))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=3, width=0.15)+
  # scale_color_manual(values=sample_cols)+
  labs(x="", y="% cells CD1d+")+
  theme_cowplot()+
  ylim(c(0,100))+
  theme(axis.text.x=element_text(size=20),
        axis.title.y=element_text(size=20))
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/flow_boxplot_tec.jpeg", width=3, height=5)

