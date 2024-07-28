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




# **************
# FUNCTIONS ####
# **************

plot_DotPlot <- function(seurobj, group, features, scaling=T){
  # get plot data
  p <- Seurat::DotPlot(
    seurobj,
    group.by=group,
    features=features,
    scale=scaling
  )
  # plot in personalized way
  p <- ggplot(p$data, aes(x=id, y=features.plot, fill=avg.exp.scaled, size=pct.exp))+
    geom_point(color="black", shape=21)+
    # scale_fill_gradient2(low=scales::muted("blue"), high=scales::muted("red"), name="z-score\nnormalized\navg expression")+
    scale_size_continuous(range=c(0,6), limits=c(0,100), name="%cells\nexpressing\ngene")+
    theme_bw()+
    theme(
      axis.text.y=element_text(face="italic"),
      axis.text.x=element_text(angle=45, hjust=1)
    )+
    labs(y="", x="")
  # different color scale if scaled or not
  if(scaling==T){
    p <- p + scale_fill_gradient2(low=scales::muted("blue"), high=scales::muted("red"), name="z-score\naverage\nnormalized\nexpression")
  } else{
    p <- p + viridis::scale_fill_viridis(option="B", direction=-1, name="average\nnormalized\nexpression")
  }
  return(p)
}


# ******************
# 2. PLOT UMAPs ####
# ******************

## 2.1. HUMAN ####

seur.human <- subset(seur.human, Anno_level_3 != "Epi_GCM2") # remove parathyroid cells

# Create a new level of annotation
seur.human@meta.data$Anno_curated <- seur.human@meta.data$Anno_level_1
seur.human@meta.data$Anno_curated <- case_when(
  seur.human@meta.data$Anno_curated=="TEC" & seur.human@meta.data$Anno_level_3=="cTEC" ~ "cTEC",
  seur.human@meta.data$Anno_curated=="TEC" & seur.human@meta.data$Anno_level_3%in%c("mTEC","TEC(myo)", "TEC(neuro)") ~ "mTEC",
  # seur.human@meta.data$Anno_curated=="TEC" & seur.human@meta.data$Anno_level_3!="cTEC" & seur.human@meta.data$Anno_level_3!="mTEC" ~ "TEC_other",
  seur.human@meta.data$Anno_curated=="T"   & seur.human@meta.data$Anno_level_2=="DN" ~ "DN",
  seur.human@meta.data$Anno_curated=="T"   & seur.human@meta.data$Anno_level_2=="DP" ~ "DP",
  seur.human@meta.data$Anno_curated=="T"   & seur.human@meta.data$Anno_level_2=="SP" ~ "SP",
  seur.human@meta.data$Anno_curated=="Innate_T"   & seur.human@meta.data$Anno_level_3=="NKT" ~ "SP",
  seur.human@meta.data$Anno_curated=="Innate_T"   & seur.human@meta.data$Anno_level_3=="γδT" ~ "γδT",
  seur.human@meta.data$Anno_curated=="Endo"  ~ "Endothelial",
  seur.human@meta.data$Anno_curated=="Ery"   ~ "Erythrocyte",
  seur.human@meta.data$Anno_curated=="Mesen" ~ "Mesenchymal",
  seur.human@meta.data$Anno_curated=="Mgk"   ~ "Megakaryocyte",
  .default=seur.human@meta.data$Anno_curated
)
table(seur.human@meta.data$Anno_curated, useNA="ifany")
table(seur.human@meta.data[,c("Anno_curated", "Anno_level_1")], useNA="ifany")

# expect:
# B               5,082                          
# Endo              115
# Ery               644
# HSC               501
# Innate_lymphoid 2,176
# NKT               349 (Anno_level_1 Innate_T 2931)
# γδT             2,582 (Anno_level_1 Innate_T 2931)
# Mast              148
# Mesen          21,290
# Mgk                36
# Myeloid         4,801
# DN             42,474 (Anno_level_1 T 201,019)
# DP            108,418 (Anno_level_1 T 201,019)
# SP             50,127 (Anno_level_1 T 201,019)
# cTEC           10,156 (Anno_level_1 TEC 17,158)
# mTEC            6,448 (Anno_level_1 TEC 17,158)
# TEC_other         554 (Anno_level_1 TEC 17,158)

seur.human@meta.data$Anno_curated <- factor(seur.human@meta.data$Anno_curated,
                                            levels=c(
                                              "DN", "DP", "SP", "γδT", #"NKT",
                                              "cTEC", "mTEC", #"TEC_other",
                                              "HSC",
                                              "Innate_lymphoid",
                                              "B",
                                              "Myeloid",
                                              "Mast",
                                              "Erythrocyte",
                                              "Megakaryocyte",
                                              "Mesenchymal",
                                              "Endothelial"
                                            ))

# see how many cells there are per cluster
as.data.frame(table(seur.human$Anno_curated)) %>%
  mutate(totalcells=sum(Freq),
         percentcells=Freq*100/totalcells) %>%
  arrange(percentcells)
hu_clusters_abundant <- levels(seur.human@meta.data$Anno_curated)[!levels(seur.human@meta.data$Anno_curated) %in% c("Megakaryocyte",
                                                                                                                    "Endothelial",
                                                                                                                    "Mast",
                                                                                                                    "HSC",
                                                                                                                    "Erythrocyte")]
# remove clusters that contain less than 0.5% of all cells
# seur.human <- subset(seur.human, !Anno_curated %in% c("Megakaryocyte", "Endothelial", "Mast", "HSC", "Erythrocyte"))


celltypes_col <- c(
  "mTEC"            = "#CE3F37",
  "cTEC"            = "#8C6143",
  "TEC_other"= "#666666",
  "DN"              = "#FDCB89",
  "DP"              = "#9ecae1",
  "SP"              = "#7FC97F",
  "γδT"             = "#B6B1C9",
  # "NKT"             = "#B35C20",
  "Mesenchymal"     = "#FEE791",
  "Myeloid"         = "#F2F59A",
  "B"               = "#2171b5",
  "Endothelial"     = "#7D449D",
  "Erythrocyte"     = "#CD1588",
  "HSC"             = "#9BB5A4",
  "Innate_lymphoid" = "#EDBB99",
  "Mast"            = "#D1B3BB",
  "Megakaryocyte"   = "#9ABDA4"
  )

# UMAP
Idents(seur.human) <- "Anno_curated"
p1 <- ggrastr::rasterise(
  # ---
  SCpubr::do_DimPlot(seur.human,
                     reduction="UMAP",
                     idents.keep=hu_clusters_abundant,
                     na.value="grey90",
                     legend.position="right", legend.ncol=1, font.size=30, legend.icon.size=10)+
    scale_color_manual(values=celltypes_col),
  # ---
  layers="Point", dpi=300)
# ggsave(filename="./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots_final/hu_umap_clusters_highlight.pdf",
#        plot=p1,
#        device = cairo_pdf,
#        width=12, height=8)


# Plot CD1d expression
p2 <- ggrastr::rasterise(
  # ---
  SCpubr::do_FeaturePlot(seur.human, 
                       features = "CD1D", order = T,
                       # plot.title = "CD1D",
                       border.color = "black", border.size = 2,
                       reduction = "UMAP", pt.size = 1.2, legend.position = "right", font.size=30) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0.1, end = 0.85, direction = -1),
  # ---
  layers="Point", dpi=300)
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots_final/hu_umap_cd1d_orderT.pdf",
#        p2,
#        device = cairo_pdf,
#        width=11, height=8)

p1 | p2
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots_final/hu_umapcombined.pdf",
#        device = cairo_pdf,
#        width=20, height=8)

# Plot CD1d expression per cluster in violin plot
seur.human.abundant <- subset(seur.human, Anno_curated %in% hu_clusters_abundant)

ggrastr::rasterise(
  VlnPlot(seur.human.abundant,
          features = "CD1D", group.by = "Anno_curated", raster=F)+
    scale_fill_manual(values=celltypes_col)+
    labs(y="CD1D (normalized expression)", title="", x="")+
    theme(legend.position="none",
          axis.text.x=element_text(size=15),
          axis.title.y=element_text(size=15)),
  layers="Point", dpi=300)
# ggsave(filename="./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots_final/hu_vlnplt_cd1d.pdf",
#        device = cairo_pdf,
#        width=8, height=6)

# Plot CD1d expression per cluster in bubble plot
DotPlot(seur.human.abundant,
        features=rev(c("CD1D", "SLAMF1", "SLAMF6")),
        group.by="Anno_curated",
        cols=c("lightgrey", "darkred"),
        col.min=0,
        dot.scale=10
        )+
  coord_flip()+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  labs(x="", y="")
# ggsave(filename="./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots_final/hu_dotplot_cd1d.pdf",
#        device = cairo_pdf,
#        width=8, height=4.5)

# Plot for ANDCS poster
plot_DotPlot(seur.human.abundant, features = rev(c("CD1D", "SLAMF1", "SLAMF6")), group = "Anno_curated", scaling = F)
# ggsave("~/Desktop/Meyer-lab/Conferences/2024-07_ANDCS_Paris/fig4_cd1d_expression_human.pdf", width=10, height=5, units="cm")


# Plot all SFR expression
# Plot CD1d expression per cluster in bubble plot
DotPlot(seur.human.abundant,
        features=rev(c("CD1D",
                       "SLAMF1",
                       "LY9", # SLAMF3
                       "CD244", # SLAMF4
                       "CD84", # SLAMF5
                       "SLAMF6",
                       "SLAMF7",
                       "SLAMF8",
                       "SLAMF9"
                       )),
        group.by="Anno_curated",
        cols=c("lightgrey", "darkred"),
        col.min=0,
        dot.scale=10
)+
  coord_flip()+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  labs(x="", y="")
ggsave(filename="~/Desktop/Meyer-lab/Conferences/2024-02_CD1MR1-Hobart/presentation/hu_dotplot_slamfs.pdf",
       device = cairo_pdf,
       width=8, height=6)


## end 2.1. ####


## 2.2. MOUSE ####
seur.mouse <- subset(seur.mouse, age != "Rag1KO") # 34,073 cells

seur.mouse@meta.data$Anno_curated <- case_when(
  seur.mouse@meta.data$cell.types=="B"    ~ "B",
  seur.mouse@meta.data$cell.types %in% c("CD4+T", "CD8+T", "Treg", "αβT(entry)")  ~ "SP",
  seur.mouse@meta.data$cell.types %in% c("DN(P)", "DN(Q)")                        ~ "DN",
  seur.mouse@meta.data$cell.types %in% c("DP(P)", "DP(Q)")                        ~ "DP",
  seur.mouse@meta.data$cell.types=="Endo" ~ "Endothelial",
  seur.mouse@meta.data$cell.types=="Ery"  ~ "Erythrocyte",
  seur.mouse@meta.data$cell.types %in% c("Fb", "VSMC")                            ~ "Mesenchymal",
  seur.mouse@meta.data$cell.types %in% c("HSC", "NMP")                            ~ "HSC",
  seur.mouse@meta.data$cell.types %in% c("IELpA", "IELpB/NKT", "NK")              ~ "Innate_lymphoid",
  seur.mouse@meta.data$cell.types %in% c("DC1", "DC2", "aDC", "pDC")              ~ "Myeloid",
  seur.mouse@meta.data$cell.types %in% c("Mac", "Mono")                           ~ "Myeloid",
  seur.mouse@meta.data$cell.types %in% c("Epi_unknown", "TEC_early")              ~ "TEC_other",
  seur.mouse@meta.data$cell.types == "cTEC"  ~ "cTEC",
  seur.mouse@meta.data$cell.types == "mTEC"  ~ "mTEC",
  seur.mouse@meta.data$cell.types == "γδT"   ~ "γδT"
)

seur.mouse@meta.data$Anno_curated <- factor(seur.mouse@meta.data$Anno_curated,
                                            levels=c(
                                              "DN", "DP", "SP", "γδT",
                                              "cTEC", "mTEC", "TEC_other",
                                              "HSC",
                                              "Innate_lymphoid",
                                              "B",
                                              "Myeloid",
                                              "Erythrocyte",
                                              "Mesenchymal",
                                              "Endothelial"
                                            ))

# see how many cells there are per cluster
as.data.frame(table(seur.mouse$Anno_curated)) %>%
  mutate(totalcells=sum(Freq),
         percentcells=Freq*100/totalcells) %>%
  arrange(percentcells)
ms_clusters_abundant <- levels(seur.mouse@meta.data$Anno_curated)[!levels(seur.mouse@meta.data$Anno_curated) %in% c("Endothelial",
                                                                                                                    "HSC",
                                                                                                                    "Erythrocyte")]


# UMAP
p3 <- ggrastr::rasterise(
  # ---
  SCpubr::do_DimPlot(seur.mouse,
                     reduction="umap",
                     group.by="Anno_curated",
                     idents.keep=ms_clusters_abundant,
                     na.value="grey90",
                     legend.position="right", legend.ncol=1, font.size=30, legend.icon.size=10)+
    scale_color_manual(values=celltypes_col),
  # ---
  layers="Point", dpi=300)

p4 <- ggrastr::rasterise(
  # ---
  SCpubr::do_FeaturePlot(seur.mouse, 
                       features = "Cd1d1", order = T,
                       border.color = "black", border.size = 2,
                       reduction = "umap", pt.size = 1.2, legend.position = "right", font.size=30) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0.1, end = 0.85, direction = -1),
  # ---
  layers="Point", dpi=300)

p3 | p4
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots_final/ms_umapcombined.pdf",
#        device = cairo_pdf,
#        width=20, height=8)


# plot Cd1d1 expression in all clusters
seur.mouse.abundant <- subset(seur.mouse, Anno_curated %in% ms_clusters_abundant)

ggrastr::rasterise(
  VlnPlot(seur.mouse.abundant, features = "Cd1d1", group.by = "Anno_curated", raster=F)+
    scale_fill_manual(values=celltypes_col)+
    labs(y="Cd1d1 (normalized expression)", title="", x="")+
    theme(legend.position="none",
          axis.text.x=element_text(size=15),
          axis.title.y=element_text(size=15)),
  layers="Point", dpi=300)
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots_final/ms_vlnplt_cd1d1.pdf",
#        device=cairo_pdf,
#        width=8, height=6)


# Plot CD1d expression per cluster in bubble plot
DotPlot(seur.mouse.abundant,
        features=rev(c("Cd1d1", "Slamf1", "Slamf6")),
        group.by="Anno_curated",
        cols=c("lightgrey", "darkred"),
        col.min = 0,
        dot.scale=10
)+
  coord_flip()+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  labs(x="", y="")
# ggsave(filename="./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots_final/ms_dotplot_cd1d1.pdf",
#        device = cairo_pdf,
#        width=8, height=4.5)

# Plot for ANDCS poster
plot_DotPlot(seur.mouse.abundant, features = rev(c("Cd1d1", "Slamf1", "Slamf6")), group = "Anno_curated", scaling = F)
# ggsave("~/Desktop/Meyer-lab/Conferences/2024-07_ANDCS_Paris/fig4_cd1d_expression_mouse.pdf", width=10, height=5, units="cm")




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
flowdf <- read.csv("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/flow_cytometry_plots/percentages_human.csv", header=T)

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
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/flow_boxplot_human_thymocyte.jpeg", width=4, height=5)


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
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/flow_boxplot_human_tec.jpeg", width=3, height=5)


# MOUSE
flowdf.ms <- read.csv("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/flow_cytometry_plots/percentages_mouse.csv", header=T)

ggplot(flowdf.ms,
       aes(x=factor(celltype, levels=c("CD45+", "mTEC", "cTEC")), y=percent_CD1d))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=3, width=0.15)+
  # scale_color_manual(values=sample_cols)+
  labs(x="", y="% cells CD1d+", title="%CD1d+ cells")+
  theme_cowplot()+
  ylim(c(0,100))+
  theme(axis.text.x=element_text(size=20, angle=45, hjust=1),
        axis.title.y=element_text(size=20)) |
ggplot(flowdf.ms,
       aes(x=factor(celltype, levels=c("CD45+", "mTEC", "cTEC")), y=median_fluorescence_CD1d_pos_cells))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=3, width=0.15)+
  # scale_color_manual(values=sample_cols)+
  labs(x="", y="MFI (CD1d+ cells)", title="Median Fluorescence")+
  theme_cowplot()+
  ylim(c(0,1000))+
  theme(axis.text.x=element_text(size=20, angle=45, hjust=1),
        axis.title.y=element_text(size=20))
ggsave("./scripts-in-progress/human-thymus/HumanThymus_22_ParkData_CD1dMR1/plots/flow_boxplot_mouse.jpeg", width=7, height=5)
