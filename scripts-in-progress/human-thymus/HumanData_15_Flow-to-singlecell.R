# Purpose: reproduce flow data on single cell
# Author: Salomé Carcy
# Date: March 2023


# **************
# 1. IMPORT ####
# **************

# Import librairies
library(ggplot2)
library(ggpointdensity)
library(cowplot)
library(tidyverse)
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(pals)


# Import integrated seurat object
seur <- readRDS("~/Projects/HumanThymusProject/data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS") # removed low-quality cells
cols_integrated <- c("0" = "#f4c40f", "1" = "#b75347", "2" = "#d8443c", "3" = "#e09351", "4" = "#2b9b81", 
                     "5" = "#421401", "6" = "#92c051", "7" = "#9f5691", "8" = "#17154f", "9" = "#74c8c3", 
                     "10" = "#5a97c1", "11" = "gold", "12" = "#a40000", "13" = "#72bcd5", "14" = "grey50",
                     "15" = "orange", "16" = "blueviolet", "17" = "#0a2e57", "18" = "#bdbdbd")
SCpubr::do_DimPlot(seur, 
                   group.by = "new_clusters",
                   label=T,
                   label.color="black",
                   legend.position = "none",
                   repel=T,
                   colors.use=cols_integrated,
                   font.size = 24)

path.plots <- "~/Projects/HumanThymusProject/data/human-thymus/HumanData_15_Flow-to-singlecell"


# *******************
# 2. THYMIC iNKT ####
# *******************

## 2.1. Get thymic iNKT data ####

# Import thymic iNKT object
seur.thynkt <- readRDS("~/Projects/HumanThymusProject/data/human-thymus/HumanData_12_AnalysisByLineage/thymus.nkt_02_28_23.RDS")

# Get the counts
counts <- GetAssayData(seur.thynkt, slot="data")
counts.df <- counts[c("CD4", "CD8A", "EOMES", "KLRB1", "ZBTB16", "GZMK", "CCR7", "CCR9", "SELL"),] %>%
  data.frame() %>%
  t()
counts.df <- data.frame(counts.df) # 2535 cells

# Sanity checks
table(is.na(counts.df$CD4))
table(is.na(counts.df$CD8A))
table(is.na(counts.df$EOMES))
table(is.na(counts.df$KLRB1))
table(is.na(counts.df$ZBTB16))
table(is.na(counts.df$GZMK))
table(is.na(counts.df$CCR7))
table(is.na(counts.df$CCR9))



## 2.2. Flow with CCR9 ####

# Get all the conditions
counts.df <- counts.df %>%
  mutate("tag_cd4cd8" = ifelse(CD4 > 1 & CD8A > 1, "DP", ifelse(CD4 > 1 & CD8A < 1, "CD4+", "DN/CD8+"))) %>%
  mutate("tag_cd161eomes" = ifelse(EOMES < 1 & KLRB1 < 1, "EOMES-CD161-",
                                   ifelse(EOMES < 1 & KLRB1 >= 1, "EOMES-CD161+",
                                          ifelse(EOMES >= 1 & KLRB1 < 1, "EOMES+CD161-",
                                                 ifelse(EOMES >= 1 & KLRB1 >=1, "EOMES+CD161+", "NA"))))) %>%
  mutate("tag_GzmkPlzf" = ifelse(ZBTB16 < 1 & GZMK < 1, "ZBTB16-GZMK-",
                                 ifelse(ZBTB16 < 1 & GZMK >= 1, "ZBTB16-GZMK+",
                                        ifelse(ZBTB16 >= 1 & GZMK < 1, "ZBTB16+GZMK-",
                                               ifelse(ZBTB16 >= 1 & GZMK >=1, "ZBTB16+GZMK+", "NA"))))) %>%
  mutate("tag_Ccr7Ccr9" = ifelse(CCR7 < 1 & CCR9 < 1, "CCR7-CCR9-",
                                 ifelse(CCR7 < 1 & CCR9 >= 1, "CCR7-CCR9+",
                                        ifelse(CCR7 >= 1 & CCR9 < 1, "CCR7+CCR9-",
                                               ifelse(CCR7 >= 1 & CCR9 >=1, "CCR7+CCR9+", "NA"))))) %>%
  mutate("window1"= tag_cd4cd8,
         "window2"= paste(tag_cd4cd8, tag_cd161eomes, sep=" > "),
         "window3"= paste(window2, tag_GzmkPlzf, sep=" > "),
         "window4" = paste(window3, tag_Ccr7Ccr9, sep=" > "))


# CD4 vs CD8
ggplot(data.frame(counts.df), aes(x=CD8A, y=CD4)) +
  # geom_point(aes(color=tag_cd4cd8))+
  geom_pointdensity() +
  scale_color_viridis_c()+
  theme_cowplot()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
# ggsave("~/Downloads/cd4cd8.jpeg", width=6, height=5)


# CD161 vs EOMES
ggplot(counts.df %>% filter(tag_cd4cd8=="CD4+"), aes(x=EOMES, y=KLRB1)) +
  # geom_point(aes(color=tag_cd161eomes))+
  geom_pointdensity() +
  scale_color_viridis_c()+
  theme_cowplot()+
  labs(y="CD161", title="CD4+") +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
# ggsave("~/Downloads/EomesKlrb1_cd4.jpeg", width=6, height=5)
ggplot(counts.df %>% filter(tag_cd4cd8=="DN/CD8+"), aes(x=EOMES, y=KLRB1)) +
  # geom_point(aes(color=tag_cd161eomes))+
  geom_pointdensity() +
  scale_color_viridis_c()+
  theme_cowplot()+
  labs(y="CD161", title="DN/CD8+") +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
# ggsave("~/Downloads/EomesKlrb1_cd8-dn.jpeg", width=6, height=5)


# GZMK vs PLZF
table(counts.df$window2)
ggplot(counts.df %>% filter(window2=="CD4+ > EOMES-CD161-"), aes(x=ZBTB16, y=GZMK)) +
  # geom_point(aes(color=tag_GzmkPlzf))+
  geom_pointdensity() +
  scale_color_viridis_c()+
  theme_cowplot()+
  labs(x="PLZF", title="CD4+ | EOMES-CD161- (819 cells)") +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
# ggsave("~/Downloads/PlzfGzmk_Cd4.jpeg", width=6, height=5)
ggplot(counts.df %>% filter(window2=="DN/CD8+ > EOMES-CD161-"), aes(x=ZBTB16, y=GZMK)) +
  # geom_point(aes(color=tag_GzmkPlzf))+
  geom_pointdensity() +
  scale_color_viridis_c()+
  theme_cowplot()+
  labs(x="PLZF", title="DN/CD8+ | EOMES-CD161- (1115 cells)") +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
# ggsave("~/Downloads/PlzfGzmk_Cd8-dn.jpeg", width=6, height=5)
ggplot(counts.df %>% filter(window2=="DN/CD8+ > EOMES+CD161+"), aes(x=ZBTB16, y=GZMK)) +
  # geom_point(aes(color=tag_GzmkPlzf))+
  geom_pointdensity() +
  scale_color_viridis_c()+
  theme_cowplot()+
  labs(x="PLZF", title="DN/CD8+ | EOMES+CD161+ (56 cells)") +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
# ggsave("~/Downloads/PlzfGzmk_Cd8-dn-eomespos-cd161pos.jpeg", width=6, height=5)


# CCR7 vs CCR9
table(counts.df$window3)
ggplot(counts.df %>% filter(window3=="CD4+ > EOMES-CD161- > ZBTB16-GZMK-"), aes(x=CCR7, y=CCR9)) +
  # geom_point(aes(color=tag_Ccr7Ccr9))+
  geom_pointdensity() +
  scale_color_viridis_c()+
  theme_cowplot()+
  labs(title="CD4+ | EOMES-CD161- | ZBTB16-GZMK- (686 cells)")
# ggsave("~/Downloads/Ccr7Ccr9_Cd4-Plzfneg-Gzmkneg.jpeg", width=6, height=5)

ggplot(counts.df %>% filter(window3=="CD4+ > EOMES-CD161- > ZBTB16+GZMK-"), aes(x=CCR7, y=CCR9)) +
  # geom_point(aes(color=tag_Ccr7Ccr9))+
  geom_pointdensity() +
  scale_color_viridis_c()+
  theme_cowplot()+
  labs(title="CD4+ | EOMES-CD161- | ZBTB16+GZMK- (120 cells)")
# ggsave("~/Downloads/Ccr7Ccr9_Cd4-Plzfpos-Gzmkneg.jpeg", width=6, height=5)

ggplot(counts.df %>% filter(window3=="DN/CD8+ > EOMES-CD161- > ZBTB16-GZMK-"), aes(x=CCR7, y=CCR9)) +
  # geom_point(aes(color=tag_Ccr7Ccr9))+
  geom_pointdensity() +
  scale_color_viridis_c()+
  theme_cowplot()+
  labs(title="DN/CD8+ | EOMES-CD161- | ZBTB16-GZMK- (938 cells)")
# ggsave("~/Downloads/Ccr7Ccr9_Cd8DN-Plzfneg-Gzmkneg.jpeg", width=6, height=5)

ggplot(counts.df %>% filter(window3=="DN/CD8+ > EOMES-CD161- > ZBTB16+GZMK-"), aes(x=CCR7, y=CCR9)) +
  # geom_point(aes(color=tag_Ccr7Ccr9))+
  geom_pointdensity() +
  scale_color_viridis_c()+
  theme_cowplot()+
  labs(title="DN/CD8+ | EOMES-CD161- | ZBTB16+GZMK- (160 cells)")
# ggsave("~/Downloads/Ccr7Ccr9_Cd8DN-Plzfpos-Gzmkneg.jpeg", width=6, height=5)

# Get it back on UMAP
table(rownames(seur@meta.data[seur@meta.data$group.ident=="NKT_Thymus",]) == rownames(counts.df))
seur@meta.data[seur@meta.data$group.ident=="NKT_Thymus","flow_window1"] <- counts.df$window1
seur@meta.data[seur@meta.data$group.ident=="NKT_Thymus","flow_window2"] <- counts.df$window2
seur@meta.data[seur@meta.data$group.ident=="NKT_Thymus","flow_window3"] <- counts.df$window3
seur@meta.data[seur@meta.data$group.ident=="NKT_Thymus","flow_window4"] <- counts.df$window4

colnames(seur@meta.data)
table(seur@meta.data$flow_window3, useNA="ifany")

SCpubr::do_DimPlot(seur, 
                   group.by = "flow_window1", shuffle=F,
                   label=T,
                   label.color="black",
                   legend.position = "none",
                   repel=T,
                   # colors.use=cols_integrated,
                   font.size = 24)


## 2.3. Flow with CCR7xPLZF ####

# Get all the conditions
counts.df <- counts.df %>%
  mutate("tag_cd161eomes" = ifelse(EOMES < 1 & KLRB1 < 1, "EOMES-CD161-",
                                   ifelse(EOMES < 1 & KLRB1 >= 1, "EOMES-CD161+",
                                          ifelse(EOMES >= 1 & KLRB1 < 1, "EOMES+CD161-",
                                                 ifelse(EOMES >= 1 & KLRB1 >=1, "EOMES+CD161+", "NA"))))) %>%
  mutate("tag_Ccr7Plzf" = ifelse(ZBTB16 < 1 & CCR7 < 1, "ZBTB16-CCR7-",
                                 ifelse(ZBTB16 < 1 & CCR7 >= 1, "ZBTB16-CCR7+",
                                        ifelse(ZBTB16 >= 1 & CCR7 < 1, "ZBTB16+CCR7-",
                                               ifelse(ZBTB16 >= 1 & CCR7 >=1, "ZBTB16+CCR7+", "NA"))))) %>%
  mutate("window1"= tag_cd161eomes,
         "window2"= paste(tag_cd161eomes, tag_Ccr7Plzf, sep=" > "))
# Sanity checks
# table(counts.df$tag_cd161eomes, useNA="ifany")
# table(counts.df$tag_Ccr7Plzf, useNA="ifany")
# table(counts.df$window2, useNA="ifany")

# iNKT cells
ggplot(data.frame(counts.df), aes(x=EOMES, y=KLRB1)) +
  geom_pointdensity() +
  scale_color_viridis_c()+
  theme_cowplot()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

# Window 1
ggplot(counts.df %>% filter(window1 == "EOMES-CD161-"), aes(x=CCR7, y=ZBTB16)) +
  # geom_point(aes(color=window2))+
  geom_pointdensity() +
  scale_color_viridis_c()+
  theme_cowplot()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

# Window 2
ggplot(counts.df %>% filter(window2 == "EOMES-CD161- > ZBTB16+CCR7-"), aes(x=CD8A, y=CD4)) +
  geom_pointdensity() +
  scale_color_viridis_c()+
  theme_cowplot()+
  labs(title="EOMES-CD161- > ZBTB16+CCR7-")+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
ggplot(counts.df %>% filter(window2 == "EOMES-CD161- > ZBTB16+CCR7-")) +
  geom_histogram(aes(x=SELL), bins = 50)+
  theme_cowplot()+
  labs(title="EOMES-CD161- > ZBTB16+CCR7-")

ggplot(counts.df %>% filter(window2 == "EOMES-CD161- > ZBTB16-CCR7-"), aes(x=CD8A, y=CD4)) +
  geom_pointdensity() +
  scale_color_viridis_c()+
  theme_cowplot()+
  labs(title="EOMES-CD161- > ZBTB16-CCR7-")+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
ggplot(counts.df %>% filter(window2 == "EOMES-CD161- > ZBTB16-CCR7-")) +
  geom_histogram(aes(x=SELL), bins = 50)+
  theme_cowplot()+
  labs(title="EOMES-CD161- > ZBTB16-CCR7-")

ggplot(counts.df %>% filter(window2 == "EOMES-CD161- > ZBTB16-CCR7+"), aes(x=CD8A, y=CD4)) +
  geom_pointdensity() +
  scale_color_viridis_c()+
  theme_cowplot()+
  labs(title="EOMES-CD161- > ZBTB16-CCR7+")+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
ggplot(counts.df %>% filter(window2 == "EOMES-CD161- > ZBTB16-CCR7+")) +
  geom_histogram(aes(x=SELL), bins = 50)+
  theme_cowplot()+
  labs(title="EOMES-CD161- > ZBTB16-CCR7+")


# Add the windows to seurat iNKT object
table(rownames(seur.thynkt@meta.data) == rownames(counts.df))
seur.thynkt@meta.data[,"flow_window1"] <- counts.df$window1
seur.thynkt@meta.data[,"flow_window2"] <- counts.df$window2

# Plot the coexpression
library(SCpubr)
# do_DimPlot(seur.thynkt, group.by="new_clusters_NKT")
FeaturePlot(seur.thynkt, features=c("EOMES", "KLRB1"), blend=T,
            cols=c("#e0e0e0", "#2166ac", "#b2182b"), order=T, blend.threshold = 0.1)
# ggsave(file.path(path.plots,"thyNKT_eomes-klrb1.jpeg"), width=18, height=5)
# FeaturePlot(seur.thynkt,
#             features=c("CCR7", "ZBTB16"), blend=T,
#             cols=c("#e0e0e0", "#2166ac", "#b2182b"), order=T, blend.threshold = 0.1)
# FeaturePlot(seur.thynkt,
#             features=c("CD8A", "CD4"), blend=T,
#             cols=c("#e0e0e0", "#2166ac", "#b2182b"), order=T, blend.threshold = 0.1)

# Only EOMES- CD161- cells
FeaturePlot(seur.thynkt,
            cells=rownames(seur.thynkt@meta.data[seur.thynkt@meta.data$flow_window1=="EOMES-CD161-",]),
            features=c("CCR7", "ZBTB16"), blend=T,
            cols=c("#e0e0e0", "#2166ac", "#b2182b"), order=T, blend.threshold = 0.1)+
  ggtitle("EOMES-CD161- cells")
# ggsave(file.path(path.plots,"thyNKT_eomesNeg-klrb1Neg_ccr7-plzf.jpeg"), width=18, height=5)
FeaturePlot(seur.thynkt,
            cells=rownames(seur.thynkt@meta.data[seur.thynkt@meta.data$flow_window1=="EOMES-CD161-",]),
            features=c("CCR9", "ZBTB16"), blend=T,
            cols=c("#e0e0e0", "#2166ac", "#b2182b"), order=T, blend.threshold = 0.1)+
  ggtitle("EOMES-CD161- cells")
# ggsave(file.path(path.plots,"thyNKT_eomesNeg-klrb1Neg_ccr9-plzf.jpeg"), width=18, height=5)

# Only CCR7- PLZF+ cells
FeaturePlot(seur.thynkt,
            cells=rownames(seur.thynkt@meta.data[seur.thynkt@meta.data$flow_window2=="EOMES-CD161- > ZBTB16+CCR7-",]),
            features=c("CD8A", "CD4"), blend=T,
            cols=c("#e0e0e0", "#2166ac", "#b2182b"), order=F, blend.threshold = 0.1)+
  ggtitle("ZBTB16+CCR7- cells")
# ggsave(file.path(path.plots,"thyNKT_eomesNeg-klrb1Neg_ccr7Neg-plzfPos_cd8a-cd4.jpeg"), width=18, height=5)
# Only CCR7- PLZF- cells
FeaturePlot(seur.thynkt,
            cells=rownames(seur.thynkt@meta.data[seur.thynkt@meta.data$flow_window2=="EOMES-CD161- > ZBTB16-CCR7-",]),
            features=c("CD8A", "CD4"), blend=T,
            cols=c("#e0e0e0", "#2166ac", "#b2182b"), order=F, blend.threshold = 0.1)+
  ggtitle("ZBTB16-CCR7- cells")
# ggsave(file.path(path.plots,"thyNKT_eomesNeg-klrb1Neg_ccr7Neg-plzfNeg_cd8a-cd4.jpeg"), width=18, height=5)
# Only CCR7+ PLZF- cells
FeaturePlot(seur.thynkt,
            cells=rownames(seur.thynkt@meta.data[seur.thynkt@meta.data$flow_window2=="EOMES-CD161- > ZBTB16-CCR7+",]),
            features=c("CD8A", "CD4"), blend=T,
            cols=c("#e0e0e0", "#2166ac", "#b2182b"), order=F, blend.threshold = 0.1)+
  ggtitle("ZBTB16-CCR7+ cells")
# ggsave(file.path(path.plots,"thyNKT_eomesNeg-klrb1Neg_ccr7Pos-plzfNeg_cd8a-cd4.jpeg"), width=18, height=5)

# SELL expression
Idents(seur.thynkt) <- "flow_window1"
RidgePlot(seur.thynkt,
          idents="EOMES-CD161-",
          features = "SELL", group.by="flow_window2", same.y.lims=T, y.max=100)+
  theme(legend.position="none")
VlnPlot(seur.thynkt,
          idents="EOMES-CD161-",
          features = "SELL", group.by="flow_window2")+
  xlab("")+
  theme(legend.position="none")
# ggsave(file.path(path.plots,"thyNKT_eomesNeg-klrb1Neg_Sell-vlnplot.jpeg"), width=8, height=6)

# PLZF expression
VlnPlot(seur.thynkt,
        # idents="EOMES-CD161-",
        features = "ZBTB16", group.by="flow_window1")+
  xlab("")+
  theme(legend.position="none")

counts.df %>%
  as_tibble() %>%
  filter(window1=="EOMES-CD161-") %>%
  filter(ZBTB16>0)
