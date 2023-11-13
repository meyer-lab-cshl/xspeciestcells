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
library(Seurat)

# Import human data
seur.nkt  <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.nkt.RDS")
seur.mait <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.mait.RDS")
seur.gdt  <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.gd.RDS")
seur.ms <- readRDS("./data/cross-species/08_innateT_cross_species/Reanalysis_wGD_filtered_seurat_MNN.rds")

# Plot DimPlot
DimPlot(seur.nkt,  group.by="cell_annot", label=T)
DimPlot(seur.mait, group.by="cell_annot", label=T)
DimPlot(seur.gdt,  group.by="cell_annot", label=T)
DimPlot(seur.ms,  group.by="Mouse_clusters", label=T)

# Cleanup a bit
seur.nkt@meta.data[,33:75]  <- NULL
seur.mait@meta.data[,33:75] <- NULL
seur.gdt@meta.data[,33:75]  <- NULL
seur.ms@meta.data[,22:301]  <- NULL

# ORTHOLOGS table
ortholog.df <- read.csv("./data/cross-species/03_BiomartTable/big_ass_ortholog_table.csv")
length(unique(ortholog.df$ms_symbol_data)) # 17,085 ms genes
length(unique(ortholog.df$hu_symbol)) # 17,152 hu genes




# ****************
# 3. ANALYSIS ####
# ****************

#___________________________
## 3.1. Cleanup metadata ####
colnames(seur.ms@meta.data)[21] <- "cell_annot"
# DimPlot(seur.ms,  group.by="cell_annot", label=T) # sanity check

## /end ####


#___________________________
## 3.2. Prepare ortholog table ####

# First, check whether genes can all be found in the ortholog table (many can't be found because I removed genes with orthology confidence=0)
table(unique(rownames(seur.ms)) %in% unique(ortholog.df$ms_symbol_bmt)) # 9,039 not
table(unique(rownames(seur.nkt)) %in% ortholog.df$hu_symbol) # 4,757 not

# Subset the ortholog table to only genes that we can "translate"
dictionary <- ortholog.df %>%
  as_tibble() %>%
  select(ms_ensemblID, version_msgene, ms_symbol_bmt, hu_ensemblID, hu_symbol, hu_homology_type, hu_orthology_confidence) %>%
  ## Intersection
  filter(ms_symbol_bmt %in% unique(rownames(seur.ms)) & hu_symbol %in% unique(rownames(seur.nkt))) %>%
  ## Remove any symbols that are NAs
  filter(!is.na(ms_symbol_bmt)) %>%
  filter(!is.na(hu_symbol)) %>%
  ## Keep only 1:1 orthologs
  group_by(ms_symbol_bmt) %>% filter(n_distinct(hu_symbol) == 1) %>% ungroup() %>%
  group_by(hu_symbol) %>% filter(n_distinct(ms_symbol_bmt) == 1) %>% ungroup() %>%
  ## remove duplicate rows
  distinct()
dim(dictionary) # 11,074 genes

# Subset seurat objects to genes that we have in the dictionnary (can't do comparison on genes we can't translate between species)
seur.nkt  <- seur.nkt[unique(dictionary$hu_symbol),]
seur.mait <- seur.mait[unique(dictionary$hu_symbol),]
seur.gdt  <- seur.gdt[unique(dictionary$hu_symbol),]
seur.ms   <- seur.ms[unique(dictionary$ms_symbol_bmt),] # keep only genes that we'll be able to translate


## /end ####


#___________________________
## 3.2. Get HVGs ####

# Human HVGs
hu.hvg <- unique(c(VariableFeatures(seur.nkt), VariableFeatures(seur.mait), VariableFeatures(seur.gdt)))
length(hu.hvg) # 2578 genes

# Mouse HVGs
ms.hvg <- VariableFeatures(FindVariableFeatures(seur.ms, nfeatures=length(hu.hvg)))
length(ms.hvg) # 2578 genes

# Translate the mouse HVGs into "human gene" language
ms.hvg.translated <- unique(pull(dictionary %>% filter(ms_symbol_bmt %in% ms.hvg), # not all ms HVGs are found in ortholog.df
                          hu_symbol))
hu.hvg.translated <- unique(pull(dictionary %>% filter(hu_symbol %in% hu.hvg), # not all hu HVGs are found in ortholog.df
                          hu_symbol))
total.hvg <- unique(union(ms.hvg.translated, hu.hvg.translated))
length(total.hvg) # 4229 genes


# VennDiagram::venn.diagram(
#   x = list(VariableFeatures(seur.nkt), VariableFeatures(seur.mait), VariableFeatures(seur.gdt), ms.hvg.translated),
#   category.names = c("hvg_humanNKT" , "hvg_humanMAIT " , "hvg_humanGD", "hvg_mouseall"),
#   filename = './scripts-in-progress/cross-species/08_innateT_cross_species/plots/hvg_venn.png',
#   output=TRUE
# )

## /end ####


#___________________________
## 3.3. Get counts ####

# counts
hu.nkt.counts  <- seur.nkt[["RNA"]]@counts
hu.mait.counts <- seur.mait[["RNA"]]@counts
hu.gdt.counts  <- seur.gdt[["RNA"]]@counts
hu.counts <-  cbind(hu.nkt.counts, hu.mait.counts, hu.gdt.counts)
ncol(hu.counts) == ncol(hu.nkt.counts)+ncol(hu.mait.counts)+ncol(hu.gdt.counts)
dim(hu.counts)
ms.counts      <- seur.ms[["RNA"]]@counts
dim(ms.counts)


# Keep only ms and hu genes that have 1:1 orthologs
# table(unique(rownames(ms.counts)) %in% dictionary$ms_symbol_bmt) # sanity check all TRUE
# table(unique(rownames(hu.counts)) %in% dictionary$hu_symbol) # sanity check all TRUE
ms.counts <- ms.counts[rownames(ms.counts) %in% dictionary$ms_symbol_bmt,]
hu.counts <- hu.counts[rownames(hu.counts) %in% dictionary$hu_symbol,]
nrow(ms.counts)==nrow(hu.counts) # should be TRUE


# Translate the mouse genes in count table into "human gene"
ms.dict <- dictionary %>%
  filter(ms_symbol_bmt %in% rownames(ms.counts)) %>%
  select(ms_symbol_bmt, hu_symbol, hu_orthology_confidence) %>%
  # distinct() %>%
  # group_by(ms_symbol_bmt) %>% filter(n_distinct(hu_symbol)>1)
  distinct(ms_symbol_bmt, .keep_all=T)
ms.dict <- ms.dict[match(rownames(ms.counts), ms.dict$ms_symbol_bmt),]
table(ms.dict$ms_symbol_bmt == rownames(ms.counts)) # should be all TRUE
table(is.na(ms.dict$hu_symbol)) # should have no NAs
# Translate
rownames(ms.counts) <- ms.dict$hu_symbol


# Verify mouse genes (in "human symbols") correspond to human genes
table(rownames(ms.counts) %in% rownames(hu.counts))

## /end ####


#___________________________
## 3.4. Get metadata ####

# get metadata
hu.metadata <- rbind(seur.nkt@meta.data[,c("cell.ident", "cell_annot")],
                     seur.mait@meta.data[,c("cell.ident", "cell_annot")],
                     seur.gdt@meta.data[,c("cell.ident", "cell_annot")])
# table(hu.metadata$cell.ident, useNA="ifany")
ms.metadata <- seur.ms@meta.data[,c("Cell.type", "cell_annot")]
colnames(ms.metadata)[1] <- "cell.ident"

# add study
ms.metadata$study <- "Mouse"
hu.metadata$study <- "Human"

head(ms.metadata)
head(hu.metadata)

## /end ####


#___________________________
## 3.5. Merge everything into one seurat object ####

# sanity checks
# table(rownames(hu.metadata)==colnames(hu.counts), useNA="ifany") # 10k cells
# table(rownames(ms.metadata)==colnames(ms.counts), useNA="ifany") # 45k cells

ms.seur <- CreateSeuratObject(counts=ms.counts, meta.data=ms.metadata)
hu.seur <- CreateSeuratObject(counts=hu.counts, meta.data=hu.metadata)
seur.total <- merge(ms.seur, hu.seur)

# convert to .h5ad
SeuratDisk::SaveH5Seurat(seur.total,
                         filename = "data/cross-species/08_innateT_cross_species/cluster/input/innateT_ms_hu_full_seu_gene_names.h5Seurat")
SeuratDisk::Convert("data/cross-species/08_innateT_cross_species/cluster/input/innateT_ms_hu_full_seu_gene_names.h5Seurat", dest = "h5ad")

# export HVG list into .csv file
hvg.df <- data.frame("features"=rownames(seur.total),
                     "highly_variable"=ifelse(rownames(seur.total)%in%total.hvg, T, F))
table(hvg.df$highly_variable, useNA="ifany")
table(hvg.df[hvg.df$highly_variable==T, "features"]%in%total.hvg, useNA="ifany") # sanity check
write.csv(hvg.df, "data/cross-species/08_innateT_cross_species/cluster/input/list_hvg.csv")

## /end ####



#___________________________
## 3.6. Bubbleplot (post-Elzar) ####

library(RColorBrewer)
library(reshape2)
library(ggrepel)
library(patchwork)

# import matrix
mtn <- read.csv("./data/cross-species/08_innateT_cross_species/cluster/output/pymtn_crosspecies_innateT_slowversion_2023-11-13.csv", row.names=1)
mtn <- read.csv("./data/cross-species/08_innateT_cross_species/cluster/output/pymtn_crosspecies_innateT_slowversion_mousebylineage_2023-11-13.csv", row.names=1)
mtn <- read.csv("./data/cross-species/08_innateT_cross_species/cluster/output/pymtn_crosspecies_innateT_slowversion_allgenes_2023-11-13.csv", row.names=1)
mtn <- as.matrix(mtn)

# heatmap
gplots::heatmap.2(mtn,
          # trace
          trace="none",
          # dendrogram
          # Rowv=FALSE,
          # Colv=FALSE,
          # dendrogram="none",
          # superimpose a density histogram on color key
          density.info="none",
          # color scale
          col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100)),
          # breaks=seq(0,1,length=101),
          key.xlab = "AUROC",
          key.title="",
          keysize = 1.2,
          # text labels
          main="Mouse vs Human innate T (4229 HVGs)",
          cexRow=0.6,
          cexCol=0.6,
          # margins
          margins=c(9,9))

# melt
mtn.df <- melt(mtn)

# Define plotting function
bubble_plot <- function(df, order_x, order_y, label_x, label_y, auroc_min){
  
  print(head(df))
  
  # x barplot
  bp.x <- ggplot(data=df %>% select(var_x, ncells_x) %>% distinct(),
                 aes(x=factor(var_x, levels=order_x), y=ncells_x))+
    geom_bar(stat="identity", fill="#bdbdbd") + theme_cowplot()+
    scale_x_discrete(position="top")+
    labs(y="#cells")+
    theme(axis.text = element_text(size=15),
          axis.text.x=element_text(angle=45, hjust=0),
          axis.ticks.y = element_blank(),
          axis.title.y=element_text(size=20),
          axis.title.x = element_blank(),
          axis.line.x=element_blank(),
          legend.position = "none",
          plot.margin = margin(10,30,10,10))
  
  # PROPORTION OF MOUSE GDT CELLS IN EACH CLUSTER
  bp.y <- ggplot(data=df%>% select(var_y, ncells_y) %>% distinct(),
                 aes(x=factor(var_y, levels=order_y), y=ncells_y))+
    geom_bar(stat="identity", fill="#bdbdbd") +
    scale_x_discrete(position="top") +
    labs(y="#cells")+ coord_flip() + theme_cowplot()+
    theme(axis.title.y = element_blank(),
          axis.text = element_text(size=15),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(size=20),
          axis.line.y=element_blank(),
          legend.position = "none")
  
  # BUBBLE PLOT
  hm.clean <- ggplot(df, aes(x=factor(var_x, levels=order_x),
                             y=factor(var_y, levels=order_y))) +
    geom_point(aes(size = auroc, color= auroc))+
    geom_text(data=df %>% filter(auroc>auroc_min) %>% mutate(across("auroc", \(x) round(x,2))), aes(label=auroc), color="white")+
    scale_size_continuous(limits=c(0,1), breaks=seq(0,1, by=0.2), range = c(1, 15))+
    scale_color_gradient2(low="#2166ac", mid="white", high="#a50f15", midpoint=0.5, limits=c(0,1), name="AUROC", breaks=seq(0,1, by=0.2))+
    labs(x=label_x,y=label_y, size="AUROC")+
    theme_cowplot()+
    theme(legend.position="bottom", legend.key.width = unit(0.8, 'cm'),
          axis.text = element_text(size=15), axis.text.x = element_text(angle=45, hjust=1), axis.title=element_text(size=20))
  
  # COMBINE
  p <- (bp.x+plot_spacer() + plot_layout(widths = c(5, 1))) / (hm.clean + bp.y + plot_layout(widths = c(5, 1))) + plot_layout(heights = c(1, 5))
  return(p)
}



## 3.6.1. Human x Mouse bubble plot ####
ms.seur$cell_annot_lineage <- paste0(ms.seur$cell.ident, "_", ms.seur$cell_annot)
mtn.df_mshu <- mtn.df %>%
  filter(str_detect(Var1,"Human")) %>%
  mutate(Var1 = gsub("Human\\|", "", Var1)) %>%
  filter(str_detect(Var2, "Mouse")) %>%
  mutate(Var2 = gsub("Mouse\\.", "", Var2)) %>%
  left_join(as.data.frame(table(hu.seur$cell_annot)), by="Var1") %>%
  rename(var_x=Var1, Var1=Var2, auroc=value, ncells_x=Freq) %>%
  left_join(as.data.frame(table(ms.seur$cell_annot_lineage)), by="Var1") %>%
  rename(var_y=Var1, ncells_y=Freq) %>%
  # keep min 100 cells
  filter(ncells_x>=100 & ncells_y>=100)

order_mouse <- rev(c("immature_GzmA_Gd", "immature_DP", "stage_0_signaling", "cycling_I", "cycling_II", "type_I", "type_II", "type_III"))
order_mouse_2 <- rev(paste0(c("GD_", "MAIT_", "NKT_"),
                            c(rep("immature_GzmA_Gd", 3), rep("immature_DP", 3), rep("stage_0_signaling", 3), rep("cycling_I", 3), rep("cycling_II", 3), rep("type_I", 3), rep("type_II", 3), rep("type_III", 3))))
bubble_plot(df=mtn.df_mshu,
            order_x = unique(mtn.df_mshu$var_x),
            order_y = order_mouse_2,
            label_x = "Human clusters",
            label_y = "Mouse clusters",
            auroc_min= 0.6)
ggsave("./scripts-in-progress/cross-species/08_innateT_cross_species/plots/innateT_allgenes11072_mshu.pdf", width=17, height=14)


## 3.6.2. Human x Human bubble plot ####
mtn.df_huhu <- mtn.df %>%
  filter(str_detect(Var1,"Human")) %>%
  mutate(Var1 = gsub("Human\\|", "", Var1)) %>%
  filter(str_detect(Var2, "Human")) %>%
  mutate(Var2 = gsub("Human\\.", "", Var2)) %>%
  left_join(as.data.frame(table(hu.seur$cell_annot)), by="Var1") %>%
  rename(var_x=Var1, Var1=Var2, auroc=value, ncells_x=Freq) %>%
  left_join(as.data.frame(table(hu.seur$cell_annot)), by="Var1") %>%
  rename(var_y=Var1, ncells_y=Freq) %>%
  # keep min 100 cells
  filter(ncells_x>=100 & ncells_y>=100)

bubble_plot(df=mtn.df_huhu,
            order_x = unique(mtn.df_huhu$var_x),
            order_y = rev(unique(mtn.df_huhu$var_x)),
            label_x = "Human clusters",
            label_y = "Human clusters",
            auroc_min= 0.6)
ggsave("./scripts-in-progress/cross-species/08_innateT_cross_species/plots/innateT_allgenes11072_huhu.pdf", width=15, height=15)


## 3.6.3. Mouse x Mouse bubble plot ####
mtn.df_msms <- mtn.df %>%
  filter(str_detect(Var1,"Mouse")) %>%
  mutate(Var1 = gsub("Mouse\\|", "", Var1)) %>%
  filter(str_detect(Var2, "Mouse")) %>%
  mutate(Var2 = gsub("Mouse\\.", "", Var2)) %>%
  left_join(as.data.frame(table(ms.seur$cell_annot_lineage)), by="Var1") %>%
  rename(var_x=Var1, Var1=Var2, auroc=value, ncells_x=Freq) %>%
  left_join(as.data.frame(table(ms.seur$cell_annot_lineage)), by="Var1") %>%
  rename(var_y=Var1, ncells_y=Freq) %>%
  # keep min 100 cells
  filter(ncells_x>=100 & ncells_y>=100)

bubble_plot(df=mtn.df_msms,
            order_x = rev(order_mouse_2),
            order_y = order_mouse_2,
            label_x = "Mouse clusters",
            label_y = "Mouse clusters",
            auroc_min= 0.8)
ggsave("./scripts-in-progress/cross-species/08_innateT_cross_species/plots/innateT_allgenes11072_msms.pdf", width=18, height=18)
