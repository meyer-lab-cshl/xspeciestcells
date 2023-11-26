###
# Purpose: Metaneighbor on ms vs hu MAIT
# Date: May 4th 2024
# Author: Salomé Carcy
###




# Import librairies
library(Seurat)
library(MetaNeighbor)
library(SummarizedExperiment)
library(gplots)
library(dendextend)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(tidyverse)
setwd("~/Projects/HumanThymusProject/")

# Import data
seur.ms <- readRDS("./data/cross-species/00_Reproduce_UMAPs/ms_mait_seurobj.rds")
# DimPlot(seur.ms)
cols_mait <- c("MAIT0"     = "#fee391",
               "Cluster 7" = "#bf812d",
               "MAIT1"     = "#80cdc1",
               "MAIT17a"   = "#084594",
               "MAIT17b"   = "#6baed6",
               "CyclingS"  = "#f768a1",
               "CyclingG2M"= "#8c6bb1")
SCpubr::do_DimPlot(seur.ms, 
                   group.by = "cell_type",
                   label=T,
                   label.color="black",
                   legend.position = "none",
                   repel=T,
                   plot.title = "Mouse",
                   colors.use=cols_mait,
                   font.size = 24)
# ggsave("./data/cross-species/04_Metaneighbor_mait/ms_mait_umap.jpeg", width=6, height=6)


seur.hu <- readRDS("./data/human-thymus/HumanData_12_AnalysisByLineage/thymus.MAIT_03_02_23.RDS")
colors_clusters_MAIT <- c("MAIT_c0" = "#d8443c", "MAIT_c1" = "#e09351", "MAIT_c2" = "gold", "MAIT_c3" = "#74c8c3", "MAIT_c5" = "#5a97c1",
                         "MAIT_c4" = "#a40000", "MAIT_c6" = "#fec44f")
SCpubr::do_DimPlot(seur.hu, 
                   group.by = "new_clusters_MAIT",
                   label=T,
                   label.color="black",
                   legend.position = "none",
                   repel=T,
                   plot.title = "Human",
                   colors.use=colors_clusters_MAIT,
                   font.size = 24)
# ggsave("./data/cross-species/04_Metaneighbor_mait/hu_mait_umap.jpeg", width=6, height=6)


ortholog.df <- read.csv("~/Projects/HumanThymusProject/data/cross-species/03_BiomartTable/big_ass_ortholog_table.csv")
length(unique(ortholog.df$ms_symbol_data)) # 17,085 ms genes
length(unique(ortholog.df$hu_symbol)) # 17,152 hu genes



# Get counts, HVGs and metadata
ms.hvg <- VariableFeatures(FindVariableFeatures(seur.ms))
ms.counts <- seur.ms[["RNA"]]@counts
ms.metadata <- seur.ms@meta.data

hu.hvg <- VariableFeatures(seur.hu)
hu.counts <- seur.hu[["RNA"]]@counts
hu.metadata <- seur.hu@meta.data


# First, check whether genes can all be found in the ortholog table (many can't be found because I removed genes with orthology confidence=0)
table(unique(rownames(ms.counts)) %in% ortholog.df$ms_symbol_data) # 1,037 not
table(unique(rownames(hu.counts)) %in% ortholog.df$hu_symbol) # 4747 not


# Subset the ortholog table to only genes that we can "translate"
dictionary <- ortholog.df %>%
  as_tibble() %>%
  ## Intersection
  filter(ms_symbol_data %in% unique(rownames(ms.counts)) & hu_symbol %in% unique(rownames(hu.counts))) %>%
  ## Union
  # filter(ms_symbol_data %in% unique(rownames(ms.counts)) | hu_symbol %in% unique(rownames(hu.counts))) %>%
  # filter(hu_orthology_confidence == 1) %>%
  ## Remove any symbols that are NAs
  filter(!is.na(ms_symbol_data)) %>%
  filter(!is.na(hu_symbol)) %>%
  ## Keep only 1:1 orthologs
  group_by(ms_symbol_data) %>% filter(n_distinct(hu_symbol) == 1) %>% ungroup() %>%
  group_by(hu_symbol) %>% filter(n_distinct(ms_symbol_data) == 1) %>% ungroup()


# Translate the mouse HVGs into "human gene" language
ms.hvg.translated <- pull(dictionary %>% filter(ms_symbol_data %in% ms.hvg), # not all ms HVGs are found in ortholog.df
                          hu_symbol) # 1604 ms HVG
hu.hvg.translated <- pull(dictionary %>% filter(hu_symbol %in% hu.hvg), # not all hu HVGs are found in ortholog.df
                          hu_symbol) # 1022 hu HVG
total.hvg <- unique(union(ms.hvg.translated, hu.hvg.translated))
length(total.hvg) # 2,240 genes


# Keep only ms and hu genes that have 1:1 orthologs
table(unique(rownames(ms.counts)) %in% dictionary$ms_symbol_data) # 9,788 genes should have a translation
table(unique(rownames(hu.counts)) %in% dictionary$hu_symbol) # 9,788 genes should have a translation
ms.counts <- ms.counts[rownames(ms.counts) %in% dictionary$ms_symbol_data,]
hu.counts <- hu.counts[rownames(hu.counts) %in% dictionary$hu_symbol,]


# Translate the mouse genes in count table into "human gene"
ms.dict <- dictionary %>%
  filter(ms_symbol_data %in% rownames(ms.counts)) %>%
  select(ms_symbol_data, hu_symbol, hu_orthology_confidence) %>%
  # distinct() %>%
  # group_by(ms_symbol_data) %>% filter(n_distinct(hu_symbol)>1)
  distinct(ms_symbol_data, .keep_all=T)
ms.dict <- ms.dict[match(rownames(ms.counts), ms.dict$ms_symbol_data),]
table(ms.dict$ms_symbol_data == rownames(ms.counts)) # all true
table(is.na(ms.dict$hu_symbol)) # no NAs
# Translate
rownames(ms.counts) <- ms.dict$hu_symbol
# table(rownames(ms.counts)%in%rownames(hu.counts)) # sanity check


# Merge everything into one
ms.metadata$study <- "Mouse"
ms.metadata <- ms.metadata[,c("cell_type", "study")]
colnames(ms.metadata)[1] <- "clusters_MAIT"
head(ms.metadata)

hu.metadata$study <- "Human"
hu.metadata <- hu.metadata[,c("new_clusters_MAIT", "study")]
colnames(hu.metadata)[1] <- "clusters_MAIT"
head(hu.metadata)

ms.counts <- ms.counts[rownames(hu.counts),]
table(rownames(ms.counts)==rownames(hu.counts), useNA="ifany")
table(colnames(ms.metadata)==colnames(hu.metadata), useNA="ifany")
tot.counts   <- cbind(ms.counts, hu.counts)
tot.metadata <- rbind(ms.metadata, hu.metadata)
table(colnames(tot.counts)==rownames(tot.metadata), useNA="ifany")

# Convert seurat count matrix to SummarizedExperiment object
count <- SummarizedExperiment(assays=tot.counts,
                              colData=tot.metadata)

# ms.seur <- CreateSeuratObject(counts=ms.counts, meta.data=ms.metadata)
# hu.seur <- CreateSeuratObject(counts=hu.counts, meta.data=hu.metadata)
# seur.total <- merge(ms.seur, hu.seur) # 9788 genes
# 
# # Sanity checks
# head(seur.total@meta.data)
# table(seur.total$clusters_MAIT, useNA="ifany")
# table(seur.total$study, useNA="ifany")


# Convert seurat count matrix to SummarizedExperiment object
# count <- SummarizedExperiment(assays=seur.total@assays[["RNA"]]@counts,
#                               colData=seur.total$clusters_MAIT)




# _______________________
# METANEIGHBOR
# Run metaneighbor
mtn <- MetaNeighborUS(var_genes=total.hvg,
                      dat=count,
                      study_id=count$study,
                      cell_type=count$clusters_MAIT,
                      fast_version=T)
# saveRDS(mtn, "./data/cross-species/04_Metaneighbor_mait/mait_ms-hu_mtnslowversion.rds")
mtn <- readRDS("./data/cross-species/04_Metaneighbor_mait/mait_ms-hu_mtnslowversion.rds")

# Heatmap
mtn.sub <- mtn[1:7,8:14]
# mtn.sub <- mtn.sub[rev(rownames(mtn.sub)),]

jpeg("./data/cross-species/04_Metaneighbor_mait/mait_ms-hu_fastversion_fulltree.jpeg", width=1500, height=1500, res=200)
heatmap.2(mtn,
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
          breaks=seq(0,1,length=101),
          key.xlab = "AUROC",
          key.title="",
          keysize = 1.2,
          # text labels
          main="Mouse vs Human MAIT",
          cexRow=0.6,
          cexCol=0.6,
          # margins
          margins=c(6,6)
          )
dev.off()


# Bubble plot
library(reshape2)
library(ggrepel)
mtn.df <- melt(mtn)
mtn.df <- mtn.df %>%
  filter(str_detect(Var1,"Human")) %>%
  mutate(Var1 = gsub("Human\\|", "", Var1)) %>%
  filter(str_detect(Var2, "Mouse")) %>%
  mutate(Var2 = gsub("Mouse\\|", "", Var2))


# Get number of cells per human cluster and per mouse cluster
mtn.df <- mtn.df %>%
  as_tibble() %>%
  # add nb of cells per human cluster
  left_join(as.data.frame(table(seur.hu$new_clusters_MAIT)), by="Var1") %>%
  dplyr::rename(ncells_human=Freq, human=Var1, Var1=Var2) %>%
  mutate(totalcells_human = dim(seur.hu)[2],
         propcells_human = ncells_human*100/totalcells_human) %>%
  # add nb of cells per mouse cluster
  left_join(as.data.frame(table(seur.ms$cell_type)), by="Var1") %>%
  dplyr::rename(ncells_mouse=Freq, mouse=Var1, auroc=value) %>%
  mutate(totalcells_mouse = dim(seur.ms)[2],
         propcells_mouse = ncells_mouse*100/totalcells_mouse)
# saveRDS(mtn.df, "./data/cross-species/04_Metaneighbor_mait/mait_ms-hu_mtnslowversion_DF.rds")


# PROPORTION OF HUMAN MAIT CELLS IN EACH CLUSTER
bp.x <- ggplot(data=mtn.df %>% select(human,propcells_human) %>% distinct(),
               aes(x=factor(human, levels=paste0("MAIT_c", 0:6)), y=propcells_human))+
  geom_bar(stat="identity", fill="#bdbdbd") + theme_cowplot()+
  scale_x_discrete(position="top")+
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100))+
  labs(y="%cells")+
  theme(axis.text.y = element_text(size=15),
        axis.ticks.y = element_blank(),
        axis.title.y=element_text(size=20),
        axis.text.x = element_text(size = 15, angle=45, hjust=0),
        axis.title.x = element_blank(), axis.line.x=element_blank(),
        legend.position = "none")


# PROPORTION OF MOUSE MAIT CELLS IN EACH CLUSTER
order_mouse <- c("MAIT0", "Cluster 7", "CyclingG2M", "CyclingS", "MAIT1", "MAIT17a", "MAIT17b")
bp.y <- ggplot(data=mtn.df %>% select(mouse,propcells_mouse) %>% distinct(),
               aes(x=factor(mouse, levels=rev(order_mouse)), y=propcells_mouse))+
  geom_bar(stat="identity", fill="#bdbdbd") +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100))+
  scale_x_discrete(position="top") +
  labs(y="%cells")+ coord_flip() + theme_cowplot()+
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size=15),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size=20),
        axis.line.y=element_blank(),
        legend.position = "none")



# BUBBLE PLOT
# library(scales)
hm.clean <- ggplot(mtn.df, aes(x=factor(human, levels=paste0("MAIT_c", 0:6)),
                               y=factor(mouse, levels=rev(order_mouse)))) +
  geom_point(aes(size=auroc, color= auroc))+
  geom_text(data=mtn.df %>% filter(auroc>0.6) %>% mutate(across("auroc", \(x) round(x,2))), aes(label=auroc), color="white")+
  scale_size_continuous(limits=c(0,1), breaks=seq(0.2,0.8, by=0.2), range = c(1, 15))+
  scale_colour_gradient2(low="#d9d9d9", mid="white", high="#a50f15", midpoint=0.5, limits=c(0,1), name="AUROC", breaks=seq(0,1, by=0.2))+
  labs(x="Human clusters",y="Mouse clusters", size="AUROC")+
  theme_cowplot()+
  theme(legend.position="bottom", legend.key.width = unit(0.8, 'cm'),
        axis.text = element_text(size=15), axis.title=element_text(size=20), axis.text.x=element_text(angle=45, hjust=1))


# COMBINE
library(patchwork)
(bp.x+plot_spacer() + plot_layout(widths = c(5, 1))) / (hm.clean + bp.y + plot_layout(widths = c(5, 1))) + plot_layout(heights = c(1, 5))
# ggsave("./data/cross-species/04_Metaneighbor_mait/mait_ms-hu_metaneighbor_fastversion_bubbleplot.jpeg", width=9, height=8)





# _______________________
# METANEIGHBOR TROUBLESHOOTING FAST/SLOW
# get slow mtn
mtn_slow <- readRDS("./data/cross-species/04_Metaneighbor_mait/mait_ms-hu_mtnslowversion.rds")

# plot fast mtn
htm_fast <- heatmap.2(mtn,
          # trace
          trace="none",
          # superimpose a density histogram on color key
          density.info="none",
          # color scale
          col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100)),
          breaks=seq(0,1,length=101),
          key.xlab = "AUROC",
          key.title="",
          keysize = 1.2,
          # text labels
          main="Mouse vs Human MAIT",
          cexRow=0.6,
          cexCol=0.6,
          # margins
          margins=c(6,6)
)

cluster_order <- rownames(mtn)[htm_fast$rowInd]

# plot slow mtn in same order
mtn_slow <- mtn_slow[rev(cluster_order), cluster_order]
jpeg("./data/cross-species/04_Metaneighbor_mait/mait_ms-hu_slowversion_fulltree_orderedlikefastversion.jpeg", width=1500, height=1500, res=200)
heatmap.2(mtn_slow,
          # trace
          trace="none",
          # dendrogram
          dendrogram="none",
          Rowv=F,
          Colv=F,
          # superimpose a density histogram on color key
          density.info="none",
          # color scale
          col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100)),
          breaks=seq(0,1,length=101),
          key.xlab = "AUROC",
          key.title="",
          keysize = 1.2,
          # text labels
          main="Mouse vs Human MAIT",
          cexRow=0.6,
          cexCol=0.6,
          # margins
          margins=c(6,6)
)
dev.off()



# RUN METANEIGHBOR BETWEEN BATCHES

# summarized experiment object (only human)
table(rownames(hu.metadata)==rownames(seur.hu@meta.data), useNA="ifany")
hu.metadata$batch <- seur.hu@meta.data$Batch
head(hu.metadata)

se.hu <- SummarizedExperiment(assays=hu.counts,
                              colData=hu.metadata)

# run metaneighbor (fast/slow)
mtn_hu_fast <- MetaNeighborUS(var_genes=total.hvg,
                              dat=se.hu,
                              study_id=se.hu$batch,
                              cell_type=se.hu$clusters_MAIT,
                              fast_version=T)
mtn_hu_slow <- MetaNeighborUS(var_genes=total.hvg,
                              dat=se.hu,
                              study_id=se.hu$batch,
                              cell_type=se.hu$clusters_MAIT,
                              fast_version=F)
# reorder
order_mait_test <- c(rep("|MAIT_c0", 3),
                     rep("|MAIT_c1", 3),
                     rep("|MAIT_c2", 3),
                     rep("|MAIT_c3", 3),
                     rep("|MAIT_c4", 3),
                     rep("|MAIT_c5", 3),
                     rep("|MAIT_c6", 3))
order_mait_test <- paste0(rep(c("B", "C", "D"), 7), order_mait_test)
mtn_hu_fast <- mtn_hu_fast[order_mait_test, order_mait_test]
mtn_hu_slow <- mtn_hu_slow[order_mait_test, order_mait_test]

# plot fast
heatmap.2(mtn_hu_fast,
          # trace
          trace="none",
          # dendrogram
          dendrogram="none",
          Rowv=F,
          Colv=F,
          # superimpose a density histogram on color key
          density.info="none",
          # color scale
          col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100)),
          breaks=seq(0,1,length=101),
          key.xlab = "AUROC",
          key.title="",
          keysize = 1.2,
          # text labels
          main="fast version",
          cexRow=0.6,
          cexCol=0.6,
          # margins
          margins=c(6,6)
)

# plot slow
heatmap.2(mtn_hu_slow,
          # trace
          trace="none",
          # dendrogram
          # dendrogram="none",
          # Rowv=F,
          # Colv=F,
          # superimpose a density histogram on color key
          density.info="none",
          # color scale
          col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100)),
          breaks=seq(0,1,length=101),
          key.xlab = "AUROC",
          key.title="",
          keysize = 1.2,
          # text labels
          main="slow version",
          cexRow=0.6,
          cexCol=0.6,
          # margins
          margins=c(6,6)
)
