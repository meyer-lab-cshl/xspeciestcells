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

# Import data
seur.ms <- readRDS("./data/cross-species/00_Reproduce_UMAPs/ms_gdt_seurobj_lee.rds")
# order_mouse <- c("cKIT+ DN1", "DN2", "DN3", "Pre-selected GD", "Post-selected GD", "Pan GD (mainly CD24+)", "CD122+ GD", "CD24- GD")
order_mouse <- c("Tγδp", "immature Tγδ1/17", "immature Tγδ17", "Tγδ1", "Tγδ17")
cols_gdt <- brewer.pal(length(unique(seur.ms$gd_clusters)), "Blues")
names(cols_gdt) <- unique(seur.ms$gd_clusters)
SCpubr::do_DimPlot(seur.ms, 
                   group.by = "gd_clusters",
                   label=T,
                   label.color="black",
                   legend.position = "none",
                   repel=T,
                   plot.title = "Mouse",
                   colors.use=cols_gdt,
                   font.size = 24)
# ggsave("./data/cross-species/04_Metaneighbor_gdt/ms_umap_lee.jpeg", width=6, height=6)


seur.hu <- readRDS("./data/human-thymus/HumanData_12_AnalysisByLineage/thymus.GD.03_09_23.RDS")
colors_clusters_GDT <- c("GD_c0" = "#d8443c", "GD_c1" = "#e09351", "GD_c2" = "gold", "GD_c3"="#bcbddc", "GD_c4" = "#74c8c3", "GD_c5" = "#6a51a3",
                         "GD_c6" = "#addd8e", "GD_c7" = "#bdbdbd")
# Idents(seur.hu) <- seur.hu$new_clusters_GDT
SCpubr::do_DimPlot(seur.hu, 
                   group.by = "new_clusters_GD",
                   label=T,
                   label.color="black",
                   legend.position = "none",
                   repel=T,
                   plot.title = "Human",
                   colors.use=colors_clusters_GDT,
                   font.size = 24)
# ggsave("./data/cross-species/04_Metaneighbor_gdt/hu_umap2.jpeg", width=6, height=6)


ortholog.df <- read.csv("./data/cross-species/03_BiomartTable/big_ass_ortholog_table.csv")
length(unique(ortholog.df$ms_symbol_data)) # 17,085 ms genes
length(unique(ortholog.df$hu_symbol)) # 17,152 hu genes



# Get counts, HVGs and metadata
ms.hvg <- VariableFeatures(FindVariableFeatures(seur.ms, nfeatures=2000))
length(ms.hvg)
ms.counts <- seur.ms[["RNA"]]@counts
ms.metadata <- seur.ms@meta.data

hu.hvg <- VariableFeatures(seur.hu)
length(hu.hvg)
hu.counts <- seur.hu[["RNA"]]@counts
hu.metadata <- seur.hu@meta.data


# First, check whether genes can all be found in the ortholog table (many can't be found because I removed genes with orthology confidence=0)
table(unique(rownames(ms.counts)) %in% unique(ortholog.df$ms_ensemblID)) # 12,401 not
table(unique(rownames(hu.counts)) %in% ortholog.df$hu_symbol) # 4,747 not


# Subset the ortholog table to only genes that we can "translate"
dictionary <- ortholog.df %>%
  as_tibble() %>%
  ## Intersection
  filter(ms_ensemblID %in% unique(rownames(ms.counts)) & hu_symbol %in% unique(rownames(hu.counts))) %>%
  ## Remove any symbols that are NAs
  filter(!is.na(ms_ensemblID)) %>%
  filter(!is.na(hu_symbol)) %>%
  ## Keep only 1:1 orthologs
  group_by(ms_ensemblID) %>% filter(n_distinct(hu_symbol) == 1) %>% ungroup() %>%
  group_by(hu_symbol) %>% filter(n_distinct(ms_ensemblID) == 1) %>% ungroup()
dim(dictionary) # 11,619 genes


# Translate the mouse HVGs into "human gene" language
ms.hvg.translated <- pull(dictionary %>% filter(ms_ensemblID %in% ms.hvg), # not all ms HVGs are found in ortholog.df
                          hu_symbol)
hu.hvg.translated <- pull(dictionary %>% filter(hu_symbol %in% hu.hvg), # not all hu HVGs are found in ortholog.df
                          hu_symbol)
total.hvg <- unique(union(ms.hvg.translated, hu.hvg.translated))
length(total.hvg) # 2,002 genes (if 2000 HVG each) / 4872 (if 5000 HVGs each)
table(total.hvg %in% unique(ms.hvg.translated), useNA="ifany") # 1,230 come from mouse
table(total.hvg %in% unique(hu.hvg.translated), useNA="ifany") # 1,164 come from human


# Keep only ms and hu genes that have 1:1 orthologs
table(unique(rownames(ms.counts)) %in% dictionary$ms_ensemblID) # 11,331 genes should have a translation
table(unique(rownames(hu.counts)) %in% dictionary$hu_symbol) # 11,331 genes should have a translation
ms.counts <- ms.counts[rownames(ms.counts) %in% dictionary$ms_ensemblID,]
hu.counts <- hu.counts[rownames(hu.counts) %in% dictionary$hu_symbol,]
nrow(ms.counts)==nrow(hu.counts) # should be TRUE


# Translate the mouse genes in count table into "human gene"
ms.dict <- dictionary %>%
  filter(ms_ensemblID %in% rownames(ms.counts)) %>%
  select(ms_ensemblID, hu_symbol, hu_orthology_confidence) %>%
  # distinct() %>%
  # group_by(ms_ensemblID) %>% filter(n_distinct(hu_symbol)>1)
  distinct(ms_ensemblID, .keep_all=T)
ms.dict <- ms.dict[match(rownames(ms.counts), ms.dict$ms_ensemblID),]
table(ms.dict$ms_ensemblID == rownames(ms.counts)) # should be all TRUE
table(is.na(ms.dict$hu_symbol)) # should have no NAs
# Translate
rownames(ms.counts) <- ms.dict$hu_symbol


# Merge everything into one
ms.metadata$study <- "Mouse"
ms.metadata <- ms.metadata[,c("gd_clusters", "study")]
colnames(ms.metadata)[1] <- "clusters_GDT"
head(ms.metadata)

hu.metadata$study <- "Human"
hu.metadata <- hu.metadata[,c("new_clusters_GD", "study")]
colnames(hu.metadata)[1] <- "clusters_GDT"
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

# _______________________
# METANEIGHBOR
# Run metaneighbor
mtn <- MetaNeighborUS(var_genes=total.hvg,
                      dat=count,
                      study_id=count$study,
                      cell_type=count$clusters_GDT,
                      fast_version=T)
# saveRDS(mtn, "./data/cross-species/04_Metaneighbor_gdt/gdt_mslee-hu_mtnslowversion_2023-08-03.rds")
mtn <- readRDS("./data/cross-species/04_Metaneighbor_gdt/gdt_mslee-hu_mtnslowversion_2023-08-03.rds")

# Heatmap
mtn.sub <- mtn[1:5,6:13]
mtn.sub <- mtn.sub[,paste0("Human|GD_c", 0:7)]
mtn.sub <- mtn.sub[paste0("Mouse|", order_mouse),]

# jpeg("./data/cross-species/04_Metaneighbor_gdt/gd_hvg_test/fulltree_subcluters_hvg2000_intersect392.jpeg",
#      width=1500, height=1500, res=200)
jpeg("./data/cross-species/04_Metaneighbor_gdt/gdt_mslee-hu_fastversion_fulltree.jpeg",
     width=1500, height=1500, res=200)
heatmap.2(mtn, # mtn.sub[,order(colnames(mtn.sub))],
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
          main="Mouse vs Human GDT (union: 2002 HVGs)",
          cexRow=0.6,
          cexCol=0.6,
          # margins
          margins=c(9,9))
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
  dplyr::rename(auroc=value) %>%
  # add nb of cells per human cluster
  left_join(as.data.frame(table(seur.hu$new_clusters_GD)), by="Var1") %>%
  dplyr::rename(human=Var1, ncells_human=Freq) %>%
  mutate(totalcells_human = dim(seur.hu)[2],
         propcells_human = ncells_human*100/totalcells_human) %>%
  # add nb of cells per mouse cluster
  dplyr::rename(Var1=Var2) %>%
  left_join(as.data.frame(table(seur.ms$gd_clusters)), by="Var1") %>%
  dplyr::rename(mouse=Var1, ncells_mouse=Freq) %>%
  mutate(totalcells_mouse = dim(seur.ms)[2],
         propcells_mouse = ncells_mouse*100/totalcells_mouse)
# saveRDS(mtn.df, "./data/cross-species/04_Metaneighbor_gdt/gdt_mslee-hu_mtnslowversion_DF_2023-08-03.rds")


# PROPORTION OF HUMAN GDT CELLS IN EACH CLUSTER
bp.x <- ggplot(data=mtn.df %>% select(human,propcells_human) %>% distinct(),
               aes(x=factor(human, levels=paste0("GD_c", 0:7)), y=propcells_human))+
  geom_bar(stat="identity", fill="#bdbdbd") + theme_cowplot()+
  scale_x_discrete(position="top")+
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100))+
  labs(y="%cells")+
  theme(axis.text = element_text(size=15),
        axis.text.x=element_text(angle=45, hjust=0),
        axis.ticks.y = element_blank(),
        axis.title.y=element_text(size=20),
        axis.title.x = element_blank(),
        axis.line.x=element_blank(),
        legend.position = "none")

# PROPORTION OF MOUSE GDT CELLS IN EACH CLUSTER
bp.y <- ggplot(data=mtn.df%>% select(mouse,propcells_mouse) %>% distinct(),
               aes(x=factor(mouse, levels=rev(order_mouse)), y=propcells_mouse))+
               # aes(x=factor(mouse, levels=sort(unique(seur.ms$gd_clusters), decreasing=T)), y=propcells_mouse))+
  geom_bar(stat="identity", fill="#bdbdbd") +
  scale_x_discrete(position="top") +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100))+
  labs(y="%cells")+ coord_flip() + theme_cowplot()+
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size=15),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size=20),
        axis.line.y=element_blank(),
        legend.position = "none")

# BUBBLE PLOT
# library(scales)
hm.clean <- ggplot(mtn.df, aes(x=factor(human, levels=paste0("GD_c", 0:7)),
                               y=factor(mouse, levels=rev(order_mouse)))) +
                               # y=factor(mouse, levels=sort(unique(seur.ms$gd_clusters), decreasing=T)))) +
  geom_point(aes(size = auroc, color= auroc))+
  geom_text(data=mtn.df %>% filter(auroc>0.65) %>% mutate(across("auroc", \(x) round(x,2))), aes(label=auroc), color="white")+
  scale_size_continuous(limits=c(0,1), breaks=seq(0,1, by=0.2), range = c(1, 15))+
  scale_color_gradient2(low="#d9d9d9", mid="white", high="#a50f15", midpoint=0.5, limits=c(0,1), name="AUROC", breaks=seq(0,1, by=0.2))+
  labs(x="Human clusters",y="Mouse clusters", size="AUROC")+
  theme_cowplot()+
  theme(legend.position="bottom", legend.key.width = unit(0.8, 'cm'),
        axis.text = element_text(size=15), axis.text.x = element_text(angle=45, hjust=1), axis.title=element_text(size=20))

# COMBINE
library(patchwork)
(bp.x+plot_spacer() + plot_layout(widths = c(5, 1))) / (hm.clean + bp.y + plot_layout(widths = c(5, 1))) + plot_layout(heights = c(1, 5))
# ggsave("./data/cross-species/04_Metaneighbor_gdt/gdt_mslee-hu_metaneighbor_fastversion_bubbleplot.jpeg", width=11, height=9)
# ggsave("./data/cross-species/04_Metaneighbor_gdt/gd_hvg_test/bubbleplot_subclusters_hvg2000_intersect392.jpeg", width=12, height=9)





# _______________________
# METANEIGHBOR TROUBLESHOOTING FAST/SLOW
# get slow mtn
mtn_slow <- readRDS("./data/cross-species/04_Metaneighbor_gdt/gdt_mslee-hu_mtnslowversion_2023-08-03.rds")

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
jpeg("./data/cross-species/04_Metaneighbor_gdt/gdt_mslee-hu_slowversion_fulltree_orderedlikefastversion.jpeg", width=1500, height=1500, res=200)
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
          main="Mouse vs Human GD",
          cexRow=0.6,
          cexCol=0.6,
          # margins
          margins=c(7,7)
)
dev.off()

