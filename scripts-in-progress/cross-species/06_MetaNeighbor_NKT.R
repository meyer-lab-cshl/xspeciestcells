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
seur.ms <- readRDS("./data/cross-species/00_Reproduce_UMAPs/ms_nkt_seurobj.rds")
seur.ms@meta.data$cell_type <- as.character(seur.ms@meta.data$cell_type)
# DimPlot(seur.ms)
cols_nkt <- c("#810f7c", "#8856a7", "#8c96c6", "#b3cde3", "#edf8fb")
names(cols_nkt) <- unique(seur.ms$cell_type)
SCpubr::do_DimPlot(seur.ms, 
                   group.by = "cell_type",
                   label=T,
                   label.color="black",
                   legend.position = "none",
                   repel=T,
                   plot.title = "Mouse",
                   colors.use=cols_nkt,
                   font.size = 24)
# ggsave("./data/cross-species/04_Metaneighbor/ms_umap.jpeg", width=6, height=6)


# seur.hu <- readRDS("./data/human-thymus/HumanData_12_AnalysisByLineage/seurat_thymus-NKT_2023-02-27.rds")
seur.hu <- readRDS("./data/human-thymus/HumanData_12_AnalysisByLineage/thymus.nkt_02_28_23.RDS")
colors_clusters_NKT <- c("0" = "#d8443c", "1" = "#e09351", "2" = "gold", "3" = "#74c8c3", "4" = "#5a97c1",
                         "5" = "#a40000", "6" = "#72bcd5")
Idents(seur.hu) <- seur.hu$new_clusters_NKT
SCpubr::do_DimPlot(seur.hu, 
                   group.by = "new_clusters_NKT",
                   label=T,
                   label.color="black",
                   legend.position = "none",
                   repel=T,
                   plot.title = "Human",
                   colors.use=colors_clusters_NKT,
                   font.size = 24)


ortholog.df <- read.csv("./data/cross-species/03_BiomartTable/big_ass_ortholog_table.csv")
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
table(unique(rownames(ms.counts)) %in% ortholog.df$ms_symbol_data) # 13,969 not
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
                          hu_symbol)
hu.hvg.translated <- pull(dictionary %>% filter(hu_symbol %in% hu.hvg), # not all hu HVGs are found in ortholog.df
                          hu_symbol)
total.hvg <- unique(union(ms.hvg.translated, hu.hvg.translated))
length(total.hvg) # 2,368 genes


# Keep only ms and hu genes that have 1:1 orthologs
table(unique(rownames(ms.counts)) %in% dictionary$ms_symbol_data) # 12,201 genes should have a translation
table(unique(rownames(hu.counts)) %in% dictionary$hu_symbol) # 12,201 genes should have a translation
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


# Merge everything into one
ms.metadata$study <- "Mouse"
ms.metadata <- ms.metadata[,c("cell_type", "study")]
colnames(ms.metadata)[1] <- "clusters_NKT"
head(ms.metadata)

hu.metadata$study <- "Human"
hu.metadata <- hu.metadata[,c("new_clusters_NKT", "study")]
colnames(hu.metadata)[1] <- "clusters_NKT"
hu.metadata$clusters_NKT <- paste0("iNKT_c", hu.metadata$clusters_NKT)
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
# seur.total <- merge(ms.seur, hu.seur)
# 
# # Sanity checks
# head(seur.total@meta.data)
# table(seur.total$clusters_NKT, useNA="ifany")
# table(seur.total$study, useNA="ifany")
# 
# 
# # Convert seurat count matrix to SummarizedExperiment object
# count <- SummarizedExperiment(assays=seur.total@assays[["RNA"]]@counts,
#                               colData=seur.total$clusters_NKT)

# _______________________
# METANEIGHBOR
# Run metaneighbor
mtn <- MetaNeighborUS(var_genes=total.hvg,
                      dat=count,
                      study_id=count$study,
                      cell_type=count$clusters_NKT,
                      fast_version=T)
# saveRDS(mtn, "./data/cross-species/04_Metaneighbor_nkt/nkt_ms-hu_mtnslowversion.rds")
mtn <- readRDS("./data/cross-species/04_Metaneighbor_nkt/nkt_ms-hu_mtnslowversion.rds")

# Heatmap
mtn.sub <- mtn[1:5,6:12]
mtn.sub <- mtn.sub[rev(rownames(mtn.sub)),]

jpeg("./data/cross-species/04_Metaneighbor_nkt/nkt_ms-hu_fastversion_fulltree.jpeg",
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
          main="Mouse vs Human iNKT",
          cexRow=0.6,
          cexCol=0.6,
          # margins
          margins=c(6,6))
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
  # dplyr::rename(human=Var1, mouse=Var2, auroc=value) %>%
  # add nb of cells per human cluster
  # mutate(Var1=gsub("c", "", human)) %>%
  left_join(as.data.frame(table(seur.hu$new_clusters_NKT)) %>% mutate(Var1=paste0("iNKT_c", Var1)), by="Var1") %>%
  dplyr::rename(human=Var1, Var1=Var2, auroc=value, ncells_human=Freq) %>%
  mutate(totalcells_human = dim(seur.hu)[2],
         propcells_human = ncells_human*100/totalcells_human) %>%
  # add nb of cells per mouse cluster
  left_join(as.data.frame(table(seur.ms$cell_type)), by="Var1") %>%
  dplyr::rename(mouse=Var1, ncells_mouse=Freq) %>%
  mutate(totalcells_mouse = dim(seur.ms)[2],
         propcells_mouse = ncells_mouse*100/totalcells_mouse) #%>%
  # rename human clusters
  # mutate(human=paste0("iNKT_", human))
# saveRDS(mtn.df, "./data/cross-species/04_Metaneighbor_nkt/nkt_ms-hu_mtnslowversion_DF.rds")


# PROPORTION OF HUMAN NKT CELLS IN EACH CLUSTER
bp.x <- ggplot(data=mtn.df%>% select(human,propcells_human) %>% distinct(),
               aes(x=factor(human, levels=paste0("iNKT_c", 0:6)), y=propcells_human))+
  geom_bar(stat="identity", fill="#bdbdbd") + theme_cowplot()+
  scale_x_discrete(position="top")+
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100))+
  labs(y="%cells")+
  theme(axis.text = element_text(size=15),
        axis.text.x = element_text(size = 15, angle=45, hjust=0),
        axis.ticks.y = element_blank(),
        axis.title.y=element_text(size=20),
        axis.title.x = element_blank(),
        axis.line.x=element_blank(),
        legend.position = "none")

# PROPORTION OF MOUSE NKT CELLS IN EACH CLUSTER
order_mouse <- c("Stage0", "iNKTp", "iNKT2", "iNKT17", "iNKT1")
bp.y <- ggplot(data=mtn.df%>% select(mouse,propcells_mouse) %>% distinct(),
               aes(x=factor(mouse, levels=rev(order_mouse)), y=propcells_mouse))+
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
hm.clean <- ggplot(mtn.df, aes(x=factor(human, levels=paste0("iNKT_c", 0:6)),
                               y=factor(mouse, levels=rev(order_mouse)))) +
  geom_point(aes(size = auroc, color= auroc))+
  geom_text(data=mtn.df %>% filter(auroc>0.65) %>% mutate(across("auroc", \(x) round(x,2))), aes(label=auroc), color="white")+
  scale_size_continuous(limits=c(0,1), breaks=seq(0,1, by=0.2), range = c(1, 15))+
  scale_color_gradient2(low="#d9d9d9", mid="white", high="#a50f15", midpoint=0.5, limits=c(0,1), name="AUROC", breaks=seq(0,1, by=0.2))+
  labs(x="Human clusters",y="Mouse clusters", size="AUROC")+
  theme_cowplot()+
  theme(legend.position="bottom", legend.key.width = unit(0.8, 'cm'),
        axis.text = element_text(size=15), axis.title=element_text(size=20), axis.text.x=element_text(angle=45, hjust=1))

# COMBINE
library(patchwork)
(bp.x+plot_spacer() + plot_layout(widths = c(5, 1))) / (hm.clean + bp.y + plot_layout(widths = c(5, 1))) + plot_layout(heights = c(1, 5))
# ggsave("./data/cross-species/04_Metaneighbor_nkt/nkt_ms-hu_metaneighbor_fastversion_bubbleplot.jpeg", width=9, height=9.5)



# _______________________
# METANEIGHBOR TROUBLESHOOTING FAST/SLOW
# get slow mtn
mtn_slow <- readRDS("./data/cross-species/04_Metaneighbor_nkt/nkt_ms-hu_mtnslowversion.rds")
colnames(mtn_slow)[6:12] <- paste0("Human|iNKT_", substr(colnames(mtn_slow)[6:12], 7, 9))
rownames(mtn_slow)[6:12] <- paste0("Human|iNKT_", substr(rownames(mtn_slow)[6:12], 7, 9))

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
jpeg("./data/cross-species/04_Metaneighbor_nkt/nkt_ms-hu_slowversion_fulltree_orderedlikefastversion.jpeg", width=1500, height=1500, res=200)
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
          main="Mouse vs Human iNKT",
          cexRow=0.6,
          cexCol=0.6,
          # margins
          margins=c(6,6)
)
dev.off()




# _____________________
# CORRELATION ####

# Get average expression per gene df
avgexp.sub <- averageGE(seuratobj=seur.total, geneslist=total.hvg)
head(avgexp.sub)
colnames(avgexp.sub)[1:5] <- paste0("Mouse|", colnames(avgexp.sub)[1:5])
colnames(avgexp.sub)[6:10] <- paste0("Human|", colnames(avgexp.sub)[6:10])

#7a:  Correlation
Corr.Coeff.Table = cor(avgexp.sub, method="spearman")

#7b:  Shuffle data
shuffled.cor.list = list()
pb   <- txtProgressBar(1, 100, style=3)
nPermutations <- 1000

# Calculate correlations on shuffled expression table
for (i in 1:nPermutations){
  # # Shuffle gene exp across all columns
  # shuffled = apply(avgexp.sub,1,sample)
  # shuffled = t(shuffled)
  # Shuffle gene exp by dataset
  shuffled_ms  = apply(avgexp.sub[,str_detect(colnames(avgexp.sub),"Mouse")],1,sample)
  shuffled_hu  = apply(avgexp.sub[,str_detect(colnames(avgexp.sub),"Human")],1,sample)
  # Bind back together avg exp table (shuffled)
  shuffled = cbind(t(shuffled_ms),t(shuffled_hu))
  if(ncol(shuffled)!=ncol(avgexp.sub)){ print("PROBLEM") }
  # Correlation on shuffled table
  shuffled.cor = cor(shuffled, method="spearman")
  shuffled.cor.list[[i]] = shuffled.cor
  # rm(list=c('shuffled', 'shuffled.cor'))
  rm(list=c('shuffled_ms','shuffled_hu', 'shuffled', 'shuffled.cor'))
  if ((i %% 10) ==0){
    setTxtProgressBar(pb, (i*100)/nPermutations)
  }
}

# Make empty p.value and corr mean matrixes
p.value.table = matrix(ncol=ncol(avgexp.sub), nrow = ncol(avgexp.sub))
rownames(p.value.table) = colnames(avgexp.sub)
colnames(p.value.table) = colnames(avgexp.sub)
shuffled.mean.table = p.value.table

# Get the mean "random" correlation from all permutations, and pvalue comparing the "observed" correlation with the "random" one
a = combn(1:ncol(avgexp.sub),2)
for (i in 1:ncol(a)){
  cor.scores = sapply(shuffled.cor.list,"[",a[1,i],a[2,i])
  # mean correlation score on random data
  shuffled.mean.table[a[1,i],a[2,i]] = mean(cor.scores)
  shuffled.mean.table[a[2,i],a[1,i]] = mean(cor.scores)
  # empirical p.val
  p.value = mean(abs(cor.scores)>=abs(Corr.Coeff.Table[a[1,i],a[2,i]]))
  p.value.table[a[1,i],a[2,i]] = p.value
  p.value.table[a[2,i],a[1,i]] = p.value
  rm(list=c('cor.scores','p.value'))
  setTxtProgressBar(pb, (i*100)/ncol(a))
}


# Plot
# Get the correlation and ptable with only 2 species (mouse x other)
comp_table <- Corr.Coeff.Table[1:5,6:10]
p_table <- p.value.table[1:5,6:10]

# Prepare colors
col1 <- colorRampPalette(c("darkblue", "white","darkred"))


heatmap.2(comp_table,
          # trace
          trace="none",
          # superimpose a density histogram on color key
          density.info="none",
          # color scale
          col=colorRampPalette(c("darkblue", "white","darkred")),
          breaks=seq(0.4,0.6,length=101),
          # text labels
          main="",
          cexRow=0.6,
          cexCol=0.6,
          # margins
          margins=c(8,6))

# Plot
jpeg("~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/02_CorrelationComparion/ms_hu_corrplot_laurent-orthologs.jpeg", width=1500, height=1500, res=300)
corrplot(comp_table,
         # 2 species
         order="hclust",
         hclust.method="average",
         tl.col="black",
         
         # other parameters
         tl.pos="lt", method="color",
         col=col1(200),
         # cl.lim=c(min(comp_table.intersect),max(comp_table.intersect)),
         is.corr=F, tl.cex=0.7, sig.level=0.05,
         p.mat=p_table, # significance
         insig="pch", pch=19, pch.cex=0.25, pch.col="black",
         # addrect = 3,
         main= "",mar=c(3,1,5,1),cl.align.text="l") #%>%
# corrRect(c(1,6,13, 20))
# corrRect(namesMat=c("ms_Stage0", "hu_iNKT_4", "ms_iNKT1", "hu_iNKT_3"))
dev.off()