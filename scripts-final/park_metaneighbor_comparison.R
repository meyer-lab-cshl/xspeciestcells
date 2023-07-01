# Purpose: Comparison of our thymic clusters to Park thymic clusters
# Author: Salomé Carcy
# Date: June 2023


# 1.IMPORT ---------------------- 
## 1.1. Import librairies -------
library(ggplot2)
# library(gplots)
library(cowplot)
library(tidyverse)
library(dplyr)
library(Seurat)
# library(SeuratData)
# library(SeuratDisk)
# library(RColorBrewer)
# library(harmony)
# library(corrplot)
# library(pals)
# library(VennDiagram)
library(MetaNeighbor)
library(SummarizedExperiment)


## 1.2. Import data -------
# get universal color scales
source("./scripts-final/colors_universal.R")

# import our data
seur.gapin <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
print(seur.gapin) # 78,607 cells
SCpubr::do_DimPlot(seur.gapin, 
                   group.by = "new_clusters",
                   label=T,
                   label.color="black",
                   legend.position = "none",
                   repel=T,
                   colors.use=cols_integrated,
                   font.size = 24)

# import park thymocyte data
# downloaded the seurat object from https://cellxgene.cziscience.com/collections/de13e3e2-23b6-40ed-a413-e9e12d7d3910
# I may have updated the cell_type metadata column with the metadata from anndata object (not sure anymore)
seur.park.thymocytes <- readRDS("./data/human-thymus/HumanData_16_ParkIntegration/park_thymus_seu_gene_names.rds")
print(seur.park.thymocytes) # 76,994 cells
SCpubr::do_DimPlot(seur.park.thymocytes,
                   group.by = "cell_type",
                   label=T,
                   label.color="black",
                   legend.position = "none",
                   repel=T,
                   colors.use=cols_park,
                   font.size = 24)




# 2.METANEIGHBOR ---------------------- 
## 2.1. Clean-up & prepare data -------

# Prepare metadata
Idents(seur.gapin) <- "Tissue"
seur.thymus <- subset(seur.gapin, ident="Thymus")
SCpubr::do_DimPlot(seur.thymus, 
                   group.by = "new_clusters",
                   label=T,
                   label.color="black",
                   legend.position = "none",
                   repel=T,
                   colors.use=cols_integrated,
                   font.size = 24)
# ggsave("./data/human-thymus/HumanData_16_ParkIntegration/gapin_umap_thymus.jpeg", width=10, height=10)
# call the clusters "c0", "c1", etc.
colnames(seur.thymus@meta.data)[17] <- "cell_type"
seur.thymus@meta.data[,"cell_type"] <- paste0("c",seur.thymus@meta.data[,"cell_type"])

seur.thymus@meta.data$study          <- "gapin_data"
seur.park.thymocytes@meta.data$study <- "park_data"


# Keep only genes that are found in both datasets
genes_to_keep <- intersect(rownames(seur.thymus), rownames(seur.park.thymocytes))
length(genes_to_keep) # 13,988 genes in common
seur.thymus          <- seur.thymus[genes_to_keep,]
seur.park.thymocytes <- seur.park.thymocytes[genes_to_keep,]


# Get HVGs from each dataset to run metaneighbor on
seur.thymus          <- FindVariableFeatures(seur.thymus, nfeatures = 2000)
seur.park.thymocytes <- FindVariableFeatures(seur.park.thymocytes, nfeatures = 2000)
ggVennDiagram::ggVennDiagram(
  x=list(VariableFeatures(seur.thymus), VariableFeatures(seur.park.thymocytes)),
  category.names = c("Gapin HVGs" , "Park HVGs"),
  label_alpha=0,
  label="count"
) # 894 common HVGs

# get union of HVGs from our data & park data
hvg.union <- unique(c(VariableFeatures(seur.thymus), VariableFeatures(seur.park.thymocytes)))
length(hvg.union) # 3,106 genes
# table(hvg.union %in% rownames(seur.thymus)) # sanity check all TRUE
# table(hvg.union %in% rownames(seur.park.thymocytes)) # sanity check all TRUE


# Merge seurat objects
seur.merged <- merge(seur.thymus, y=seur.park.thymocytes)
# sanity checks
# ncol(seur.thymus) + ncol(seur.park.thymocytes) # total 114,363 cells
# length(union(rownames(seur.thymus), rownames(seur.park.thymocytes))) # total 13,988 genes
# print(seur.merged)


# Convert seurat count matrix to SummarizedExperiment object
table(seur.merged$cell_type, useNA="ifany")
table(seur.merged$study, useNA="ifany")
count <- SummarizedExperiment(assays=seur.merged@assays[["RNA"]]@counts,
                              colData=seur.merged@meta.data[,c("cell_type", "study")])

# Export to .rds files to run on Elzar (HPC)
# VariableFeatures(seur.merged) <- hvg.union
# saveRDS(seur.merged, "./data/human-thymus/HumanData_16_ParkIntegration/metaneighbor_elzar/seuratobj_merged.rds")
# saveRDS(count, "./data/human-thymus/HumanData_16_ParkIntegration/metaneighbor_elzar/summarizedexpobj_merged.rds")



## 2.2. Run metaneighbor -------

mtn <- MetaNeighborUS(var_genes=hvg.union,
                      dat=count,
                      study_id=seur.merged$study,
                      cell_type=seur.merged$cell_type,
                      one_vs_best = FALSE,
                      fast_version=FALSE)

# checkpoint save
saveRDS(mtn, "./data/human-thymus/HumanData_16_ParkIntegration/mtn_park_2023-06-28.rds")
mtn <- readRDS("./data/human-thymus/HumanData_16_ParkIntegration/mtn_park.rds")

# Heatmap
jpeg("~/Projects/HumanThymusProject/data/human-thymus/HumanData_16_ParkIntegration/metaneighbor_thymocytesonly_onevsbest.jpeg",
     width=1000, height=1000, res=150)
heatmap.2(mtn2,
          # trace
          trace="none",
          # Rowv=FALSE,
          # Colv=FALSE,
          # superimpose a density histogram on color key
          density.info="none",
          # color scale
          col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100)),
          breaks=seq(0,1,length=101),
          #na.rm=F,
          na.color="grey",
          # text labels
          main="Comparison with Park thymocytes",
          cexRow=0.6,
          cexCol=0.6,
          # margins
          margins=c(7,7))
dev.off()


# Reformat AUROC into df
library(reshape2)
library(ggrepel)
mtn.df <- melt(mtn2)
mtn.df <- mtn.df %>%
  filter(str_detect(Var1,"gapin_data")) %>%
  mutate(Var1 = gsub("gapin_data\\|", "", Var1)) %>%
  filter(str_detect(Var2, "park_data")) %>%
  mutate(Var2 = gsub("park_data\\|", "", Var2))

interesting_clusters <- mtn.df %>%
  group_by(Var1) %>%
  top_n(1, value)

ggplot(mtn.df, aes(x=factor(Var1, level=paste0("c", 0:17)),
                   y=value,
                   color=factor(Var2, level=names(cols_park))))+
  geom_hline(yintercept=0.9, linetype="dashed", color="grey")+
  geom_hline(yintercept=0.5, linetype="dashed", color="grey")+
  geom_point()+
  geom_text_repel(data = interesting_clusters,
                  aes(label = Var2), nudge_y=0.1, show.legend=FALSE)+
  scale_color_manual(values=cols_park, name="Park clusters")+
  ylim(c(0,1.2))+
  labs(x="Gapin clusters", y="AUROC")+
  theme_cowplot()
# ggsave("~/Projects/HumanThymusProject/data/human-thymus/HumanData_16_ParkIntegration/metaneighbor_thymocytesonly_ggplot.jpeg", width=8, height=5)

# Get number of cells per gapin cluster and per park cluster
mtn.df <- mtn.df %>%
  left_join(as.data.frame(table(seur.thymus$cell_type)), by="Var1") %>%
  dplyr::rename(ncells_gapin=Freq) %>%
  left_join(as.data.frame(table(seur.park.thymocytes$cell_type)) %>% dplyr::rename(Var2=Var1), by="Var2") %>%
  dplyr::rename(ncells_park=Freq)

bp.x <- ggplot(data=mtn.df, aes(x=factor(Var1, levels=paste0("c", 0:17)), y=ncells_gapin))+
  geom_bar(stat="identity", fill="#bdbdbd") + theme_cowplot()+
  scale_x_discrete(position="top")+
  labs(y="#cells")+
  theme(#axis.title.y = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 15), axis.title.x = element_blank(), axis.line.x=element_blank(),
    legend.position = "none") #+
# scale_fill_distiller(name = "Value", palette = "Greys", direction = 1)

# Change order of Park clusters
order_park <- c("DN(early)", "DN(P)", "DN(Q)", "DP(P)", "DP(Q)", "αβT(entry)", "T(agonist)", "CD8αα(I)", "CD8αα(II)",
                "γδT", "Treg(diff)", "Treg", "CD4+T", "CD8+T", "Th17", "NKT")
bp.y <- ggplot(data=mtn.df, aes(x=factor(Var2, levels=rev(order_park)), y=ncells_park))+
  geom_bar(stat="identity", fill="#bdbdbd") +
  scale_x_discrete(position="top") + labs(y="#cells")+ coord_flip() + theme_cowplot()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size=15), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.y=element_blank(),
        legend.position = "none")# +
# scale_fill_distiller(name = "Value", palette = "Greys", direction = 1)

# Bubble plot
library(scales)
hm.clean <- ggplot(mtn.df, aes(x=factor(Var1, levels=paste0("c", 0:17)),
                               y=factor(Var2, levels=rev(order_park)), size = value,
                               color=value)) +
  geom_point(alpha=0.7)+
  scale_colour_gradient2(low="#d9d9d9", mid="white", high="#a50f15", midpoint=0.5, limits=c(0,1), name="AUROC")+
  labs(x="",y="Park clusters", size="AUROC")+
  theme_cowplot()+
  theme(legend.position="bottom", legend.key.width = unit(1, 'cm'))
