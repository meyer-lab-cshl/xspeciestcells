# Purpose: Check co-expression of SCENIC regulons
# Author: Salomé Carcy
# Date: 2023-10-25




# **************
# 1. IMPORT ####
# **************

# Import librairies
library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)
library(Seurat)

# Import data
# seur.human <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_08_28_23.RDS")
seur.human <- readRDS("./data/clean_data/seurat_human_integrated_object_23_12_01.rds")
DimPlot(seur.human, group.by="new_clusters", reduction="UMAP_50")



# *****************
# 2. FUNCTIONS ####
# *****************





# ****************
# 3. ANALYSIS ####
# ****************

#___________________________
## 3.1. Re-run clustering with high resolution ####
seur.highclus <- seur.human
seur.highclus <- FindNeighbors(seur.highclus, reduction="harmony")
seur.highclus <- FindClusters(seur.highclus, reduction="harmony", resolution=75)
table(seur.highclus$RNA_snn_res.75)

## /end ####


#___________________________
## 3.2. Group cells per cluster, get average expression of each gene by cluster ####
counts.norm <- seur.highclus@assays$RNA@data
metadata    <- seur.highclus@meta.data

# Subset metadata only to information we're interested in (tissue, cell identity, batch, and clusters)
# metadata <- metadata[,c("Tissue", "cell.ident", "Batch", "new_clusters", "RNA_snn_res.75", "Donor")]
metadata <- metadata[,c("tissue", "tcell_lineage", "batch_id", "clusters_integrated_data", "RNA_snn_res.75", "donor_id")]
colnames(metadata) <- c("tissue_id", "cell_id", "batch_id", "cluster_id", "cluster_highres", "donor_id")
head(metadata) # sanity check

# AGGREGATE COUNTS
counts.agg <- matrix(0,nrow=nrow(counts.norm), ncol=length(unique(metadata$cluster_highres)))
rownames(counts.agg) <- rownames(counts.norm)
colnames(counts.agg) <- sort(unique(metadata$cluster_highres))
counts.agg[1:5,1:5]

for(i in sort(unique(metadata$cluster_highres))){
  print(i)
  # subset
  cells.temp <- rownames(metadata[metadata$cluster_highres==i,])
  counts.temp <- counts.norm[,cells.temp]
  # get gene averages
  genesavg.temp <- rowMeans(counts.temp)
  # print(genesavg.temp[1:5])
  counts.agg[,i] <- genesavg.temp
}
counts.agg[1:5,1:5]
dim(counts.agg)
## /end ####


#___________________________
## 3.3. Spearman correlation between all genes ####
cor_matrix_allgenes <- cor(t(counts.agg), method="spearman") # takes about 5min to compute

# make lower triangle full of NAs
cor_matrix_allgenes[lower.tri(cor_matrix_allgenes)] <- NA
dim(cor_matrix_allgenes)
cor_matrix_allgenes[1:5,1:5]


# normalize coexpression network like cococonet
cor_matrix_allgenes_norm <- reshape2::melt(cor_matrix_allgenes) %>%
  as_tibble() %>%
  rename(spearman_cor=value) %>%
  filter(!is.na(spearman_cor))

cor_matrix_allgenes_norm <- cor_matrix_allgenes_norm %>%
  # rank standardize correlation
  arrange(Var1, spearman_cor) %>%
  group_by(Var1) %>%
  mutate(rank=rank(spearman_cor),
         rank_max=max(rank)) %>%
  # normalize by dividing by maximum rank
  mutate(rank_normalized=rank/rank_max) %>%
  ungroup() %>%
  select(-rank_max)

## /end ####


#___________________________
## 3.4. Spearman correlation between TFs and their gene in regulon ####

# Import regulons
regulons <- read.csv("./data/human-PBMC/HumanData_32_Coexpression/nes2_maskdropouts_full_obj_seurat_filtered_harmony_08_28_23_raw_counts_100_multirun_reg100-target80_regulons_pruned_final.csv",
                     header=F, fill=T, col.names=1:3373)
regulons[1:20,1:5]
dim(regulons)
table(regulons$X2)

rownames(regulons) <- regulons$X1
regulons[,1:2] <- NULL
regulons[regulons==""] <- NA
regulons[1:17,1:10]


# Get spearman correlation for each regulon
# df.cor <- data.frame()
# for(tf in rownames(regulons)){
#   print(tf)
#   # get the genes in that regulon
#   genes.temp <- regulons[tf,]
#   genes.temp <- genes.temp[!is.na(genes.temp)]
#   # print(length(genes.temp))
#   # if a gene is not found in seurat object
#   if(length(genes.temp[!genes.temp %in% rownames(seur.highclus)]) != 0){
#     print(paste0("genes", genes.temp[!genes.temp %in% rownames(seur.highclus)], "not found in seurat object"))
#   }
#   # get the spearman correlation
#   df.temp <- data.frame("regulon"=tf,
#                         "gene"=genes.temp)
#   df.temp$spearman_cor <- apply(df.temp, 1, function(x) stats::cor(counts.agg[x[1],], counts.agg[x[2],], method="spearman"))
#   # add to big dataframe
#   df.cor <- rbind(df.cor, df.temp)
# }
# table(df.cor$regulon)
# sanity check
# table(is.na(df.cor))
# df.cor %>%
#   filter(regulon=="ZBTB16")

# # Plot their spearman correlation
# hist(df.cor$spearman_cor, breaks=100)
# df.cor %>%
#   filter(spearman_cor==1)
# 
# df.cor %>%
#   ggplot(aes(y=regulon, x=spearman_cor))+
#   geom_boxplot(outlier.shape=NA)+
#   geom_point(size=0.1)
# ggsave("./data/human-PBMC/HumanData_32_Coexpression/coexpression_per_regulon.jpeg", width=5, height=40)
# 
# 
# # plot spearman correlation of regulons kept by Laurent
# df.cor %>%
#   filter(regulon %in% c("RUNX2", "NR1D2", "CREM", "NFE2L2", "FLI1", "MYBL1", "RORA", "XBP1", "CEBPD", "FOSL2", "BHLHE40", "MAF", "PRDM1", "EOMES")) %>%
#   ggplot(aes(y=regulon, x=spearman_cor))+
#   geom_boxplot(outlier.shape=NA)+
#   geom_point(size=0.1)


# Get TF-target genes in a long DF
regulons_londf <- data.frame()
for(tf in rownames(regulons)){
  print(tf)
  genes.temp <- regulons[tf,]
  genes.temp <- genes.temp[!is.na(genes.temp)]
  # if a gene is not found in seurat object
  if(length(genes.temp[!genes.temp %in% rownames(seur.highclus)]) != 0){
    print(paste0("genes", genes.temp[!genes.temp %in% rownames(seur.highclus)], "not found in seurat object"))
  }
  df.temp <- data.frame("regulon"=tf, "targetgene"=genes.temp)
  regulons_londf <- rbind(regulons_londf, df.temp)
}
# sanity checks
# table(regulons_londf$regulon)
# table(is.na(regulons_londf))
# regulons_londf %>%
#   filter(regulon=="ZBTB16")

# subset the correlation matrix to our regulons
regulons_londf <- regulons_londf %>%
  rowwise() %>%
  mutate(genepair=paste0(sort(c(regulon, targetgene)), collapse="|"))

test <- cor_matrix_allgenes_norm
test$Var1 <- as.character(test$Var1)
test$Var2 <- as.character(test$Var2)


head(test)
test <- cor_matrix_allgenes_norm %>%
  mutate(Var1=as.character(Var1),
         Var2=as.character(Var2)) %>%
  rowwise() %>%
  mutate(genepair=paste0(sort(c(Var1, Var2)), collapse="|")) %>%
  head(10)

head(cor_matrix_allgenes_norm)



## /end ####


#___________________________
## 3.5. Check the null distribution of spearman correlations in CoCoCoNet ####
counts.agg[1:5,1:5]
nrow(counts.agg)

# get matrix from coconet: https://labshare.cshl.edu/shares/gillislab/resource/CoCoBLAST/CoCoCoNet/networks/human_HC_AggNet.hdf5
matrix.coconet <- rhdf5::H5Fopen("./data/human-PBMC/HumanData_32_Coexpression/human_HC_AggNet.hdf5")
table(matrix.coconet$col==matrix.coconet$row) # rows and cols are same

coconet.matrix <- matrix.coconet$agg
rownames(coconet.matrix) <- matrix.coconet$row
colnames(coconet.matrix) <- matrix.coconet$col
coconet.matrix[1:5,1:5]

# keep only genes of interest
regulons.genes <- unique(unique(df.cor$gene), unique(df.cor$regulon)) # 10,399 genes
coconet.genes  <- matrix.coconet$col

ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
genemap <- biomaRt::getBM(attributes=c('ensembl_gene_id', "external_gene_name"),
                          filters = 'external_gene_name',
                          values = rownames(seur.human), # get all genes in our seurat object
                          mart = ensembl) %>%
  dplyr::rename(gene_id = ensembl_gene_id,
                symbol = external_gene_name) %>%
  distinct() %>%
  filter(symbol %in% regulons.genes)
table(regulons.genes %in% genemap$symbol) # 264 genes from our regulons not found
table(genemap$gene_id %in% coconet.genes) # 9,944 of genes within our regulons have a corresponding ENSEMBL ID

genemap <- genemap %>%
  filter(gene_id %in% coconet.genes)
length(unique(genemap$gene_id))
length(unique(genemap$symbol))
genemap <- genemap %>%
  filter(!duplicated(symbol))


# subset coconet matrix to genes of interest
coconet.matrix.sub <- coconet.matrix[genemap$gene_id, genemap$gene_id]
# table(colnames(coconet.matrix.sub)==genemap$gene_id, useNA="ifany")
rownames(coconet.matrix.sub) <- genemap$symbol
colnames(coconet.matrix.sub) <- genemap$symbol
# dim(coconet.matrix.sub)
# coconet.matrix.sub[1:5,1:5]



# get spearman correlation for each regulon (now with only genes from genemap)
regulons2 <- regulons[rownames(regulons) %in% genemap$symbol,]
# df.cor %>%
#   filter(regulon %in% rownames(regulons)[!rownames(regulons) %in% genemap$symbol]) %>%
#   ggplot(aes(y=regulon, x=spearman_cor))+
#   geom_boxplot(outlier.shape=NA)+
#   geom_point(size=0.1) # the regulons we're removing are not really the ones that have a high coexpression with their target genes anyway...
regulons2[1:5,1:5]

df.cor2 <- data.frame()
for(tf in rownames(regulons2)){
  print(tf)
  # get the genes in that regulon
  genes.temp <- regulons2[tf,]
  genes.temp <- genes.temp[!is.na(genes.temp)]
  genes.temp <- genes.temp[genes.temp %in% genemap$symbol]
  # print(length(genes.temp))
  # if a gene is not found in seurat object
  if(length(genes.temp[!genes.temp %in% rownames(seur.highclus)]) != 0){
    print(paste0("genes", genes.temp[!genes.temp %in% rownames(seur.highclus)], "not found in seurat object"))
  }
  if(length(genes.temp)>=1){
    # get the spearman correlation
    df.temp <- data.frame("regulon"=tf,
                          "gene"=genes.temp)
    df.temp$spearman_cor <- apply(df.temp, 1, function(x) stats::cor(counts.agg[x[1],], counts.agg[x[2],], method="spearman"))
    df.temp$coconet_cor <- apply(df.temp, 1, function(x) coconet.matrix.sub[x[1], x[2]])
    # add to big dataframe
    df.cor2 <- rbind(df.cor2, df.temp)
  }
}
table(df.cor2$regulon)
# sanity check
# table(is.na(df.cor2))
# df.cor2 %>%
  # filter(regulon=="ZBTB16")
apply(test, 1, function(x) coconet.matrix.sub[x[1], x[2]])



# Plot their spearman correlation
hist(df.cor2$spearman_cor, breaks=100, main="correlations btw TF-targets in our data")
hist(df.cor2$coconet_cor, breaks=100, main="correlations btw TF-targets across human tissues (cococonet data)")
df.cor2 %>%
  ggplot(aes(y=regulon, x=spearman_cor))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(size=0.1)+
  xlim(c(-1,1))+
  labs(title="correlation on our data (pseudobulked to small clusters)") |
  df.cor2 %>%
  ggplot(aes(y=regulon, x=coconet_cor))+
  geom_boxplot(outlier.shape=NA, color="#756bb1")+
  geom_point(size=0.1)+
  xlim(c(-1,1))+
  labs(x="spearman_cor", y="", title="correlation across human tissues (cococonet data)")
ggsave("./data/human-PBMC/HumanData_32_Coexpression/coexpression_per_regulon_withbaseline.jpeg", width=10, height=35)


# plot spearman correlation of regulons kept by Laurent
df.cor %>%
  filter(regulon %in% c("RUNX2", "NR1D2", "CREM", "NFE2L2", "FLI1", "MYBL1", "RORA", "XBP1", "CEBPD", "FOSL2", "BHLHE40", "MAF", "PRDM1", "EOMES")) %>%
  ggplot(aes(y=regulon, x=spearman_cor))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(size=0.1)





## /end ####

