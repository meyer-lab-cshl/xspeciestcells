###
# Purpose: Perform DE analysis btw lineages on T cells in clusters 15-17
# Date: July 2023
# Author: Salomé Carcy
###


# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(ggplot2)
library(tidyverse)
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(cowplot)
source("~/Projects/HumanThymusProject/scripts-final/colors_universal.R") # get color palettes

## 1.2. Data ####
# path.plots <- "~/Projects/HumanThymusProject/data/human-thymus/HumanData_19_CellIdentity"
seur.human <- readRDS("~/Projects/HumanThymusProject/data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
print(seur.human) # 78,607 cells (it's the whole seurat object)

# Quick visualization
DimPlot(seur.human, reduction="UMAP_50", group.by="new_clusters", cols = cols_integrated)




# *******************
# 2. DE ANALYSIS ####
# *******************

# Question: Is there any DE genes between NKT/MAIT/GDT in clusters 15-17 (GEP4)?
setwd("~/Projects/HumanThymusProject/scripts-in-progress/human-PBMC/HumanData_23_DEanalysisPBMCclusters15_17/")

# ********************************
## 2.1. Prepare DESeq2 object ####

# Create SingleCellExperiment object
counts   <- seur.human@assays$RNA@counts
metadata <- seur.human@meta.data

# Subset metadata only to information we're interested in (tissue, cell identity, batch, and clusters)
metadata <- metadata[,c("Tissue", "cell.ident", "Batch", "new_clusters", "Donor")]
colnames(metadata) <- c("tissue_id", "cell_id", "batch_id", "cluster_id", "donor_id")
head(metadata) # sanity check

# Create a "groups" df that will (1) keep only cells of interest; (2) keep the columns of interest (that will define how we group counts)
groups <- metadata %>%
  rownames_to_column("cell") %>%
  # keep only cells of interest (only batches/donors E5 and I11 contain all 5 lineages, so more thorough for DE analysis)
  filter(tissue_id=="PBMC") %>%
  filter(cluster_id %in% 15:17) %>%
  filter(donor_id %in% c(5, 11)) %>%
  filter(batch_id %in% c("E", "I")) %>%
  # keep only groups with at least 10 cells (with higher threshold we lose a lot of groups)
  group_by(cell_id, batch_id) %>%
  filter(n()>40) %>%
  ungroup() %>%
  # keep only columns of interest
  column_to_rownames("cell") %>%
  select(cell_id, batch_id)
head(groups) # sanity check


# AGGREGATE COUNTS
count.agg <- t(counts[, which(colnames(counts) %in% rownames(groups))]) # keep only cells of interest (defined in "groups") & put cells in rows (genes as columns)
nrow(groups) == nrow(count.agg) # verify nb of cells in "groups" (rows) is same nb of cells in "count.agg" (rows)
count.agg <- Matrix.utils::aggregate.Matrix(count.agg, groupings = groups, fun = "sum") # aggregate counts based on columns in "groups"

# Sanity checks
dim(count.agg)[2] == nrow(seur.human) # same nb of genes in seurat object & count.agg
count.agg[, 1:6]
table(colSums(count.agg) == 0) # check if any gene has total count of 0

# Final counts for DESeq (groups as columns, genes as rows)
counts.deseq <- data.frame(t(count.agg))
counts.deseq <- counts.deseq[rowSums(counts.deseq)!=0,] # remove genes that have total count of 0


# PREPARE METADATA DF FOR DESEQ
# Have metadata df where rownames are the colnames of "counts.deseq" (the groups)
metadf.deseq <- groups %>%
  distinct() %>%
  unite(sample_id, remove=FALSE) %>%
  dplyr::rename(lineage_id=cell_id)
rownames(metadf.deseq) <- metadf.deseq$sample_id

# Check that the row names of metadf.deseq are the same as the column names of counts.deseq in order to use as input to DESeq2
nrow(metadf.deseq)==ncol(counts.deseq)
metadf.deseq <- metadf.deseq[match(colnames(counts.deseq), rownames(metadf.deseq)),] # reorder rows in metadf.deseq to match order of columns in counts.deseq
all(rownames(metadf.deseq) == colnames(counts.deseq))


# PREPARE COLOR SCALES
cols_batchid <- brewer.pal(length(unique(metadf.deseq$batch_id)), "Greys")
names(cols_batchid) <- unique(metadf.deseq$batch_id)



# CREATE DESEQ2 OBJECT
dds <- DESeqDataSetFromMatrix(counts.deseq, 
                              colData = metadf.deseq,
                              design = ~ batch_id + lineage_id)


# ****************************
## 2.2. PCA on the groups ####

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
# head(assay(rld))

# Plot PCA
# DESeq2::plotPCA(rld, intgroup = "sample_id")
# DESeq2::plotPCA(rld, intgroup = "batch_id")
# DESeq2::plotPCA(rld, intgroup = "lineage_id")

# Correct for batch effect
counts_batchcorrect <- limma::removeBatchEffect(x=assay(rld),
                                                batch=metadf.deseq$batch_id,
                                                design=model.matrix(~ lineage_id, metadf.deseq))

# Run PCA on batch-corrected counts
rv <- matrixStats::rowVars(counts_batchcorrect) # variance of each gene
select_rv <- order(rv, decreasing = TRUE)[seq_len(500)] # get the positions of the top 500 most variable genes?...
pca <- prcomp(t(counts_batchcorrect[select_rv,])) # run pca on top 500 HVG
percentVar <- pca$sdev^2/sum(pca$sdev^2)
counts_pca <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                         batch_id = metadf.deseq$batch_id,
                         lineage_id  = metadf.deseq$lineage_id)
ggplot(counts_pca, aes(x = PC1, y = PC2, color = lineage_id, shape=batch_id)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  coord_fixed(expand=TRUE)+
  scale_color_manual(values=cols_lineages)


# # Compute pairwise correlation values on the batch-effect-corrected matrix
# rld_cor <- cor(counts_batchcorrect[select_rv,])
# # Plot heatmap
# pheatmap(rld_cor,
#          annotation     = metadf.deseq,
#          annotation_row = metadf.deseq,
#          annotation_colors = list(lineage_id = cols_lineages, batch_id=cols_batchid))



# ***************************
## 2.3. Run DESeq2 (LRT) ####

# Run DESeq2 differential expression analysis
dds <- DESeq(dds, test="LRT", reduced=~batch_id)
plotDispEsts(dds) # Plot dispersion estimates


# Output results for contrast batch_id + lineage_id VS batch_id
res <- results(dds, alpha = 0.05) # alpha is FDR
# Only need to shrink log2FC if thresholding on log2FC (but here we'll threshold on padj)
# res <- lfcShrink(dds, type="ashr",
#                  # contrast = c("lineage_id", "MAIT", "GD"),
#                  res=res)
print(res)

# Keep only significant DE genes
genes.sig <- res %>%
  data.frame() %>%
  filter(padj<0.05)
dim(genes.sig)
# sum(res$padj[!is.na(res$padj)]<0.01) # 39 genes with padj<0.01 (and with min cells >10)

# Compare nb genes in common with the clusters 13-14 comparison
# genes.sig.1314 <- read.csv("~/Projects/HumanThymusProject/scripts-in-progress/human-PBMC/HumanData_23_DEanalysisPBMCclusters13_14/plots/deg_clust13-14_batchE5I11_padj0_05.csv", row.names = 1)
# table(rownames(genes.sig) %in% genes.sig.1314$gene, useNA="ifany") # 23 common, 29 different


# ***************************************
## 2.4. Plot heatmap on all lineages ####

# Get corrected counts (that we batch corrected with limma earlier)
counts.correc.sig <- t(counts_batchcorrect[rownames(genes.sig),]) # keep only the genes of interest
# counts.correc.sig <- counts.correc.sig[grep("GD|MAIT|NKT", rownames(counts.correc.sig), value=T),] # if want to plot only innate T cells
counts.correc.sig[,1:5]
dim(counts.correc.sig)

# Set a color palette
heat_colors <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(100))


# Run pheatmap using the metadata data frame for the annotation
# pdf("./plots/All_heatmap_clust15-17_batchE5I11_min40cells_padj0_05.pdf", width=10, height=12)
pheatmap::pheatmap(t(counts.correc.sig),
                   color = heat_colors,
                   scale = "row", # z-score
                   clustering_method="ward.D2",
                   cluster_rows = T,
                   # cutree_rows = 5,
                   cluster_cols = T,
                   border_color = NA,
                   # Columns (cell groups)
                   show_colnames = T,
                   fontsize_col = 10,
                   annotation_col =metadf.deseq %>% select(c(lineage_id, batch_id)),
                   annotation_colors = list(lineage_id = cols_lineages,
                                            batch_id=cols_batchid),
                   # Rows (genes)
                   show_rownames=T,
                   fontsize_row=8,
                   # title
                   main="PBMC clusters 15-17 (batches E5, I11, >40cells)")
# dev.off()


# # ************************************
# ## 2.5. Identifying gene clusters ####
# 
# # Tutorial: https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
# 
# # Use degPatterns function to show gene clusters across sample groups
# # BiocManager::install("DEGreport")
# library(DEGreport)
# all(rownames(metadf.deseq) == rownames(counts.correc.sig)) # sanity check
# metadf.deseq$lineage_id <- factor(metadf.deseq$lineage_id, levels=unique(metadf.deseq$lineage_id))
# clusters <- degPatterns(t(counts.correc.sig),
#                         metadata = metadf.deseq,
#                         minc=5, # minimum nb of genes per group
#                         time = "lineage_id",
#                         col=NULL)
# # ggsave(plot=clusters$plot + labs(title="PBMC clusters 15-17 (batches E5, I11): DE genes btw all lineages"),
#        # filename="./plots/degpatterns_clust15-17_batchEI_padj0_05.jpeg", width=8, height=6)
# 
# # Extract gene lists
# cluster_groups <- clusters$df
# # write.csv(cluster_groups, "./plots/degpatterns_clust13-14_padj0_01_listgenespergroup.csv")
# clusters$df %>% filter(cluster==1) # CD4 and NKT group
# clusters$df %>% filter(cluster==3) # CD8-GD-MAIT group
# clusters$df %>% filter(cluster==5) # CD8-MAIT group
# 
# write.csv(genes.sig %>%
#             rownames_to_column("gene") %>%
#             dplyr::select(gene, padj) %>%
#             left_join(cluster_groups, by=join_by("gene"=="genes")) %>%
#             dplyr::rename(degPatterns_group=cluster),
#           "./plots/deg_clust15-17_batchE5I11_padj0_05.csv")
