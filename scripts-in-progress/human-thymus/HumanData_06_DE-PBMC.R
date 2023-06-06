###
# Purpose: Perform DE analysis on PBMC data (group cells by lineage & cluster)
# Date: June 2023
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
source("./scripts-final/colors_universal.R") # get color palettes

## 1.2. Data ####
# path.plots <- "~/Projects/HumanThymusProject/data/human-thymus/HumanData_19_CellIdentity"
seur.human <- readRDS("~/Projects/HumanThymusProject/data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
print(seur.human) # 78,607 cells (it's the whole seurat object)

# Quick visualization
DimPlot(seur.human, reduction="UMAP_50", group.by="new_clusters", cols = cols_integrated)




# *******************
# 2. DE ANALYSIS ####
# *******************

# Question: in PBMCs, do cells differ by cell lineage (NKT, MAIT, etc.) or by cell state?... And which genes define their identity?

# ********************************
## 2.1. Prepare DESeq2 object ####

# Create SingleCellExperiment object
counts   <- seur.human@assays$RNA@counts
metadata <- seur.human@meta.data

# Subset metadata only to information we're interested in (tissue, cell identity, batch, and clusters)
metadata <- metadata[,c("Tissue", "cell.ident", "Batch", "new_clusters")]
colnames(metadata) <- c("tissue_id", "cell_id", "batch_id", "cluster_id")
head(metadata) # sanity check

# Create a "groups" df that will (1) keep only cells of interest; (2) keep the columns of interest (that will define how we group counts)
groups <- metadata %>%
  rownames_to_column("cell") %>%
  # keep only cells of interest
  filter(tissue_id=="PBMC") %>%
  # filter(cluster_id %in% 9:17) %>%
  # keep only groups with at least 100 cells
  group_by(cluster_id, cell_id, batch_id) %>%
  filter(n()>100) %>%
  ungroup() %>%
  # keep only columns of interest
  column_to_rownames("cell") %>%
  select(cell_id, batch_id, cluster_id) %>%
  mutate(cluster_id = factor(cluster_id)) # remove levels that are absent
head(groups) # sanity check


# AGGREGATE COUNTS
count.agg <- t(counts[, which(colnames(counts) %in% rownames(groups))]) # keep only cells of interest (defined in "groups") & put cells in rows (genes as columns)
nrow(groups) == nrow(count.agg) # verify nb of cells in "groups" (rows) is same nb of cells in "count.agg" (rows)
count.agg <- Matrix.utils::aggregate.Matrix(count.agg, groupings = groups, fun = "sum") # aggregate counts based on columns in "groups"

# Sanity checks
dim(count.agg)[2] == nrow(seur.human) # same nb of genes in seurat object & count.agg
count.agg[1:6, 1:6]
table(colSums(count.agg) == 0) # check if any gene has total count of 0

# Final counts for DESeq (groups as columns, genes as rows)
counts.deseq <- data.frame(t(count.agg))
counts.deseq <- counts.deseq[rowSums(counts.deseq)!=0,] # remove genes that have total count of 0


# PREPARE METADATA DF FOR DESEQ
# Have metadata df where rownames are the colnames of "counts.deseq" (the groups)
metadf.deseq <- groups %>%
  distinct() %>%
  unite(sample_id, remove=FALSE)# %>%
# mutate(gene_program=ifelse(cluster_id %in% c(9,10,11), "naive",
#                            ifelse(cluster_id==12, "Tcm_GEP6",
#                                   ifelse(cluster_id %in% c(13,14), "Th17_GEP1",
#                                          ifelse(cluster_id %in% 15:17, "Temra_GEP4", NA)))))
rownames(metadf.deseq) <- metadf.deseq$sample_id

# Check that the row names of metadf.deseq are the same as the column names of counts.deseq in order to use as input to DESeq2
nrow(metadf.deseq)==ncol(counts.deseq)
metadf.deseq <- metadf.deseq[match(colnames(counts.deseq), rownames(metadf.deseq)),] # reorder rows in metadf.deseq to match order of columns in counts.deseq
all(rownames(metadf.deseq) == colnames(counts.deseq))


# PREPARE COLOR SCALES
matcol <- metadf.deseq %>%
  select(cell_id, batch_id, cluster_id) %>%
  relocate(batch_id) %>%
  dplyr::rename(lineage_id=cell_id) %>%
  mutate(cell_state=ifelse(cluster_id %in% c(3, 6, 9:11), "Tnaive",
                           ifelse(cluster_id==12, "Tcm",
                                  ifelse(cluster_id %in% 13:14, "Th17",
                                         ifelse(cluster_id %in% 15:17, "Temra",
                                            ifelse(cluster_id==7, "Treg", "?"))))))
cols_cellstate <- c("Tnaive"= "#b3e2cd",
                   "Tcm"   = "#f4cae4",
                   "Th17"  = "#cbd5e8",
                   "Temra" = "#fdcdac",
                   "Treg" = "#fbb4ae")

cols_batchid <- brewer.pal(5, "Greys")
names(cols_batchid) <- unique(matcol$batch_id)



# CREATE DESEQ2 OBJECT
dds <- DESeqDataSetFromMatrix(counts.deseq, 
                              colData = metadf.deseq,
                              design = ~ batch_id + lineage_id + cluster_id)


# ****************************
## 2.2. PCA on the groups ####

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
# head(assay(rld))

# Plot PCA
# DESeq2::plotPCA(rld, intgroup = "sample_id")
DESeq2::plotPCA(rld, intgroup = "batch_id")
DESeq2::plotPCA(rld, intgroup = "lineage_id")

# Correct for batch effect
counts_batchcorrect <- limma::removeBatchEffect(x=assay(rld),
                                                batch=metadf.deseq$batch_id,
                                                design=model.matrix(~ lineage_id + cluster_id, metadf.deseq))

# Re-run PCA
rv <- matrixStats::rowVars(counts_batchcorrect) # variance of each gene
select_rv <- order(rv, decreasing = TRUE)[seq_len(500)] # get the positions of the top 500 most variable genes?...
pca <- prcomp(t(counts_batchcorrect[select_rv,])) # run pca on top 500 HVG
percentVar <- pca$sdev^2/sum(pca$sdev^2)
counts_pca <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                         batch_id = metadf.deseq$batch_id,
                         lineage_id  = metadf.deseq$lineage_id,
                         cluster_id = metadf.deseq$cluster_id)

# Plot PCA on batch-corrected counts
ggplot(counts_pca, aes(x = PC1, y = PC2, color = lineage_id, shape=batch_id)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  coord_fixed(expand=TRUE)+
  scale_color_manual(values=cols_lineages)
ggplot(counts_pca, aes(x = PC1, y = PC2, color = cluster_id, shape=batch_id)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  coord_fixed(expand=TRUE)+
  scale_color_manual(values=cols_integrated)


# # Compute pairwise correlation values on the batch-effect-corrected matrix
# rld_cor <- cor(counts_batchcorrect[select_rv,])
# # Plot heatmap
# pheatmap(rld_cor,
#          annotation     = matcol,
#          annotation_row = matcol,
#          annotation_colors = list(lineage_id = cols_lineages, batch_id=col_batchid, cluster_id=cols_integrated))



# ***************************
## 2.3. Run DESeq2 (LRT) ####

# Run DESeq2 differential expression analysis
dds <- DESeq(dds, test="LRT", reduced=~batch_id)
dds <- dds[which(mcols(dds)$fullBetaConv),] # remove 3 genes that didn't converge

# Plot dispersion estimates
plotDispEsts(dds)

# Output results for contrast batch_id + lineage_id + cluster_id VS batch_id
res <- results(dds,
               # contrast = contrast,
               alpha = 0.05)
res <- lfcShrink(dds, type="ashr",
                 # contrast = contrast,
                 res=res)
print(res)

# Keep only significant DE genes (don't filter on log2FC because that's based on the contrast)
genes.sig <- res %>%
  data.frame() %>%
  filter(padj < 0.001)
dim(genes.sig) # 4796 genes



# ***********************
## 2.4. Plot heatmap ####

# Get corrected counts (that we batch corrected with limma earlier)
counts.correc.sig <- t(counts_batchcorrect[rownames(genes.sig),]) # keep only the genes of interest
counts.correc.sig[1:5,1:5]
dim(counts.correc.sig) # 4796 genes

# Set a color palette
heat_colors <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(100))


# Run pheatmap using the metadata data frame for the annotation
jpeg("./data/human-thymus/HumanData_19_CellIdentity/LRT_BatchCellCluster_vs_Batch/heatmap_PBMC_LRTest_batch-cellid-clusterid_limma_corrected_adjpval01_allclusters.jpeg", height = 3000, width = 3500, res = 300)
# svg("./data/human-thymus/HumanData_19_CellIdentity/LRT_BatchCellCluster_vs_Batch/heatmap_PBMC_LRTest_batch-cellid-clusterid_limma_corrected_adjpval01.svg", height = 10, width = 11)
phtmap <- pheatmap(counts.correc.sig,
         color = heat_colors,
         scale = "column", # z-score
         cluster_rows = T,
         cluster_cols = T,
         # clustering_method="complete",
         border_color = NA,
         # Create gaps to separate CD4/NKT
         # gaps_row=c(2,6, 11, 14),
         # cutree_rows = 5,
         cutree_cols = 10,
         # Change rows text format
         show_rownames = T,
         fontsize_row = 10,
         annotation_row=matcol,
         annotation_colors = list(lineage_id = cols_lineages,
                                  batch_id=cols_batchid,
                                  cluster_id=cols_integrated[levels(metadf.deseq$cluster_id)],
                                  cell_state=cols_cellstate),
         # Change columns
         angle_col=45,
         # fontsize_col = 1,
         show_colnames=F,
         # title
         main="PBMC - LRT test - batchid+cellid+clusterid vs batchid (batch corrected, all groups with >100 cells)")
dev.off()


# Get row dendrogram and flip it
row_dend <- phtmap[[1]]
row_order <- rownames(counts.correc.sig[phtmap$tree_row[["order"]],])
plot(phtmap$tree_row)
plot(row_dend)

bloc1 <- row_order[1:11]
bloc2 <- row_order[54:65]
# bloc3 <- row_order[35:53]
bloc3a <- row_order[52:53]
bloc3b <- row_order[35:51]
# bloc4 <- row_order[23:34]
bloc4a <- row_order[23]
bloc4b <- row_order[33:34]
bloc4c <- row_order[31:32]
bloc4d <- row_order[24:30]
bloc5 <- row_order[12:22]
# row_order_new <- c(bloc1, bloc2, bloc3, bloc4, bloc5)
row_order_new <- c(bloc1, bloc2, bloc3a, bloc3b, bloc4a, bloc4b, bloc4c, bloc4d, bloc5)
# length(row_order)==length(row_order_new)
# table(unique(row_order) %in% unique(row_order_new))

row_dend_new <- dendextend::rotate(row_dend, order=row_order_new)
# plot(row_dend_new)


# Re-plot pheatmap
# jpeg("./data/human-thymus/HumanData_19_CellIdentity/LRT_BatchCellCluster_vs_Batch/heatmap_PBMC_LRTest_batch-cellid-clusterid_limma_corrected_adjpval001_allclusters3.jpeg", height = 3000, width = 3500, res = 300)
pheatmap(counts.correc.sig,
         color = heat_colors,
         scale = "column", # z-score
         cluster_rows = as.hclust(row_dend_new),
         cluster_cols = T,
         # clustering_method="complete",
         border_color = NA,
         # Create gaps to separate CD4/NKT
         # gaps_row=c(2,6, 11, 14),
         # cutree_rows = 5,
         cutree_cols = 6,
         # Change rows text format
         show_rownames = T,
         fontsize_row = 10,
         annotation_row=matcol,
         annotation_colors = list(batch_id=cols_batchid,
                                  lineage_id = cols_lineages,
                                  cluster_id=cols_integrated[levels(metadf.deseq$cluster_id)],
                                  cell_state=cols_cellstate),
         # Change columns
         angle_col=45,
         # fontsize_col = 1,
         show_colnames=F,
         # title
         main="PBMC - LRT test - batchid+cellid+clusterid vs batchid (batch corrected, all groups with >100 cells)")
# dev.off()





