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

# Question: Is there any DE genes between NKT/MAIT/GDT in clusters 13-14 (GEP1)?

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
  # keep only groups with at least 100 cells
  group_by(cell_id, batch_id, cluster_id) %>%
  filter(n()>100) %>%
  ungroup() %>%
  # keep only columns of interest
  column_to_rownames("cell") %>%
  select(cell_id, batch_id, cluster_id) %>%
  mutate(cluster_id=as.character(cluster_id))
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
  unite(sample_id, remove=FALSE) %>%
  dplyr::rename(lineage_id=cell_id)
rownames(metadf.deseq) <- metadf.deseq$sample_id

# Check that the row names of metadf.deseq are the same as the column names of counts.deseq in order to use as input to DESeq2
nrow(metadf.deseq)==ncol(counts.deseq)
metadf.deseq <- metadf.deseq[match(colnames(counts.deseq), rownames(metadf.deseq)),] # reorder rows in metadf.deseq to match order of columns in counts.deseq
all(rownames(metadf.deseq) == colnames(counts.deseq))


# PREPARE COLOR SCALES
# matcol <- metadf.deseq %>%
#   select(lineage_id, batch_id, cluster_id) %>%
#   relocate(batch_id) %>%
#   mutate(cell_state=ifelse(cluster_id %in% c(3, 6, 9:11), "Tnaive",
#                            ifelse(cluster_id==12, "Tcm",
#                                   ifelse(cluster_id %in% 13:14, "Th17",
#                                          ifelse(cluster_id %in% 15:17, "Temra",
#                                                 ifelse(cluster_id==7, "Treg", "?"))))))
# tabl(matcol$cell_state)
cols_batchid <- brewer.pal(5, "Greys")
names(cols_batchid) <- unique(metadf.deseq$batch_id)



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
                         lineage_id  = metadf.deseq$lineage_id)

# Plot PCA on batch-corrected counts
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
dds <- DESeq(dds, test="LRT", reduced=~batch_id+cluster_id)
# dds <- DESeq(dds, test="Wald")
plotDispEsts(dds) # Plot dispersion estimates


# Output results for contrast batch_id + lineage_id VS batch_id
res <- results(dds,
               # contrast = c("lineage_id", "MAIT", "GD"),
               alpha = 0.05)
# res <- lfcShrink(dds, type="ashr",
#                  # contrast = c("lineage_id", "MAIT", "GD"),
#                  res=res)
print(res)

# Keep only significant DE genes
genes.sig <- res %>%
  data.frame() %>%
  filter(padj<0.05)
dim(genes.sig) # 799 genes



# ***************************************
## 2.4. Plot heatmap on all lineages ####

# Get corrected counts (that we batch corrected with limma earlier)
counts.correc.sig <- t(counts_batchcorrect[rownames(genes.sig),]) # keep only the genes of interest
# counts.correc.sig <- counts.correc.sig[grep("GD|MAIT|NKT", rownames(counts.correc.sig), value=T),]
counts.correc.sig[,1:5]
dim(counts.correc.sig)

# Set a color palette
heat_colors <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(100))


# Run pheatmap using the metadata data frame for the annotation
pdf("./plots/All_heatmap_padj0_05.pdf", width=10, height=12)
pheatmap::pheatmap(t(counts.correc.sig),
                   color = heat_colors,
                   scale = "row", # z-score
                   # clustering_method="ward.D2",
                   cluster_rows = T,
                   cluster_cols = T,
                   border_color = NA,
                   # Columns (cell groups)
                   show_colnames = T,
                   fontsize_col = 10,
                   annotation_col =metadf.deseq %>% select(c(lineage_id, batch_id, cluster_id)),
                   annotation_colors = list(lineage_id = cols_lineages,
                                            batch_id=cols_batchid,
                                            cluster_id=cols_integrated),
                   # Rows (genes)
                   show_rownames=F,
                   # fontsize_row=8,
                   # title
                   main="PBMC clusters 13-14: DE genes btw all lineages")
dev.off()


