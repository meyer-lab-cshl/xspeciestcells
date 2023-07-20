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
metadata <- metadata[,c("Tissue", "cell.ident", "Batch", "new_clusters", "Donor")]
colnames(metadata) <- c("tissue_id", "cell_id", "batch_id", "cluster_id", "donor_id")
head(metadata) # sanity check

# Create a "groups" df that will (1) keep only cells of interest; (2) keep the columns of interest (that will define how we group counts)
groups <- metadata %>%
  rownames_to_column("cell") %>%
  # keep only cells of interest
  filter(tissue_id=="PBMC") %>%
  filter(batch_id %in% c("E", "I") & donor_id %in% c("5", "11")) %>%
  # create new variable that integrates batch and donor
  mutate(batchdonor_id=paste0(batch_id, donor_id)) %>%
  # keep only groups with at least 100 cells
  # group_by(cell_id, batch_id, cluster_id, donor_id) %>%
  group_by(cell_id, batchdonor_id, cluster_id) %>%
  filter(n()>50) %>%
  ungroup() %>%
  # keep only columns of interest
  column_to_rownames("cell") %>%
  dplyr::select(cell_id, batchdonor_id, cluster_id) %>%
  # adapt a few variables
  mutate(cluster_id=as.character(cluster_id)) %>%
  dplyr::rename(lineage_id=cell_id)
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
  unite(sample_id, remove=FALSE)
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
# cols_batchdonorid <- pals::brewer.greys(12)
# names(cols_batchdonorid) <- unique(metadf.deseq$batchdonor_id)
cols_batchdonorid <- c("#BFBFBF", "#393939")
names(cols_batchdonorid) <- unique(metadf.deseq$batchdonor_id)


# CREATE DESEQ2 OBJECT
dds <- DESeqDataSetFromMatrix(counts.deseq, 
                              colData = metadf.deseq,
                              design = ~ batchdonor_id + lineage_id + cluster_id)


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
                                                batch=metadf.deseq$batchdonor_id,
                                                design=model.matrix(~ lineage_id + cluster_id, metadf.deseq))
                                                # design=model.matrix(~ lineage_id + cluster_id, metadf.deseq))

# Re-run PCA
rv <- matrixStats::rowVars(counts_batchcorrect) # variance of each gene
select_rv <- order(rv, decreasing = TRUE)[seq_len(500)] # get the positions of the top 500 most variable genes?...
pca <- prcomp(t(counts_batchcorrect[select_rv,])) # run pca on top 500 HVG
percentVar <- pca$sdev^2/sum(pca$sdev^2)
counts_pca <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                         batch_id = metadf.deseq$batchdonor_id,
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
#          annotation     = metadf.deseq,
#          annotation_row = metadf.deseq,
#          annotation_colors = list(lineage_id = cols_lineages, batch_id=cols_batchid))



# ***************************
## 2.3. Run DESeq2 (LRT) ####

# Run DESeq2 differential expression analysis
dds <- DESeq(dds, test="LRT", reduced=~batchdonor_id+cluster_id)
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
  filter(padj<0.01)
dim(genes.sig) # 632 genes



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
# pdf("./plots/allPBMCs_batchE-I_heatmap_padj0_01.pdf", width=12, height=10)
pheatmap::pheatmap(t(counts.correc.sig),
                   color = heat_colors,
                   scale = "row", # z-score
                   clustering_method="ward.D2",
                   cluster_rows = T,
                   cluster_cols = T,
                   border_color = NA,
                   # Columns (cell groups)
                   show_colnames = T,
                   fontsize_col = 10,
                   annotation_col = metadf.deseq %>% dplyr::select(c(lineage_id, batchdonor_id, cluster_id)),
                   annotation_colors = list(lineage_id = cols_lineages,
                                            batchdonor_id=cols_batchdonorid,
                                            cluster_id=cols_integrated[names(cols_integrated) %in% unique(metadf.deseq$cluster_id)]),
                   # Rows (genes)
                   show_rownames=T,
                   fontsize_row=6,
                   # title
                   main="PBMC: DE genes btw all lineages")
# dev.off()


# Use degPatterns function to show gene clusters across sample groups
library(DEGreport)
all(rownames(metadf.deseq) == colnames(t(counts.correc.sig))) # sanity check
metadf.deseq$lineage_id <- factor(metadf.deseq$lineage_id, levels=unique(metadf.deseq$lineage_id))
clusters <- degPatterns(t(counts.correc.sig),
                        metadata = metadf.deseq,
                        minc=5, # minimum nb of genes per group
                        time = "lineage_id",
                        col=NULL)



# **************************************************
## 2.5. Cell lineage-specific upregulated genes ####


# Function to contrast
cellTypeSignature <- function(cell, padj_max=0.05, log2FC_min=1){
  
  # All cell types
  allcells <- c("CD8", "CD4", "MAIT", "NKT", "GD")
  allothercells <- allcells[allcells != cell]
  
  # Initialize
  list.up <- list()
  genesup <- c()
  
  # Loop
  for (cellid in allothercells){
    print(cellid)
    # Initialize
    out       <- NULL
    result.df <- NULL
    # Set contrast
    contrast <- c("lineage_id", cell, cellid)
    # Get results
    out <- results(dds, contrast = contrast, alpha = 0.05)
    out <- lfcShrink(dds, type="ashr", contrast = contrast, res=out)
    # Get log2FC
    result.df <- out %>%
      data.frame() %>%
      rownames_to_column(var="gene") %>%
      as_tibble() %>%
      filter(padj < padj_max) %>%
      filter(log2FoldChange > log2FC_min)
    print(paste0("# upregulated genes:", nrow(result.df)))
    # Save DF in list
    list.up[[paste0(cell, "vs", cellid)]] <- result.df
    
    # Get common upregulated genes
    length_list <- length(list.up)
    if (length_list==1){ genesup <- list.up[[1]]$gene }
    else if (length_list>1){
      genesup <- intersect(genesup, list.up[[length_list]]$gene)
    }
    print(paste0("# common upregulated genes:", length(genesup)))
  }
  
  return(list("list"=list.up, "genesup"=genesup))
}




# Get CD8-upregulated genes
cd8up <- cellTypeSignature(cell="CD8", log2FC_min=0.5)
length(cd8up$genesup) # 11

# Get CD4-upregulated genes
cd4up <- cellTypeSignature(cell="CD4", log2FC_min=0.5)
cd4up$genesup # 18

# Get MAIT-upregulated genes
maitup <- cellTypeSignature(cell="MAIT", log2FC_min=0.5)
length(maitup$genesup) # 22

# Get NKT-upregulated genes
nktup <- cellTypeSignature(cell="NKT", log2FC_min=0.5)
nktup$genesup # 5

# Get GDT-upregulated genes
gdup <- cellTypeSignature(cell="GD", log2FC_min=0.5)
gdup$genesup # 7 (most difference with CD4)



# Extract normalized counts (from DESeq2 object) of the significant genes
# counts.norm.sig2 <- t(data.frame(counts(dds, normalized = TRUE))[c(cd8up$genesup, cd4up$genesup, maitup$genesup, nktup$genesup, gdup$genesup),])
counts.norm.sig2 <- assay(rld)[c(cd8up$genesup, cd4up$genesup, maitup$genesup, nktup$genesup, gdup$genesup),]
# counts.norm.sig2[1:5,1:5]
counts.corr.sig2 <- t(limma::removeBatchEffect(counts.norm.sig2, metadf.deseq$batch_id))
counts.corr.sig2[1:5,1:5]
dim(counts.corr.sig2) # 63 genes

# Run pheatmap using the metadata data frame for the annotation
# jpeg(file.path(path.plots, "LRT_BatchCell_vs_Batch/heatmap_PBMC_LRTest_cellidgenes_limma_corrected.jpeg"), height = 1500, width = 2500, res = 300)
pheatmap(counts.corr.sig2,
         color = heat_colors,
         scale = "column", # z-score
         cluster_rows = F,
         cluster_cols = T,
         border_color = NA,
         # Create gaps to separate CD4/NKT
         gaps_row=c(2, 6, 11, 14),
         cutree_rows = 5,
         cutree_cols = 5,
         # Change rows text format
         show_rownames = T,
         fontsize_row = 10,
         # annotation_row=matcol,
         # Change columns
         angle_col=45,
         fontsize_col = 4,
         show_colnames=T,
         # title
         main="PBMC - LRT test - cell type signature genes")
# dev.off()
