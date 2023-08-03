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




# *************************
# 2. DE ANALYSIS (LRT) ####
# *************************

# Question: Are there genes that explain more lineage than anything else (batch, cluster/cell state)?

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
  # filter(batch_id %in% c("E", "F", "I")) %>%
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

# Re-run PCA
ddsPCA <- function(counts.batchcorrect, metadf){
  rv <- matrixStats::rowVars(counts.batchcorrect) # variance of each gene
  select_rv <- order(rv, decreasing = TRUE)[seq_len(500)] # get the positions of the top 500 most variable genes?...
  pca <- prcomp(t(counts.batchcorrect[select_rv,])) # run pca on top 500 HVG
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  counts_pca <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                           batch_id = metadf$batchdonor_id,
                           lineage_id  = metadf$lineage_id,
                           cluster_id = metadf$cluster_id)
  return(counts_pca)
}
counts_pca <- ddsPCA(counts.batchcorrect = counts_batchcorrect, metadf=metadf.deseq)

# Plot PCA on batch-corrected counts
ggplot(counts_pca, aes(x = PC1, y = PC2, color = lineage_id, shape=batch_id)) +
  geom_point(size = 4) +
  # xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  # ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  coord_fixed(expand=TRUE)+
  scale_color_manual(values=cols_lineages)
ggplot(counts_pca, aes(x = PC1, y = PC2, color = cluster_id, shape=batch_id)) +
  geom_point(size = 4) +
  # xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  # ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
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
print(res)

# Keep only significant DE genes
genes.sig <- res %>%
  data.frame() %>%
  filter(padj<0.01)
dim(genes.sig) # 209 genes (FDR 0.05), 122 genes (FDR 0.01)



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
# ggsave("./plots/allPBMCs_batchE-I_DEgenescluster_padj0_01.jpeg", clusters$plot + labs(title="Batches E-I | padj<0.01"), width=8, height=8)




# *********************************************
# 3. DE ANALYSIS LINEAGE-SPECIFIC DE GENES ####
# *********************************************

# *******************
## 3.1. Function ####
cellTypeSignature <- function(cell.1, padj_max=0.05, log2FC_min=0.5, shrinkage=T, DEmethod="LRT"){
  
  # All cell lineages
  allcells <- c("CD8", "CD4", "MAIT", "NKT", "GD")
  allothercells <- allcells[allcells != cell.1]
  
  # Initialize
  genes.up <- list()
  
  # Loop
  for (cell.2 in allothercells){
    print(paste0(cell.1, "vs", cell.2))
    
    # 1. Create DDS object with only contrast of interest
    metadf.temp <- metadf.deseq %>% filter(lineage_id %in% c(cell.1, cell.2)) %>% mutate(lineage_id=factor(lineage_id, levels=unique(lineage_id)))
    counts.temp <- counts.deseq[,rownames(metadf.temp)]
    # table(colnames(counts.temp)==rownames(metadf.temp))
    dds.temp <- DESeqDataSetFromMatrix(counts.temp, colData = metadf.temp, design = ~ batchdonor_id + lineage_id + cluster_id)
    
    # 2. Normalize, batch-correct counts and plot PCA
    rld.temp <- rlog(dds.temp, blind=TRUE)
    counts_batchcorrect.temp <- limma::removeBatchEffect(x=assay(rld.temp),
                                                         batch=metadf.temp$batchdonor_id,
                                                         design=model.matrix(~ lineage_id + cluster_id, metadf.temp))
    counts_pca.temp <- ddsPCA(counts.batchcorrect = counts_batchcorrect.temp, metadf=metadf.temp)
    # Plot PCA on batch-corrected counts
    ggplot(counts_pca.temp, aes(x = PC1, y = PC2, color = lineage_id, shape=batch_id)) +
      geom_point(size = 4) +
      # xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      # ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      coord_fixed(expand=TRUE)+
      scale_color_manual(values=cols_lineages)
    ggplot(counts_pca.temp, aes(x = PC1, y = PC2, color = cluster_id, shape=batch_id)) +
      geom_point(size = 4) +
      # xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      # ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      coord_fixed(expand=TRUE)+
      scale_color_manual(values=cols_integrated)
    
    # 3. Run DESeq2 (LRT test): parameters based on recommendations
    # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis
    dds.temp <- DESeq(dds.temp, test=DEmethod, reduced=~batchdonor_id+cluster_id, useT=TRUE, minmu=1e-6, minReplicatesForReplace = Inf)
    # plotDispEsts(dds.temp) # Plot dispersion estimates
    contrastvec <- c("lineage_id", cell.1, cell.2)
    res.temp <- results(dds.temp, contrast = contrastvec, alpha = 0.05)
    if(shrinkage==T){
      print("Shrinking LFC")
      res.temp <- lfcShrink(dds.temp, type="ashr", contrast = contrastvec, res=res.temp)
    }
    # print(res.temp)
    genes.sig.temp <- res.temp %>% data.frame() %>% filter(padj<padj_max & log2FoldChange>log2FC_min)
    print(paste0("# upregulated genes:", nrow(genes.sig.temp)))
    
    # 4. Save list of DE genes
    genes.up[[paste0(cell.1, "vs", cell.2)]] <- rownames(genes.sig.temp)
  } # FOR LOOP
  
  return(genes.up)
}


# *************************************************
## 3.2. Get lineage-specific upregulated genes ####

# Get upregulated genes
cd4up  <- cellTypeSignature(cell.1="CD4",  padj_max=0.05, log2FC_min = 0.1, shrinkage=F)
cd8up  <- cellTypeSignature(cell.1="CD8",  padj_max=0.05, log2FC_min = 0.1, shrinkage=F)
maitup <- cellTypeSignature(cell.1="MAIT", padj_max=0.05, log2FC_min = 0.1, shrinkage=F)
nktup  <- cellTypeSignature(cell.1="NKT",  padj_max=0.05, log2FC_min = 0.1, shrinkage=F)
gdtup  <- cellTypeSignature(cell.1="GD",   padj_max=0.05, log2FC_min = 0.1, shrinkage=F)

# Look at genes that are upregulated in at least 2 contrasts
plyr::count(Reduce(c, cd4up))  %>% filter(freq>=2) %>% pull(x) # 13Wald 11LRT
plyr::count(Reduce(c, cd8up))  %>% filter(freq>=3) %>% pull(x) # 24Wald 21LRT
plyr::count(Reduce(c, maitup)) %>% filter(freq>=2) %>% pull(x) # 17Wald 28LRT
plyr::count(Reduce(c, nktup))  %>% filter(freq>=2) %>% pull(x) # 6Wald  7LRT
plyr::count(Reduce(c, gdtup))  %>% filter(freq>=2) %>% pull(x) # 11Wald 12LRT

# Bind them
genes.lineage <- unique(c(plyr::count(Reduce(c, cd4up))  %>% filter(freq>=3) %>% pull(x),
                        plyr::count(Reduce(c, cd8up))  %>% filter(freq>=3) %>% pull(x),
                        plyr::count(Reduce(c, maitup)) %>% filter(freq>=3) %>% pull(x),
                        plyr::count(Reduce(c, nktup))  %>% filter(freq>=2) %>% pull(x),
                        plyr::count(Reduce(c, gdtup))  %>% filter(freq>=3) %>% pull(x)))


# ************************
## 3.3. Plot heatmap ####

# Get corrected counts (that we batch corrected with limma earlier)
counts.correc.sig <- t(counts_batchcorrect[genes.lineage,]) # keep only the genes of interest
counts.correc.sig[,1:5]
dim(counts.correc.sig)

# Reorder by lineage and cluster (not by batch)
counts.correc.sig <- counts.correc.sig[metadf.deseq %>% arrange(lineage_id, as.numeric(cluster_id)) %>% pull(sample_id),]
cols_batchdonorid <- pals::brewer.greys(length(unique(metadf.deseq$batchdonor_id)))
names(cols_batchdonorid) <- unique(metadf.deseq$batchdonor_id)

# Run pheatmap using the metadata data frame for the annotation
# pdf("./plots/allPBMCs_batchE-F-I_lineagespecific_upregulatedgenes_heatmap_padj0_05.pdf", width=12, height=10)
pheatmap::pheatmap(t(counts.correc.sig),
                   color = heat_colors,
                   scale = "row", # z-score
                   clustering_method="ward.D2",
                   cluster_rows = T,
                   cluster_cols = F,
                   border_color = NA,
                   # Columns (cell groups)
                   gaps_col=c(17,36,54,69),
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

