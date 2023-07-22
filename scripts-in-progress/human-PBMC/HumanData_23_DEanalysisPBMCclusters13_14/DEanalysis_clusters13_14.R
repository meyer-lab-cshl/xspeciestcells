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
  filter(cluster_id %in% 13:14) %>%
  # filter(cell_id %in% c("MAIT", "GD")) %>%
  # keep only groups with at least 100 cells
  group_by(cell_id, batch_id) %>%
  filter(n()>100) %>%
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
cols_batchid <- brewer.pal(5, "Greys")
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
  filter(padj<0.01)
dim(genes.sig) # 116 genes
sum(res$padj[!is.na(res$padj)]<0.01) # 79 genes with padj<0.01


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
# pdf("./plots/All_heatmap_clust13-14_padj0_01.pdf", width=10, height=12)
pheatmap::pheatmap(t(counts.correc.sig),
                             color = heat_colors,
                             scale = "row", # z-score
                             # clustering_method="ward.D2",
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
                             main="PBMC clusters 13-14: DE genes btw all lineages")
# dev.off()


## Plot volcano plot on MAIT vs GDT ####
interesting_genes <- result.df %>%
  # filter(gene %in% c("RUNX3", "ZBTB7B"))
  filter(abs(log2FoldChange) >= 3)

# Volcano plot
res.df <- res %>%
  data.frame()

ggplot(data.frame(res), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha=0.5) +
  geom_text_repel(data = data.frame(res) %>% rownames_to_column("gene") %>% filter(padj<0.05 & abs(log2FoldChange)>0.1),
                  aes(label = gene),
                  force = 10, max.overlaps=50,
                  nudge_y = 5) +
  xlim(c(-6,6))+
  ylim(c(0,55))
  # scale_color_manual(values=c("grey", "orange"))+
  # labs(x="log2 fold change", y="-log10 adjusted p-value", title="Thymic CD4 vs CD8 cells") +
  # theme_cowplot()+
  # theme(legend.position = "none",
  #       plot.title = element_text(size = rel(1.5), hjust = 0.5),
  #       axis.title = element_text(size = rel(1.25))) 

## end volcano plot ####

## Check IL12RB2 expression in MAIT vs GDT
# logcounts <- seur.human@assays$RNA@data
# logcounts <- logcounts["IL12RB2", ]
# cells <- rownames(groups[groups$cell_id %in% c("MAIT", "GD"),])
# logcounts <- data.frame("IL12RB2"=logcounts[cells])
# logcounts$lineage <- str_replace(substr(rownames(logcounts), start=0, stop=4), "_.*", "")
# head(logcounts)
# ggplot(logcounts, aes(x=lineage, y=IL12RB2))+
#   # geom_boxplot(outlier.shape=NA)+
#   geom_violin()+
#   geom_jitter(width=0.1)
  

df <- data.frame("lognorm"=counts_batchcorrect["IL12RB2", grep("MAIT|GD", colnames(counts_batchcorrect), value=T)],
                 "lineage"=c(rep("GD",4), rep("MAIT", 3)),
                 "batch" = c("E", "G", "H", "I", "E", "F", "I"))
ggplot(df, aes(x=lineage, y=lognorm))+
  geom_boxplot()+
  geom_point(aes(color=batch), size=4)+
  scale_color_manual(values=RColorBrewer::brewer.pal(5, "Set2"))+
  labs(x="", title="IL12RB2 expression in PBMCs in clusters 13-14")+
  theme_cowplot()+
  theme(panel.background = element_rect(fill="#f0f0f0"),
        title = element_text(size=10))
# ggsave("./plots/il12rb2_gdmait.jpeg", width=5, height=6)




# ************************************************
## 2.5. Plot heatmap on receptors of interest ####

# Plot heatmap on receptors of interest
# Check for IL12, IL18, IL15, and IL17 receptors
receptors <- c("IL18R1", "IL18RAP", "IL18BP", #IL18BP is a neutralizer of IL18
               "IL23R", # dimerizes with IL12BR1
               "IL12RB1", "IL12RB2",
               "IL17RE", "IL17RC", "IL17RA",
               "IFNGR1", "IFNGR2",
               "IL2RA",
               "IL6R")
res %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  filter(gene %in% receptors)

counts.correc.sig <- t(counts_batchcorrect[receptors,]) # keep only the genes of interest
# counts.correc.sig[,1:5]
pdf("./plots/MAITvsGDT_heatmap_allPBMC_cytokinereceptors.pdf", width=8, height=7)
pheatmap::pheatmap(t(counts.correc.sig),
                   color = heat_colors,
                   scale = "row", # z-score
                   # clustering_method="ward.D2",
                   cluster_rows = T,
                   cluster_cols = F,
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
                   main="All PBMCs: cytokine receptors")
dev.off()



# ************************************
## 2.6. Identifying gene clusters ####

# Tutorial: https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html

# Obtain normalized values for significant genes
genes.sig <- res %>%
  data.frame() %>%
  filter(padj<0.01)
dim(genes.sig) # 79 genes
counts.correc.sig <- counts_batchcorrect[rownames(genes.sig),]
# counts.correc.sig[1:5,]

# Use degPatterns function to show gene clusters across sample groups
# BiocManager::install("DEGreport")
library(DEGreport)
all(rownames(metadf.deseq) == colnames(counts.correc.sig)) # sanity check
metadf.deseq$lineage_id <- factor(metadf.deseq$lineage_id, levels=unique(metadf.deseq$lineage_id))
clusters <- degPatterns(counts.correc.sig,
                        metadata = metadf.deseq,
                        minc=5, # minimum nb of genes per group
                        time = "lineage_id",
                        col=NULL)
# ggsave(plot=clusters$plot, filename="./plots/degpatterns_clust13-14_padj0_01.jpeg", width=8, height=6)

# Extract gene lists
cluster_groups <- clusters$df
# write.csv(cluster_groups, "./plots/degpatterns_clust13-14_padj0_01_listgenespergroup.csv")
group1 <- clusters$df %>% filter(cluster==1)
group5 <- clusters$df %>% filter(cluster %in% 5:6)

# GENE SET ENRICHMENT ANALYSIS
# Create list of DE genes
isGeneSig <- rownames(res) %in% group5$genes
isGeneSig <- as.integer(isGeneSig)
# sum(isGeneSig)==nrow(group5)
names(isGeneSig) <- rownames(res)
isGeneSig[1:5]

# Weigh gene vector by length of genes
# BiocManager::install("goseq")
library(goseq)
pwf=nullp(isGeneSig,"hg19","geneSymbol")

# GO enrichment: based on https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#go-enrichment-analysis
goResults <- goseq(pwf, "hg19","geneSymbol", test.cats=c("GO:CC", "GO:BP", "GO:MF"))
goResults %>% 
  top_n(20, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=reorder(term, hitsPerc), 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="pvalue", size="#DE genes", title="DE genes in groups 5-6")
# ggsave("./plots/degpatterns_clust13-14_padj0_01_goenrich_groups5-6.jpeg", width=7, height=6)

