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
  unite(sample_id, remove=FALSE) %>%
  dplyr::rename(lineage_id=cell_id)
rownames(metadf.deseq) <- metadf.deseq$sample_id

# Check that the row names of metadf.deseq are the same as the column names of counts.deseq in order to use as input to DESeq2
nrow(metadf.deseq)==ncol(counts.deseq)
metadf.deseq <- metadf.deseq[match(colnames(counts.deseq), rownames(metadf.deseq)),] # reorder rows in metadf.deseq to match order of columns in counts.deseq
all(rownames(metadf.deseq) == colnames(counts.deseq))


# PREPARE COLOR SCALES
matcol <- metadf.deseq %>%
  select(lineage_id, batch_id, cluster_id) %>%
  relocate(batch_id) %>%
  mutate(cell_state=ifelse(cluster_id %in% c(3, 6, 9:11), "Tnaive",
                           ifelse(cluster_id==12, "Tcm",
                                  ifelse(cluster_id %in% 13:14, "Th17",
                                         ifelse(cluster_id %in% 15:17, "Temra",
                                            ifelse(cluster_id==7, "Treg", "?"))))))
# tabl(matcol$cell_state)

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
dds <- dds[which(mcols(dds)$fullBetaConv),] # remove 8 genes that didn't converge

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
  # filter(padj < 1e-3) # 0.001
  arrange(padj) %>% slice_head(n=1000)
dim(genes.sig) # 3987 genes



# ***********************
## 2.4. Plot heatmap ####

# Get corrected counts (that we batch corrected with limma earlier)
counts.correc.sig <- t(counts_batchcorrect[rownames(genes.sig),]) # keep only the genes of interest
counts.correc.sig[1:5,1:5]
dim(counts.correc.sig) # 3987 genes

# Set a color palette
heat_colors <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(100))


# Run pheatmap using the metadata data frame for the annotation
phtmap <- pheatmap::pheatmap(counts.correc.sig,
         color = heat_colors,
         scale = "column", # z-score
         clustering_method="ward.D2",
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

# Get col dendrogram and flip it
col_dend <- as.dendrogram(phtmap[[2]])
col_order <- colnames(counts.correc.sig[,phtmap$tree_col[["order"]]])
col_order_clust <- dendextend::cutree(col_dend, k = 10)[order.dendrogram(col_dend)]
# unique(col_order_clust) # order: 4  7  8  3  1  5 10  9  2  6

# Reorder
bloc1 <- col_order_clust[col_order_clust %in% c(4,7, 3, 1, 5)]
bloc2 <- col_order_clust[col_order_clust ==8]
bloc3 <- col_order_clust[col_order_clust %in% c(10,9,2,6)]
new_order <- c(bloc1, bloc2, bloc3)
table(unique(col_order) %in% unique(names(new_order)), useNA="ifany") # sanity check
col_dend_new <- dendextend::rotate(col_dend, order=names(new_order))


# Re-plot pheatmap
phtmap_reordered <- pheatmap::pheatmap(counts.correc.sig,
         color = heat_colors,
         scale = "column", # z-score
         clustering_method="ward.D2",
         cluster_rows = T,
         cluster_cols = as.hclust(col_dend_new),
         # clustering_method="complete",
         border_color = NA,
         # Create gaps to separate CD4/NKT
         # gaps_row=c(2,6, 11, 14),
         # cutree_rows = 5,
         cutree_cols = 12,
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




# ***********************************************
## 2.5. Compare some DE genes with GEP genes ####

gene_order <- colnames(counts.correc.sig[,phtmap_reordered$tree_col[["order"]]])
tree_numbers <- sort(cutree(phtmap_reordered$tree_col, k=12))
unique(unname(tree_numbers[gene_order])) # ORDER OF TREE: 4  8  3  1  5 11  9 12 10  2  6  7
tabl(tree_numbers)

# ATTEMPT AT GENE LISTS
genes_naive <- names(tree_numbers[tree_numbers==7])
length(genes_naive) # 86
genes_Tcm_GEP6 <- names(tree_numbers[tree_numbers==10])
length(genes_Tcm_GEP6) # 53
genes_Th17_GEP1 <- names(tree_numbers[tree_numbers==9])
length(genes_Th17_GEP1) # 72
genes_Temra_GEP4 <- names(tree_numbers[tree_numbers==8])
length(genes_Temra_GEP4) # 93

# Put them in dataframe
DEgenes_select_df <- rbind(data.frame("cellstate"="Tnaive", "gene"=genes_naive),
                           data.frame("cellstate"="Tcm", "gene"=genes_Tcm_GEP6),
                           data.frame("cellstate"="Th17", "gene"=genes_Th17_GEP1),
                           data.frame("cellstate"="Temra", "gene"=genes_Temra_GEP4))

# Sanity checks
# pheatmap(counts.correc.sig[,c(genes_naive, genes_Tcm_GEP6, genes_Th17_GEP1, genes_Temra_GEP4)],
#     color = heat_colors,
#     scale = "column", # z-score
#     cluster_rows = T,
#     cluster_cols = T,
#     clustering_method="ward.D2",
#     border_color = NA,
#     # Change rows text format
#     show_rownames = T,
#     fontsize_row = 10,
#     annotation_row=matcol,
#     annotation_colors = list(lineage_id = cols_lineages,
#                              batch_id=cols_batchid,
#                              cluster_id=cols_integrated[levels(metadf.deseq$cluster_id)],
#                              cell_state=cols_cellstate),    # Change columns
#     angle_col=45,
#     # fontsize_col = 1,
#     show_colnames=F,
#     # title
#     main="PBMC - LRT test - batchid+cellid+clusterid vs batchid (batch corrected)")


# COMPARE WITH LAURENT's GEPs
geps.df <- read.csv("./data/human-thymus/HumanData_17_GEPsOnParkData/genes_per_GEP_df_2023-04-07.csv", row.names=1)

tabl(genes_naive %in% geps.df$GEP_5)
tabl(genes_Tcm_GEP6 %in% geps.df$GEP_6)
tabl(genes_Th17_GEP1 %in% geps.df$GEP_1)
tabl(genes_Temra_GEP4 %in% geps.df$GEP_4)

# FUNCTION FOR PLOTTING AGAINST GEPs
compare_with_geps <- function(genesvector, color, title){
  # Create DF
  df <- geps.df %>%
    pivot_longer(everything(), names_to="GEP", values_to="gene") %>%
    filter(gene %in% genesvector) %>%
    group_by(GEP) %>%
    summarize(n_gene=n()) %>%
    ungroup() %>%
    mutate(total_gene=length(genesvector),
           prop_genes=n_gene*100/total_gene)
  # Plot
  p <- ggplot(df, aes(x=factor(GEP, level=paste0("GEP_",1:12)), y=prop_genes))+
    geom_bar(stat="identity", fill=color)+
    labs(x="", y="% genes found in each GEP", title=title)+
    ylim(0,100)+
    theme_cowplot()+
    theme(axis.text.x=element_text(angle=45, hjust=1))
  return(p)
}

# First look at overlap
plot_grid(compare_with_geps(genesvector=genes_naive, color="#b3e2cd", title="Genes from Tnaive"),
          compare_with_geps(genesvector=genes_Tcm_GEP6, color="#f4cae4", title="Genes from Tcm"),
          compare_with_geps(genesvector=genes_Th17_GEP1, color="#cbd5e8", title="Genes from Th17"),
          compare_with_geps(genesvector=genes_Temra_GEP4, color="#fdcdac", title="Genes from Temra"),
          ncol=2)
# ggsave("./data/human-thymus/HumanData_19_CellIdentity/LRT_BatchCellCluster_vs_Batch/suppfigure_plots/overlappercent_withoutnull.jpeg", width=10, height=10)


# H0: any random list of 86 genes (length of Tnaive gene list) would have 55% overlap with GEP5
# H1: "Tnaive" genes are specifically overlapping with GEP5, with a 55% overlap (more than random)


# Define function that will simulate null hypothesis
NullOverlap <- function(seuratobj=seur.pbmc, df_geneprograms=longdf, cellstates=c("Tnaive", "Tcm", "Th17", "Temra"), nbins=25, nrandom=1000){
  
  # Sanity check:
  if(!cellstates %in% unique(df_geneprograms$cellstate)){print("ERROR: some cell states are not present in df_geneprograms")}
  
  # ___________________________________________________________
  # -- 1. Get DF with all genes & binned by expression level --
  cat("\n-- 1. GET ALL GENES & BIN THEM BY EXPRESSION LEVEL --\n")
  allgenes_binned_DF <- data.frame("gene"=rownames(seuratobj),
                                   "totalexpression"=rowSums(seuratobj@assays$RNA@data))
  # nrow(allgenes_binned_DF) == nrow(seuratobj) # sanity check
  # table(rownames(allgenes_binned_DF)==allgenes_binned_DF$gene) # sanity check
  # hist(allgenes_binned_DF$totalexpression, breaks=100) # visualize distribution
  
  # Get the bins by percentiles
  allgenes_binned_DF <- allgenes_binned_DF[order(allgenes_binned_DF$totalexpression),]
  allgenes_binned_DF$bin <- ggplot2::cut_number(x = allgenes_binned_DF$totalexpression, n = nbins, labels = FALSE, right = FALSE)
  cat("Number of genes in each bin, from bin 1 to bin", nbins, ":\n")
  print(table(allgenes_binned_DF$bin, useNA="ifany"))
  # ggplot(allgenes_binned_DF)+
  #   geom_density(aes(x=totalexpression, group=bin, color=factor(bin, levels=1:10))) # sanity check
  
  # ___________________________________________________________________________________________________________
  # -- 2. Draw random ctrl genes for each GEP & get observed overlap & see if significant compared to random --
  cat("\n-- 2. DRAW RANDOM CTRL GENES FOR EACH GENE LIST AND GET PVAL --\n")
  overlaps_per_cellstate_list <- list()
  
  for(state in cellstates){
    cat("\n++ Drawing random ctrl genes for", state, "\n")
    # Find in which expression bins are our genes from the gene set
    cellstategeneset <- df_geneprograms[df_geneprograms$analysis=="DESeq2" & df_geneprograms$cellstate==state, "gene"]
    cellstategeneset_binned_DF <- allgenes_binned_DF[allgenes_binned_DF$gene %in% cellstategeneset,]
    # table(nrow(cellstategeneset_binned_DF)==length(cellstategeneset))
    
    # ++++++++++++++++++++++++++++++++++
    # Do 1000 random draws of gene lists
    ctrlgenelist <- list()
    for(i in 1:nrandom){
      # Progress bar
      if(i%%100==0){cat(paste0("\nrandom draw #", i))}
      # Loop by expression bin & sample equal number of genes from each bin
      ctrlgeneset <- c()
      for(bin in unique(cellstategeneset_binned_DF$bin)){
        # print(bin)
        ngenes_to_sample <- sum(cellstategeneset_binned_DF$bin==bin)
        # print(ngenes_to_sample)
        ctrlgenes_in_same_bin <- sample(x=allgenes_binned_DF[allgenes_binned_DF$bin==bin,"gene"], size=ngenes_to_sample, replace=FALSE)
        # print(length(ctrlgenes_in_same_bin))
        ctrlgeneset <- c(ctrlgeneset, ctrlgenes_in_same_bin)
      }
      # Sanity check
      if(length(ctrlgeneset)!=length(cellstategeneset)){cat("\nPROBLEM: control gene set of size", length(ctrlgeneset), "while cellstate gene set length is", length(cellstategeneset))}
      
      # Add ctrlgeneset to list
      ctrlgenelist[[i]] <- ctrlgeneset
    }
    # Sanity check
    cat("\n\n++ Selection of random ctrl genes over! ++\n")
    cat("Length of gene set for", state,  "is", length(unique(cellstategeneset)), "and length of the random gene sets:\n")
    print(table(sapply(ctrlgenelist, length)))
    # ++++++++++++++++++++++++++++++++++
    
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Get the % of gene overlap with every GEP
    cat("\n\n++ Compute %overlap of random gene sets with GEPs\n")
    overlap_null <- data.frame()
    othergeps <- unique(df_geneprograms[df_geneprograms$analysis=="cNMF", "cellstate"])
    # Sanity check
    # cat("There are [", length(othergeps), "] GEPS (should be 12)\n")
    for(geneprog in othergeps){
      # cat("\n-> computing %random and %observed overlap with:", geneprog)
      # get genes from GEP
      gep_genes <- df_geneprograms %>% filter(cellstate==geneprog) %>% pull(gene)
      # cat(length(gep_genes)) # sanity check
      # get % gene overlap for each random draw with the gene program (% of genes from random gene set that can be found in gene program X)
      ctrlpercent <- sapply(ctrlgenelist, function(x) sum(x %in% gep_genes)*100/length(x))
      # get % gene overlap of cell state (from DEanalysis) with the GEP (% of genes from Tnaive that can be found in GEP X)
      observedpercent <- sum(cellstategeneset %in% gep_genes)*100/length(cellstategeneset)
      # print(observedpercent)
      df.temp <- data.frame("gep"=geneprog, "randomoverlap"=ctrlpercent, "observedoverlap"=observedpercent)
      # add rows
      overlap_null <- rbind(overlap_null, df.temp)
    }
    # table(overlap_null$gep) # should all be 1000 (for the 1000 draws)
    # table(overlap_null$observedoverlap)
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # Summarize overlap
    cat("\n\n++ Summarize %overlap of random/cellstate gene sets with GEPs\n")
    final_df <- overlap_null %>%
      as_tibble() %>%
      group_by(gep) %>%
      summarize(observedoverlap=unique(observedoverlap),
                mean_randomoverlap=mean(randomoverlap),
                min_randomoverlap=min(randomoverlap),
                max_randomoverlap=max(randomoverlap),
                nb_random_above_observed=sum(randomoverlap>observedoverlap)) %>%
      mutate(pval=nb_random_above_observed/nrandom,
             padj=pval*length(unique(othergeps)),
             genesetsize=length(unique(cellstategeneset))) %>%
      arrange(-observedoverlap)
    
    # sanity check
    cat("GEP with highest overlap with", state, "is [", pull(final_df[1,"gep"]), "] with [", pull(final_df[1,"observedoverlap"]), "% ] overlap\n")
    
    # Add to list
    overlaps_per_cellstate_list[[state]] <- final_df
    cat("\n ----------------------------------------------- \n")
  } # end of GEP loop
  
  # Collapse list into one big dataframe <3
  bigdf <- bind_rows(overlaps_per_cellstate_list, .id="cellstate")
  
  return(bigdf)
}


# GET OVERLAP (RANDOM & OBSERVED)

# Get long DF with list of all genes
DEgenes_select_df$analysis <- "DESeq2"
head(DEgenes_select_df)
GEPgenes_df <- geps.df %>%
  pivot_longer(cols=colnames(geps.df), names_to="cellstate", values_to="gene") %>%
  drop_na(gene) %>%
  mutate(analysis="cNMF")
longdf <- rbind(DEgenes_select_df,GEPgenes_df)
# tabl(longdf$analysis)
# tabl(longdf$cellstate)

# subset seurat object to PBMCs only
seur.pbmc <- subset(seur.human, Tissue=="PBMC")
print(seur.pbmc) # 17,204 genes and 41,238 cells
# tabl(seur.pbmc$Tissue) # sanity check


# df_geneprograms should contain 3 columns: analysis (DESeq2, cNMF); cellstate (Tnaive, Temra, ..., GEP1, GEP2, etc.); gene
genesets_overlap <- NullOverlap(seuratobj=seur.pbmc,
                                df_geneprograms = longdf,
                                cellstates=unique(longdf[longdf$analysis=="DESeq2","cellstate"]),
                                nbins=25, nrandom=1000)
# Sanity check
# tabl(genesets_overlap$cellstate) # should all be 12 (for the 12 GEPs)


# Plot with facet
ggplot(genesets_overlap %>%
         # left_join(cols_df, by="geneprogram") %>%
         mutate(padj_toplot=ifelse(padj==0, "< 0.001", ifelse(padj<0.05, as.character(round(padj, 2)), ""))),
       aes(x=factor(gep, levels=paste0("GEP_", 1:12)), y=observedoverlap))+
  # geom_hline(yintercept = 10, color="lightgrey", linetype=2)+
  geom_bar(stat="identity", aes(fill=cellstate))+
  facet_wrap(~factor(cellstate, levels=c("Temra", "Th17", "Tcm", "Tnaive")), ncol=1, scales="free")+
  geom_point(aes(y=mean_randomoverlap), shape="_", size=3)+
  # tidytext::scale_x_reordered() +
  # ggrepel::geom_text_repel(aes(label=padj_toplot), size=4, direction="y", nudge_y=1, force=4) +
  geom_text(aes(label=padj_toplot), size=3, angle=90, hjust=-0.2)+
  ylim(c(0,100))+
  scale_fill_manual(values=cols_cellstate, name="")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10), #legend.position="none",
        strip.text.x = element_text(size = 20),
        panel.grid.major.y=element_line(colour="lightgrey", linetype=2),
        legend.position="none"
        #plot.margin=margin(10,10,10,70)
        )+
  labs(x="", y="% genes found in each GEP")
ggsave("./data/human-thymus/HumanData_19_CellIdentity/LRT_BatchCellCluster_vs_Batch/suppfigure_plots/overlappercent_withnull.pdf", width=5, height=15)

# Export heatmap
pdf("./data/human-thymus/HumanData_19_CellIdentity/LRT_BatchCellCluster_vs_Batch/suppfigure_plots/heatmap_top1000genes_with_lowest_padj.pdf", width=10, height=10)
phtmap_reordered
dev.off()

# Add nb of cells to heatmap
row_order <- rownames(counts.correc.sig[phtmap$tree_row[["order"]],])
ncells <- groups %>%
  mutate(label=paste(cell_id, batch_id, cluster_id, sep="_")) %>%
  group_by(label) %>%
  count() %>%
  ungroup()
nrow(ncells)==nrow(counts.correc.sig) # should be TRUE

ggplot(data=ncells,
       aes(x=factor(label, levels=rev(row_order)), y=n))+
  geom_bar(stat="identity", fill="#bdbdbd") +
  scale_x_discrete(position="top") +
  # scale_y_continuous(limits=c(0,100), breaks=c(0,50,100))+
  labs(y="# cells")+ coord_flip() + theme_cowplot()+
  theme(axis.text    = element_text(size=15),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size=20),
        axis.title.y = element_blank(),
        axis.line  = element_blank(),
        panel.border = element_rect(fill=NULL, colour="black", linewidth = 0.5),
        legend.position = "none")
ggsave("./data/human-thymus/HumanData_19_CellIdentity/LRT_BatchCellCluster_vs_Batch/suppfigure_plots/heatmap_top1000genes_with_lowest_padj_nbcells.pdf", height=15, width=3)
