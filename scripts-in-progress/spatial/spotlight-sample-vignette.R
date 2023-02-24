library(ggplot2)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)


## spatial data
library(TENxVisiumData)
spe <- MouseKidneyCoronal()
# Use symbols instead of Ensembl IDs as feature names
rownames(spe) <- rowData(spe)$symbol

## single cell reference
library(TabulaMurisSenisData)
sce <- TabulaMurisSenisDroplet(tissues = "Kidney")$Kidney
# clean up a bit: restrict to age, keep cells with clear cell type annotations
sce <- sce[, sce$age == "18m"]
sce <- sce[, !sce$free_annotation %in% c("nan", "CD45")]

# Processing
sce <- logNormCounts(sce)
dec <- modelGeneVar(sce)
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
hvg <- getTopHVGs(dec, n = 3000)

colLabels(sce) <- colData(sce)$free_annotation

# Get vector indicating which genes are neither ribosomal or mitochondrial
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce))

# Compute marker genes
mgs <- scoreMarkers(sce, subset.row = genes)
mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.8, ]
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)


## Important: Cell Downsampling

## In our experience we have found that for this step it is better to select the cells from each cell type from the same batch if you have a joint dataset from multiple runs. This will ensure that the model removes as much signal from the batch as possible and actually learns the biological signal. ##

# split cell indices by identity
idx <- split(seq(ncol(sce)), sce$free_annotation)
# downsample to at most 20 per identity & subset
# We are using 5 here to speed up the process but set to 75-100 for your real
# life analysis
n_cells <- 5
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
})
sce <- sce[, unlist(cs_keep)]

res <- SPOTlight(
  x = sce,
  y = spe,
  groups = as.character(sce$free_annotation),
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")

mod <- res$NMF
plotTopicProfiles(
  x = mod,
  y = sce$free_annotation,
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1) +
  theme(aspect.ratio = 1)

plotTopicProfiles(
  x = mod,
  y = sce$free_annotation,
  facet = TRUE,
  min_prop = 0.01,
  ncol = 6)

mat <- res$mat
plotCorrelationMatrix(mat)
plotInteractions(mat, which = "heatmap", metric = "prop")
plotInteractions(mat, which = "network")

ct <- colnames(mat)
mat_nn <- mat
mat_nn[mat_nn < 0.1] <- 0

# Define color palette
# (here we use 'paletteMartin' from the 'colorBlindness' package)
paletteMartin <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
           "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
           "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")
           
pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct

plotSpatialScatterpie(
  x = spe,
  y = mat,
  cell_types = colnames(mat),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.4) +
  scale_fill_manual(
    values = pal,
    breaks = names(pal))
