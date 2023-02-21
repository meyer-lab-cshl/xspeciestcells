library(SpatialExperiment)
library(ggspavis)
library(scater)
library(scran)

# https://lmweber.org/OSTA-book/quality-control.html
directory <- "data/GSE207205_RAW"
sample_ids <- c("GSM6281327", "GSM6281320", "GSM6281322",  "GSM6281321",
                "GSM6281323", "GSM6281324", "GSM6281325", "GSM6281326"
                )

samples <- file.path(directory, sample_ids, "outs")

spe <- read10xVisium(samples, sample_ids, type = "sparse", data = "raw", 
                     images = "lowres", load = FALSE)

plotSpots(spe)

# exclude background spots
spe <- spe[, colData(spe)$in_tissue == 1]
plotSpots(spe)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)

# calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
head(colData(spe), 3)

# histograms of QC metrics
par(mfrow = c(1, 4))
hist(colData(spe)$sum, xlab = "sum", main = "UMIs per spot")
hist(colData(spe)$detected, xlab = "detected", main = "Genes per spot")
hist(colData(spe)$subsets_mito_percent, xlab = "percent mitochondrial", main = "Percent mito UMIs")
#hist(colData(spe)$cell_count, xlab = "number of cells", main = "No. cells per spot")

par(mfrow = c(1, 1))

# select QC thresholds
qc_lib_size <- colData(spe)$sum < 600
qc_detected <- colData(spe)$detected < 400
qc_mito <- colData(spe)$subsets_mito_percent > 28
#qc_cell_count <- colData(spe)$cell_count > 10

# number of discarded spots for each metric
apply(cbind(qc_lib_size, qc_detected, qc_mito), 2, sum)

# combined set of discarded spots
discard <- qc_lib_size | qc_detected | qc_mito 
table(discard)

colData(spe)$discard <- discard

# check spatial pattern of discarded spots
plotQC(spe, type = "spots", discard = "discard")

# filter low-quality spots
spe <- spe[, !colData(spe)$discard]
dim(spe)

# calculate library size factors
spe <- computeLibraryFactors(spe)

summary(sizeFactors(spe))
hist(sizeFactors(spe), breaks = 20)

# calculate logcounts and store in object
spe <- logNormCounts(spe)
assayNames(spe)

# fit mean-variance relationship
dec <- modelGeneVar(spe)

# visualize mean-variance relationship
fit <- metadata(dec)
plot(fit$mean, fit$var, 
     xlab = "mean of log-expression", ylab = "variance of log-expression")
curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

# select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)


# compute PCA
set.seed(123)
spe <- runPCA(spe, subset_row = top_hvgs)

# compute UMAP on top 50 PCs
spe <- runUMAP(spe, dimred = "PCA")

# update column names for easier plotting
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# graph-based clustering
set.seed(123)
k <- 50
g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)

# store cluster labels in column 'label' in colData
colLabels(spe) <- factor(clus)

#colors <-  c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd')

colors <- c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5')

# plot clusters in spatial x-y coordinates
plotSpots(spe, annotate = "label",palette =colors)

plotDimRed(spe, type = "UMAP", 
           annotate = "label", palette =colors)

# set gene names as row names for easier plotting
rownames(spe) <- rowData(spe)$gene_name

# test for marker genes
markers <- findMarkers(spe, test = "binom", direction = "up")

# returns a list with one DataFrame per cluster
markers


