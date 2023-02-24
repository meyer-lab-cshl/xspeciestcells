###
# Purpose: Import & plot human iNKT data with SCTransform and Harmony (curtesy of Laurent Gapin)
# Date: Oct 4th 2022
# Author: Laurent Gapin
###


#### IMPORT ####

# Import librairies
library(Seurat)
library(tidyverse)
library(SeuratWrappers)
library(harmony)
library(ggplot2)
library(MetBrewer)
library(viridis)
library(scales)
library(RColorBrewer)
library(cowplot)
library(ggfortify)
library(glue)
library(VISION)


# Import data
# Laurent's paths
path1 <- "/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY11/CUTHY11BDRscRNA_seq_091621_SampleTag08_hs_NKT_RSEC_MolsPerCell.csv"
path2 <- "/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY12/CUTHY12BDRscRNA_seq_211101_SampleTag12_hs_NKT_RSEC_MolsPerCell.csv"
path3 <- "/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY13/CUThy13_220225_SampleTag05_hs_NKT/CUThy13_220225_SampleTag05_hs_NKT_RSEC_MolsPerCell.csv"
# Salome's paths
path1 <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/human_data/CUTHY11BDRscRNA_seq_091621_SampleTag08_hs_NKT_RSEC_MolsPerCell.csv"
path2 <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/human_data/CUTHY12BDRscRNA_seq_211101_SampleTag12_hs_NKT_RSEC_MolsPerCell.csv"
path3 <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/human_data/CUThy13_220225_SampleTag05_hs_NKT_RSEC_MolsPerCell.csv"
# Import data
all.data <- list()
all.data[['NKT_1_Thymus']] <- read.csv(path1,
                                       header = TRUE, skip = 8)
all.data[['NKT_2_Thymus']] <- read.csv(path2,
                                       header = TRUE, skip = 8)
all.data[['NKT_3_Thymus']] <- read.csv(path3,
                                       header = TRUE, skip = 7)




#### PRE-PROCESSING ####

# Get sample names and gene names
shared_genes <- Reduce(intersect, lapply(X = all.data, FUN = function(x){ colnames(x) }))
shared_genes <- shared_genes[!shared_genes %in% 'Cell_Index']
sample.names <- names(all.data)


# Create seurat object
all.seurat <- lapply(X = 1:length(all.data), FUN = function(x){
  sample.counts <- all.data[[x]]
  rownames(sample.counts) <- sample.counts$Cell_Index
  sample.counts$Cell_Index <- NULL
  transpose.df <- as.data.frame(t(as.matrix(sample.counts)))
  transpose.df <- transpose.df[match(shared_genes, rownames(transpose.df)), ]
  seurat.obj <- CreateSeuratObject(counts = transpose.df, project = sample.names[x])
})
merged_seurat <- merge(x = all.seurat[[1]], y = all.seurat[2:length(all.seurat)],
                       add.cell.ids = sample.names, project = 'Human_Innate') # 2661 cells with 19,693 genes


# QC
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT\\.")

VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

filtered_seurat <- subset(merged_seurat, subset = nFeature_RNA >= 500 & nFeature_RNA < 4000 & nCount_RNA >=500 & percent.mt < 25) # 2558 cells
VlnPlot(filtered_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Keep only non-ribosomal genes and genes that have a count of at least 10
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
genes.to.remove <- rownames(filtered_seurat)[grep(rownames(filtered_seurat), pattern = "^RPL|^RPS|^MT\\.")]
keep_genes[which(names(keep_genes) %in% genes.to.remove)] = FALSE
filtered_counts <- counts[keep_genes, ]

filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data) # 11,304 genes


# Import metadata
sample_metadata <- read.csv(file = "/Volumes/Samsung_T5/Human_MAIT_NKT/Data/Metadata.csv") %>% 
  dplyr::select(-sample.id, -X, -file)

rownames(sample_metadata) <- sample_metadata$orig.ident
cell_metadata <- sample_metadata[filtered_seurat$orig.ident,]
for(meta in names(cell_metadata)){
  filtered_seurat[[meta]] <- cell_metadata[[meta]]
}

filtered_seurat@meta.data


# Normalize with SCTransform
filtered_seurat <- NormalizeData(filtered_seurat)
filtered_seurat <- SCTransform(filtered_seurat, vars.to.regress = c('percent.mt'), 
                               verbose=TRUE, method = "glmGamPoi")

# Run dimension reduction (linear & non-linear)
filtered_seurat <- RunPCA(filtered_seurat)

ElbowPlot(filtered_seurat)

filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:15, min.dist = 0.3, spread = 0.5)
filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:15)
filtered_seurat <- FindClusters(filtered_seurat, resolution = 0.7)


# Plot UMAP
my_cols <- c(met.brewer("Austria", n = 7, type = "discrete"),
             met.brewer("Egypt", n = 4, type = "discrete"),
             met.brewer("Cross", n = 9, type = "discrete"),
             met.brewer("Troy", n = 2, type = "discrete"),
             met.brewer("Derain", n = 7, type = "discrete"),
             met.brewer("Hiroshige", n = 10, type = "discrete")
)

DimPlot(filtered_seurat, group.by = c("orig.ident"),
        pt.size = 1, 
        ncol = 1, 
        label = TRUE, label.size = 6, 
        cols = alpha(my_cols, 0.7)) + NoLegend()


# Run normalization with fast MNN
DefaultAssay(filtered_seurat) <- "RNA"
filtered_seurat_MNN <- RunFastMNN(object.list = SplitObject(filtered_seurat, split.by = "orig.ident"), k = 20, 
                                  auto.merge = TRUE)
filtered_seurat_MNN <- RunUMAP(filtered_seurat_MNN, reduction = "mnn", dims = 1:15)
filtered_seurat_MNN <- FindNeighbors(filtered_seurat_MNN, reduction = "mnn", dims = 1:15)
filtered_seurat_MNN <- FindClusters(filtered_seurat_MNN, resolution = 0.5)

p1 <- DimPlot(filtered_seurat_MNN, group.by = c("seurat_clusters"),
              pt.size = 1, 
              ncol = 1, 
              label = TRUE, label.size = 6, 
              cols = alpha(my_cols, 0.7)) + NoLegend()

DefaultAssay(filtered_seurat_MNN) <- "RNA"
FeaturePlot(filtered_seurat_MNN, features = c("ZBTB16", "EOMES", "GZMK", "CCR7", "SELL", "PLAC8"),
            order = T,
            pt.size = 0.7, 
            ncol = 2, 
            label = FALSE, repel = FALSE, label.size = 6) + NoLegend()


# Run integration with Harmony (from the SCTransform data)
set.seed(0229)

DefaultAssay(filtered_seurat) <- "SCT"
# filtered_seurat_Harmony <- RunHarmony(filtered_seurat, group.by.vars = c("Batch", "Method"), assay.use = "SCT",
filtered_seurat_Harmony <- RunHarmony(filtered_seurat, group.by.vars = c("orig.ident"), assay.use = "SCT",
                                      max.iter.harmony = 20) # 11 iterations
filtered_seurat_Harmony <- RunUMAP(filtered_seurat_Harmony, reduction = "harmony", 
                                   dims = 1:15)
filtered_seurat_Harmony <- FindNeighbors(filtered_seurat_Harmony, 
                                         reduction = "harmony", dims = 1:15)
filtered_seurat_Harmony <- FindClusters(filtered_seurat_Harmony, resolution = 0.5)

p2 <- DimPlot(filtered_seurat_Harmony, group.by = c("seurat_clusters"),
              pt.size = 1, 
              ncol = 1, 
              label = TRUE, label.size = 6, 
              cols = alpha(my_cols, 0.7)) + NoLegend()

DefaultAssay(filtered_seurat_Harmony) <- "RNA"
FeaturePlot(filtered_seurat_Harmony, features = c("ZBTB16", "EOMES", "SELL", "CCR7", "CD4", "ISG15"),
            order = T,
            pt.size = 0.5, 
            ncol = 2, 
            label = FALSE, repel = FALSE, label.size = 6) + NoLegend()

p1 <- p1 + theme(legend.position="top", title = element_text(size=10)) + labs(title="FastMNN")
p2 <- p2 + theme(legend.position="top", title = element_text(size=10)) + labs(title="Harmony")

cluster4 <- subset(filtered_seurat_MNN, idents = 4)
cells_cluster4 <- WhichCells(cluster4)

p3 <- DimPlot(filtered_seurat_Harmony, cells.highlight = cells_cluster4,
              pt.size = 1, 
              ncol = 1, 
              label = TRUE, label.size = 6) + NoLegend()

p4 <- FeaturePlot(filtered_seurat_Harmony, features = "ZBTB16",
                  pt.size = 0.7, 
                  ncol = 1, order = T) + NoLegend()

fig_dir <- "/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/"

pdf(paste0(fig_dir, "Umap_comparison.pdf"), width=10, height=10)
ggdraw()+
  draw_plot(p1, x=-0, y=.50, width=0.45, height=0.5)+
  draw_plot(p2, x=0.5, y=.50, width=0.45, height=0.5)+
  draw_plot(p3, x=0, y=-0, width=0.45, height=0.45)+
  draw_plot(p4, x=0.5, y=0, width=0.45, height=0.45)
dev.off()

###
signatures <- c("/Volumes/Samsung_T5/Human_MAIT_NKT/c5.go.v7.4.symbols.gmt", 
                "/Volumes/Samsung_T5/Human_MAIT_NKT/MouseNKT.signatures.gmt")

Idents(filtered_seurat_MNN) <- "RNA"
vision.obj.NKT.MNN <- Vision(filtered_seurat_MNN@assays$RNA@data, signatures = signatures, projection_methods = NULL, 
                             meta = filtered_seurat_MNN@meta.data)

vision.obj.NKT.MNN <- analyze(vision.obj.NKT.MNN)

viewResults(vision.obj.NKT.MNN)

# Display autocorrelation coefficients, p-values for signatures
head(getSignatureAutocorrelation(vision.obj.NKT.MNN))

correlations <- getSignatureAutocorrelation(vision.obj.NKT.MNN)
correlations <- rownames_to_column(correlations, var = "pathway")
correlations %>% dplyr::filter(str_detect(pathway, "NKT"))

metadata_correlation <- getMetaAutocorrelation(vision.obj.NKT.MNN)
metadata_correlation <- rownames_to_column(metadata_correlation, var = "Metadata")

umap <- filtered_seurat_MNN[["umap"]]@cell.embeddings

corr_response_Stg0 <- correlations %>% dplyr::filter(str_detect(pathway, "Stage0_signature"))
sigScores_Stg0 <- getSignatureScores(vision.obj.NKT.MNN)[, "Stage0_signature"]
p5 <- ggplot() + aes(x=umap[, 1], y=umap[, 2], color=sigScores_Stg0) + 
  geom_point(size = 1) +
  xlim(-9, 9) + ylim(-4, 4) +
  scale_colour_viridis(option = "magma", limits = c(0, 1), oob = scales::squish) +
  ggtitle(glue("Stage0_signature\n C' = {round(corr_response_Stg0[2], digit = 2)}, FDR = {round(corr_response_Stg0[4], digit = 2)}")) +
  labs(x = "", y = "") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 0.2)) + NoLegend()

corr_response_NKTp <- correlations %>% dplyr::filter(str_detect(pathway, "NKTp_signature"))
sigScores_NKTp <- getSignatureScores(vision.obj.NKT.MNN)[, "NKTp_signature"]
p6 <- ggplot() + aes(x=umap[, 1], y=umap[, 2], color=sigScores_NKTp) + 
  geom_point(size = 1) +
  xlim(-9, 9) + ylim(-4, 4) +
  scale_colour_viridis(option = "magma", limits = c(0, 1), oob = scales::squish) +
  ggtitle(glue("NKTp_signature\n C' = {round(corr_response_NKTp[2], digit = 2)}, FDR = {round(corr_response_NKTp[4], digit = 2)}")) +
  labs(x = "", y = "") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 0.2)) + NoLegend()

corr_response_NKT2 <- correlations %>% dplyr::filter(str_detect(pathway, "NKT2_signature"))
sigScores_NKT2 <- getSignatureScores(vision.obj.NKT.MNN)[, "NKT2_signature"]
p7 <- ggplot() + aes(x=umap[, 1], y=umap[, 2], color=sigScores_NKT2) + 
  geom_point(size = 1) +
  xlim(-9, 9) + ylim(-4, 4) +
  scale_colour_viridis(option = "magma", limits = c(0, 1), oob = scales::squish) +
  ggtitle(glue("NKT2_signature\n C' = {round(corr_response_NKT2[2], digit = 2)}, FDR = {round(corr_response_NKT2[4], digit = 2)}")) +
  labs(x = "", y = "") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 0.2)) + NoLegend()

corr_response_NKT1 <- correlations %>% dplyr::filter(str_detect(pathway, "NKT1_signature"))
sigScores_NKT1 <- getSignatureScores(vision.obj.NKT.MNN)[, "NKT1_signature"]
p8 <- ggplot() + aes(x=umap[, 1], y=umap[, 2], color=sigScores_NKT1) + 
  geom_point(size = 1) +
  xlim(-9, 9) + ylim(-4, 4) +
  scale_colour_viridis(option = "magma", limits = c(0, 1), oob = scales::squish) +
  ggtitle(glue("NKT1_signature\n C' = {round(corr_response_NKT1[2], digit = 2)}, FDR = {round(corr_response_NKT1[4], digit = 2)}")) +
  labs(x = "", y = "") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 0.2)) + NoLegend()

corr_response_NKT17 <- correlations %>% dplyr::filter(str_detect(pathway, "NKT17_signature"))
sigScores_NKT17 <- getSignatureScores(vision.obj.NKT.MNN)[, "NKT17_signature"]
p9 <- ggplot() + aes(x=umap[, 1], y=umap[, 2], color=sigScores_NKT17) + 
  geom_point(size = 1) +
  xlim(-9, 9) + ylim(-4, 4) +
  scale_colour_viridis(option = "magma", limits = c(0, 1), oob = scales::squish) +
  ggtitle(glue("NKT17_signature\n C' = {round(corr_response_NKT17[2], digit = 2)}, FDR = {round(corr_response_NKT17[4], digit = 2)}")) +
  labs(x = "", y = "") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 0.2)) + NoLegend()

pdf(paste0(fig_dir, "Vision_Signatures.pdf"), width=15, height=20)
ggdraw()+
  draw_plot(p5, x=-0, y=0.5, width=0.3, height=0.3)+
  draw_plot(p6, x=0.3, y=0.5, width=0.3, height=0.3)+
  draw_plot(p7, x=0.6, y=0.5, width=0.3, height=0.3)+
  draw_plot(p8, x=0.15, y=0.15, width=0.3, height=0.3)+
  draw_plot(p9, x=0.50, y=0.15, width=0.3, height=0.3)
dev.off()