library(Seurat)
library(tidyverse)
library(SeuratWrappers)
library(harmony)
library(ggplot2)
library(MetBrewer)
library(ggalluvial)
library(viridis)
library(Matrix)
library(data.table)
library(scales)
library(RColorBrewer)
library(cowplot)
library(ggfortify)
library(ggridges)
library(reticulate)
library(ggseqlogo)
library(glue)
library(vegan)
library(circlize)
library(lisi)
library(VISION)
library(dyno)
library(monocle3)
library(slingshot)
library(tradeSeq)
library(SCpubr)
library(clustree)
library(DR.SC)
library(supCPM)
library(SCPA)
library(pheatmap)

####START of the analysis
all.data <- list()
all.data[['CD4_1_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY08/CUTHY08BDRscRNA_seq_051321_SampleTag01_hs_CD4SP/CUTHY08BDRscRNA_seq_051321_SampleTag01_hs_CD4SP_RSEC_MolsPerCell.csv",
                                       header = TRUE, skip = 8)
all.data[['CD8_1_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY08/CUTHY08BDRscRNA_seq_051321_SampleTag02_hs_CD8SP/CUTHY08BDRscRNA_seq_051321_SampleTag02_hs_CD8SP_RSEC_MolsPerCell.csv",
                                       header = TRUE, skip = 8)
all.data[['CD1a_1_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY08/CUTHY08BDRscRNA_seq_051321_SampleTag04_hs_CD1a/CUTHY08BDRscRNA_seq_051321_SampleTag04_hs_CD1a_RSEC_MolsPerCell.csv",
                                        header = TRUE, skip = 8)
all.data[['CD1c_1_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY08/CUTHY08BDRscRNA_seq_051321_SampleTag03_hs_CD1c/CUTHY08BDRscRNA_seq_051321_SampleTag03_hs_CD1c_RSEC_MolsPerCell.csv",
                                        header = TRUE, skip = 8)
all.data[['CD4_2_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY11/CUTHY11BDRscRNA_seq_091621_SampleTag05_hs_CD4_RSEC_MolsPerCell.csv",
                                       header = TRUE, skip = 8)
all.data[['CD8_2_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY11/CUTHY11BDRscRNA_seq_091621_SampleTag06_hs_CD8_RSEC_MolsPerCell.csv",
                                       header = TRUE, skip = 8)
all.data[['MAIT_1_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY11/CUTHY11BDRscRNA_seq_091621_SampleTag07_hs_MAIT_RSEC_MolsPerCell.csv",
                                        header = TRUE, skip = 8)
all.data[['NKT_1_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY11/CUTHY11BDRscRNA_seq_091621_SampleTag08_hs_NKT_RSEC_MolsPerCell.csv",
                                       header = TRUE, skip = 8)
all.data[['CD4_3_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY12/CUTHY12BDRscRNA_seq_211101_SampleTag09_hs_CD4_RSEC_MolsPerCell.csv",
                                       header = TRUE, skip = 8)
all.data[['CD8_3_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY12/CUTHY12BDRscRNA_seq_211101_SampleTag10_hs_CD8_RSEC_MolsPerCell.csv",
                                       header = TRUE, skip = 8)
all.data[['MAIT_2_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY12/CUTHY12BDRscRNA_seq_211101_SampleTag11_hs_MAIT_RSEC_MolsPerCell.csv",
                                        header = TRUE, skip = 8)
all.data[['NKT_2_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY12/CUTHY12BDRscRNA_seq_211101_SampleTag12_hs_NKT_RSEC_MolsPerCell.csv",
                                       header = TRUE, skip = 8)
all.data[['GD_1_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY12/CUTHY12BDRscRNA_seq_211101_SampleTag01_hs_GD_RSEC_MolsPerCell.csv",
                                      header = TRUE, skip = 8)
all.data[['CD4_1_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/HC_55/HC55_BDRscRNA_seq_211129_SampleTag02_hs_CD4_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 8)
all.data[['CD8_1_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/HC_55/HC55_BDRscRNA_seq_211129_SampleTag03_hs_CD8_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 8)
all.data[['MAIT_1_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/HC_55/HC55_BDRscRNA_seq_211129_SampleTag04_hs_MAIT_RSEC_MolsPerCell.csv",
                                      header = TRUE, skip = 8)
all.data[['NKT_1_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/HC_55/HC55_BDRscRNA_seq_211129_SampleTag05_hs_iNKT_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 8)
all.data[['GD_1_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/HC_55/HC55_BDRscRNA_seq_211129_SampleTag06_hs_GD_RSEC_MolsPerCell.csv",
                                    header = TRUE, skip = 8)
all.data[['NKT_3_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY13/CUThy13_220225_SampleTag05_hs_NKT/CUThy13_220225_SampleTag05_hs_NKT_RSEC_MolsPerCell.csv",
                                       header = TRUE, skip = 7)
all.data[['CD4_4_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY13/CUThy13_220225_SampleTag01_hs_CD4/CUThy13_220225_SampleTag01_hs_CD4_RSEC_MolsPerCell.csv",
                                       header = TRUE, skip = 7)
all.data[['CD8_4_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY13/CUThy13_220225_SampleTag02_hs_CD8/CUThy13_220225_SampleTag02_hs_CD8_RSEC_MolsPerCell.csv",
                                       header = TRUE, skip = 7)
all.data[['GD_2_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY13/CUThy13_220225_SampleTag03_hs_GD/CUThy13_220225_SampleTag03_hs_GD_RSEC_MolsPerCell.csv",
                                      header = TRUE, skip = 7)
all.data[['MAIT_3_Thymus']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/CUTHY13/CUThy13_220225_SampleTag04_hs_MAIT/CUThy13_220225_SampleTag04_hs_MAIT_RSEC_MolsPerCell.csv",
                                        header = TRUE, skip = 7)
all.data[['NKT_2_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06/PBMC_LRS_05_06_220225_SampleTag08_hs_NKT41/PBMC_LRS_05_06_220225_SampleTag08_hs_NKT41_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)
all.data[['MAIT_2_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06/PBMC_LRS_05_06_220225_SampleTag09_hs_MAIT41/PBMC_LRS_05_06_220225_SampleTag09_hs_MAIT41_RSEC_MolsPerCell.csv",
                                      header = TRUE, skip = 7)
all.data[['GD_2_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06/PBMC_LRS_05_06_220225_SampleTag07_hs_GD41/PBMC_LRS_05_06_220225_SampleTag07_hs_GD41_RSEC_MolsPerCell.csv",
                                    header = TRUE, skip = 7)
all.data[['NKT_3_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06/PBMC_LRS_05_06_220225_SampleTag11_hs_NKT68/PBMC_LRS_05_06_220225_SampleTag11_hs_NKT68_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)
all.data[['MAIT_3_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06/PBMC_LRS_05_06_220225_SampleTag12_hs_MAIT68/PBMC_LRS_05_06_220225_SampleTag12_hs_MAIT68_RSEC_MolsPerCell.csv",
                                      header = TRUE, skip = 7)
all.data[['CD8_2_PBMC']] <- read.csv("/Volumes/Samsung_T5/KIR3DL3/KIR3DL3-02-18-22_SampleTag04_hs/KIR3DL3-02-18-22_SampleTag04_hs_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)
all.data[['GD_3_PBMC']] <- read.csv("/Volumes/Samsung_T5/KIR3DL3/KIR3DL3-02-18-22_SampleTag07_hs/KIR3DL3-02-18-22_SampleTag07_hs_RSEC_MolsPerCell.csv",
                                    header = TRUE, skip = 7)
all.data[['CD8_3_PBMC']] <- read.csv("/Volumes/Samsung_T5/KIR3DL3/KIR3DL3-02-18-22_SampleTag05_hs/KIR3DL3-02-18-22_SampleTag05_hs_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)
all.data[['GD_4_PBMC']] <- read.csv("/Volumes/Samsung_T5/KIR3DL3/KIR3DL3-02-18-22_SampleTag08_hs/KIR3DL3-02-18-22_SampleTag08_hs_RSEC_MolsPerCell.csv",
                                    header = TRUE, skip = 7)
all.data[['CD8_4_PBMC']] <- read.csv("/Volumes/Samsung_T5/KIR3DL3/KIR3DL3-02-18-22_SampleTag06_hs/KIR3DL3-02-18-22_SampleTag06_hs_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)
all.data[['GD_5_PBMC']] <- read.csv("/Volumes/Samsung_T5/KIR3DL3/KIR3DL3-02-18-22_SampleTag09_hs/KIR3DL3-02-18-22_SampleTag09_hs_RSEC_MolsPerCell.csv",
                                    header = TRUE, skip = 7)
all.data[['CD8_5_PBMC']] <- read.csv("/Volumes/Samsung_T5/KIR3DL3/3DL3-LRS08-09-11_SampleTag01_hs/3DL3-LRS08-09-11_SampleTag01_hs_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)
all.data[['GD_6_PBMC']] <- read.csv("/Volumes/Samsung_T5/KIR3DL3/3DL3-LRS08-09-11_SampleTag04_hs/3DL3-LRS08-09-11_SampleTag04_hs_RSEC_MolsPerCell.csv",
                                    header = TRUE, skip = 7)
all.data[['CD8_6_PBMC']] <- read.csv("/Volumes/Samsung_T5/KIR3DL3/3DL3-LRS08-09-11_SampleTag02_hs/3DL3-LRS08-09-11_SampleTag02_hs_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)
all.data[['GD_7_PBMC']] <- read.csv("/Volumes/Samsung_T5/KIR3DL3/3DL3-LRS08-09-11_SampleTag05_hs/3DL3-LRS08-09-11_SampleTag05_hs_RSEC_MolsPerCell.csv",
                                    header = TRUE, skip = 7)
all.data[['CD8_7_PBMC']] <- read.csv("/Volumes/Samsung_T5/KIR3DL3/3DL3-LRS08-09-11_SampleTag03_hs/3DL3-LRS08-09-11_SampleTag03_hs_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)
all.data[['GD_8_PBMC']] <- read.csv("/Volumes/Samsung_T5/KIR3DL3/3DL3-LRS08-09-11_SampleTag06_hs/3DL3-LRS08-09-11_SampleTag06_hs_RSEC_MolsPerCell.csv",
                                    header = TRUE, skip = 7)
all.data[['CD4_2_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06_2/LRS_05_06_BDscRNAseq_220506_SampleTag01_hs_CD4_LRS08/LRS_05_06_BDscRNAseq_220506_SampleTag01_hs_CD4_LRS08_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)
all.data[['CD8_8_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06_2/LRS_05_06_BDscRNAseq_220506_SampleTag02_hs_CD8_LRS08/LRS_05_06_BDscRNAseq_220506_SampleTag02_hs_CD8_LRS08_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)
all.data[['GD_9_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06_2/LRS_05_06_BDscRNAseq_220506_SampleTag03_hs_GD_LRS08/LRS_05_06_BDscRNAseq_220506_SampleTag03_hs_GD_LRS08_RSEC_MolsPerCell.csv",
                                    header = TRUE, skip = 7)
all.data[['MAIT_4_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06_2/LRS_05_06_BDscRNAseq_220506_SampleTag04_hs_MAIT_LRS08/LRS_05_06_BDscRNAseq_220506_SampleTag04_hs_MAIT_LRS08_RSEC_MolsPerCell.csv",
                                      header = TRUE, skip = 7)
all.data[['NKT_4_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06_2/LRS_05_06_BDscRNAseq_220506_SampleTag05_hs_NKT_LRS08/LRS_05_06_BDscRNAseq_220506_SampleTag05_hs_NKT_LRS08_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)
all.data[['CD4_3_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06_2/LRS_05_06_BDscRNAseq_220506_SampleTag06_hs_CD4_LRS05/LRS_05_06_BDscRNAseq_220506_SampleTag06_hs_CD4_LRS05_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)
all.data[['CD8_9_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06_2/LRS_05_06_BDscRNAseq_220506_SampleTag07_hs_CD8_LRS05/LRS_05_06_BDscRNAseq_220506_SampleTag07_hs_CD8_LRS05_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)
all.data[['GD_10_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06_2/LRS_05_06_BDscRNAseq_220506_SampleTag08_hs_GD_LRS05/LRS_05_06_BDscRNAseq_220506_SampleTag08_hs_GD_LRS05_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)
all.data[['CD4_4_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06_2/LRS_05_06_BDscRNAseq_220506_SampleTag09_hs_CD4_LRS06/LRS_05_06_BDscRNAseq_220506_SampleTag09_hs_CD4_LRS06_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)
all.data[['CD8_10_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06_2/LRS_05_06_BDscRNAseq_220506_SampleTag10_hs_CD8_LRS06/LRS_05_06_BDscRNAseq_220506_SampleTag10_hs_CD8_LRS06_RSEC_MolsPerCell.csv",
                                      header = TRUE, skip = 7)
all.data[['GD_11_PBMC']] <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/PBMC_LRS_05_06_2/LRS_05_06_BDscRNAseq_220506_SampleTag11_hs_GD_LRS06/LRS_05_06_BDscRNAseq_220506_SampleTag11_hs_GD_LRS06_RSEC_MolsPerCell.csv",
                                     header = TRUE, skip = 7)

shared_genes <- Reduce(intersect, lapply(X = all.data, FUN = function(x){ colnames(x) }))
shared_genes <- shared_genes[!shared_genes %in% 'Cell_Index']
sample.names <- names(all.data)

all.seurat <- lapply(X = 1:length(all.data), FUN = function(x){
  sample.counts <- all.data[[x]]
  rownames(sample.counts) <- sample.counts$Cell_Index
  sample.counts$Cell_Index <- NULL
  transpose.df <- as.data.frame(t(as.matrix(sample.counts)))
  transpose.df <- transpose.df[match(shared_genes, rownames(transpose.df)), ]
  seurat.obj <- CreateSeuratObject(counts = transpose.df, project = sample.names[x])
})

merged_seurat <- merge(x = all.seurat[[1]], y = all.seurat[2:length(all.seurat)],
                       add.cell.ids = sample.names, project = 'Human_Innate')

merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT\\.")

## 3. ADD METADATA ####
sample_metadata <- read.csv(file = "/Volumes/Samsung_T5/Human_MAIT_NKT/Data/Metadata.csv") %>% 
  dplyr::select(-sample.id, -X, -file)

rownames(sample_metadata) <- sample_metadata$orig.ident
cell_metadata <- sample_metadata[merged_seurat$orig.ident,]
for(meta in names(cell_metadata)){
  merged_seurat[[meta]] <- cell_metadata[[meta]]
}

merged_seurat@meta.data

## 2. REMOVE CELLS WITH LOW MITOCHONDRIAL CONTENT ####

# Look at qc measures per batch & donor
VlnPlot(merged_seurat, features="percent.mt", group.by = "Batch", split.by = "Donor") + labs(x = "Batch", fill = "Donor")

# Split seurat object to find percent.mt outliers per batch (& subset seurat object)
seur.list <- SplitObject(merged_seurat, split.by = "Batch")
seur.list <- lapply(X = seur.list, FUN = function(x) {
  threshold <- attr(scuttle::isOutlier(x@meta.data$percent.mt, type="higher"), "thresholds")["higher"]
  print(threshold)
  x <- subset(x, subset= percent.mt<threshold)
  # x@meta.data$highmt <- ifelse(x@meta.data$percent.mt > threshold, TRUE, FALSE)
  # return(x)
})
# Combine seurat object back together
filtered_seurat <- purrr::reduce(seur.list, function(x, y){
  merge(x = x, y = y, add.cell.ids = NULL, merge.data = T)
})

# Sanity check
VlnPlot(filtered_seurat, features="percent.mt", 
        group.by="Batch", split.by="Donor") + 
  labs(x = "Batch", fill = "Donor", title="post-filtering")

VlnPlot(filtered_seurat, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        pt.size = 0)

# Filter out cells with less than 500 or more than 3000 expressed genes
filtered_seurat <- subset(filtered_seurat, 
                          subset = nFeature_RNA > 500 & nFeature_RNA < 3000)

# Create violin plots to visualize the distribution of features
# in the filtered dataset
VlnPlot(filtered_seurat, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        pt.size = 0)

# Get the expression counts for each gene and filter out genes
# that are not expressed in at least 20 cells, as well as the
# mitochondrial genes
counts <- GetAssayData(object = filtered_seurat, 
                       slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 20
genes.to.remove <- rownames(filtered_seurat)[grep(rownames(filtered_seurat), pattern = "^MT\\.")]
keep_genes[which(names(keep_genes) %in% genes.to.remove)] = FALSE
filtered_counts <- counts[keep_genes, ]

# Create a new Seurat object with the filtered counts and metadata
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Normalize the data using the log normalization method
filtered_seurat <- NormalizeData(object = filtered_seurat,
                                 normalization.method = "LogNormalize", 
                                 assay = "RNA")

# Identify highly variable genes
filtered_seurat <- FindVariableFeatures(filtered_seurat, 
                                        selection.method = "vst", 
                                        nfeatures= 2000, 
                                        assay = "RNA")

# Scale the data so that it has mean 0 and variance 1
filtered_seurat <- ScaleData(filtered_seurat, 
                             assay = "RNA")

# Perform principal component analysis (PCA) on the scaled data
filtered_seurat <- RunPCA(filtered_seurat, 
                          assay = "RNA",
                          seed.use = 42, 
                          npcs = 50, 
                          weight.by.var = TRUE)

# Create an elbow plot to help choose the number of PCs to use
ElbowPlot(filtered_seurat, ndims = 50, reduction = "pca")

# Use uniform manifold approximation and projection (UMAP) to 
# visualize the data in two dimensions
pcs <- 20
filtered_seurat <- RunUMAP(filtered_seurat, 
                           dims = 1:pcs, 
                           reduction = "pca", 
                           assay = "RNA", 
                           n.neighbors = 30, 
                           seed.use = 29, 
                           reduction.name =  "initial_umap")

# Identify cell clusters using the Louvain algorithm
filtered_seurat <- FindNeighbors(filtered_seurat, 
                                 dims = 1:pcs, 
                                 assay = "RNA")
filtered_seurat <- FindClusters(filtered_seurat, 
                                resolution = 0.7, 
                                algorithm = 2, 
                                random.seed = 29)

# Identify some cell populations to be used to calculate cLISI
# DN_thymocytes
DN_thymocytes <- WhichCells(object = filtered_seurat, expression = PTCRA > 1, slot = 'data')
DN_thymocytes_df <- data.frame(DN_thymocytes, "DN")
rownames(DN_thymocytes_df) <- DN_thymocytes
DN_thymocytes_df$DN_thymocytes <- NULL
colnames(DN_thymocytes_df)[1] <- "cell_type"

# B_cells
B_cells <- WhichCells(object = filtered_seurat, expression = CD19 > 1 & IGKC > 1, slot = 'data')
B_cells_df <- data.frame(B_cells, "B")
rownames(B_cells_df) <- B_cells
B_cells_df$B_cells <- NULL
colnames(B_cells_df)[1] <- "cell_type"

# Tregs
Tregs <- WhichCells(object = filtered_seurat, expression = FOXP3 > 1, slot = 'data')
Tregs_df <- data.frame(Tregs, "Treg")
rownames(Tregs_df) <- Tregs
Tregs_df$Tregs <- NULL
colnames(Tregs_df)[1] <- "cell_type"

# MAIT
MAIT <- WhichCells(object = filtered_seurat, expression = SLC4A10 > 1 & FOXP3 <= 1, slot = 'data')
MAIT_df <- data.frame(MAIT, "MAIT")
rownames(MAIT_df) <- MAIT
MAIT_df$MAIT <- NULL
colnames(MAIT_df)[1] <- "cell_type"

# CD4
CD4 <- WhichCells(object = filtered_seurat, expression = CD4 > 1 & CD8A <= 1 & SLC4A10 <= 1 & FOXP3 <= 1 & CCR7 >1, slot = 'data')
CD4_df <- data.frame(CD4, "CD4")
rownames(CD4_df) <- CD4
CD4_df$CD4 <- NULL
colnames(CD4_df)[1] <- "cell_type"

# DP_thymocytes
DP_thymocytes <- WhichCells(object = filtered_seurat, expression = RAG1 > 1 & CD1C > 1, slot = 'data')
DP_thymocytes_df <- data.frame(DP_thymocytes, "DP")
rownames(DP_thymocytes_df) <- DP_thymocytes
DP_thymocytes_df$DP_thymocytes <- NULL
colnames(DP_thymocytes_df)[1] <- "cell_type"

# CD8AA_thymocytes
CD8AA_thymocytes <- WhichCells(object = filtered_seurat, expression = CD8A > 1 & GNG4 > 1, slot = 'data')
CD8AA_thymocytes_df <- data.frame(CD8AA_thymocytes, "CD8AA")
rownames(CD8AA_thymocytes_df) <- CD8AA_thymocytes
CD8AA_thymocytes_df$CD8AA_thymocytes <- NULL
colnames(CD8AA_thymocytes_df)[1] <- "cell_type"

# Combine the data frames for all identified populations
To_pass_to_metadata <- rbind(DP_thymocytes_df, CD8AA_thymocytes_df, DN_thymocytes_df, Tregs_df, MAIT_df, CD4_df, B_cells_df)

# Add the identified cell populations to metadata column "Test_cell
filtered_seurat <- AddMetaData(filtered_seurat, To_pass_to_metadata, col.name = "Test_cell_ident_for_cLISI")
filtered_seurat@meta.data

# Set the Idents to "Tissue"
Idents(filtered_seurat) <- "Tissue"
# Subset the data based on "Thymus"
filtered_seurat_thymus <- subset(filtered_seurat, idents = "Thymus")
# Subset the data based on "PBMC"
filtered_seurat_PBMC <- subset(filtered_seurat, idents = "PBMC")

# Set the Idents to "Test_cell_ident_for_cLISI"
Idents(filtered_seurat) <- "Test_cell_ident_for_cLISI"
# Subset the data based on the identified cell types for cLISI
filtered_seurat_cLISI <- subset(filtered_seurat, idents = c("DN", "DP", "CD8AA", "CD4", "MAIT", "Treg", "B"))

# get umap coordinates
umap_coords_PBMC <- filtered_seurat_PBMC@reductions$initial_umap@cell.embeddings
umap_coords_thymus <- filtered_seurat_thymus@reductions$initial_umap@cell.embeddings
umap_coords_cLISI <- filtered_seurat_cLISI@reductions$initial_umap@cell.embeddings

# Define the list of metadata variables to be used in the cLISI analysis
list_of_metadata <- c("Donor", "Batch", "Method", "Tissue", "cell.ident", "group.ident", "seurat_clusters", "Test_cell_ident_for_cLISI")

# Compute the before cLISI scores for the thymus and PBMC datasets
lisi_before_results_thymus <- compute_lisi(umap_coords_thymus, filtered_seurat_thymus@meta.data, list_of_metadata)
lisi_before_results_PBMC <- compute_lisi(umap_coords_PBMC, filtered_seurat_PBMC@meta.data, list_of_metadata)

# Compute the before cLISI scores for the custom cell types defined in filtered_seurat_cLISI
clisi_before_results <- compute_lisi(umap_coords_cLISI, filtered_seurat_cLISI@meta.data, list_of_metadata)

# Print a summary of the before cLISI scores for the thymus, PBMC, and custom cell types datasets
summary(lisi_before_results_thymus)
summary(lisi_before_results_PBMC)
summary(clisi_before_results)

# Define color palettes for each group
my_cols <- c(met.brewer("Austria", n = 7, type = "discrete"),
             met.brewer("Egypt", n = 4, type = "discrete"),
             met.brewer("Cross", n = 9, type = "discrete"),
             met.brewer("Troy", n = 2, type = "discrete"),
             met.brewer("Derain", n = 7, type = "discrete"),
             met.brewer("Hiroshige", n = 10, type = "discrete")
)

# Add an extra color to the palette
my_cols <- c(my_cols, "red")

# Define a new color palette
VanGogh_colors <- met.brewer("VanGogh2", type = "discrete", n = 3)

# Set directory to save plots
fig_dir <- "/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/"

# Create plot 1: plot by Method group
pdf(paste0(fig_dir, "Method_before.pdf"), width=7, height=7)
p1 <- DimPlot(filtered_seurat, group.by = c("Method"), 
              pt.size = 0.1, 
              ncol = 1, 
              label = FALSE, label.size = 0, 
              cols = alpha(my_cols, 0.7)) + NoLegend()
print(p1)
dev.off()

# Create plot 2: plot by Donor group
pdf(paste0(fig_dir, "Donor_before.pdf"), width=7, height=7)
p2 <- DimPlot(filtered_seurat, group.by = c("Donor"), 
              pt.size = 0.1, 
              ncol = 1, 
              label = FALSE, label.size = 0, 
              cols = alpha(my_cols, 0.7)) + NoLegend()
print(p2)
dev.off()

# Create plot 3: plot by Batch group
pdf(paste0(fig_dir, "Batch_before.pdf"), width=7, height=7)
p3 <- DimPlot(filtered_seurat, group.by = c("Batch"), 
              pt.size = 0.1, 
              ncol = 1, 
              label = FALSE, label.size = 0, 
              cols = alpha(my_cols, 0.7)) + NoLegend()
print(p3)
dev.off()

# Create plot 3: plot by Tissue group
pdf(paste0(fig_dir, "Tissue_before.pdf"), width=7, height=7)
p4 <- DimPlot(filtered_seurat, group.by = c("Tissue"), 
              pt.size = 0.1, 
              ncol = 1, 
              label = FALSE, label.size = 0, 
              cols = alpha(c("#a40000", "#72bcd5"), 0.3)) + NoLegend()
print(p4)
dev.off()

# Create plot 3: plot by Cluster group
pdf(paste0(fig_dir, "Clusters_before.pdf"), width=7, height=7)
p5 <- DimPlot(filtered_seurat, group.by = c("seurat_clusters"), 
              pt.size = 0.1, 
              ncol = 1, 
              label = TRUE, label.size = 6, 
              cols = alpha(my_cols, 0.7)) + NoLegend()
print(p5)
dev.off()

# Run harmony to integrate the datasets, and perform PCA on RNA data 
# while accounting for batch and method effects
filtered_seurat_harmony <- RunHarmony(filtered_seurat, 
                                      group.by.vars = c("Batch", "Method"), 
                                      assay.use = "RNA", reduction = "pca",
                                      max.iter.harmony = 30)

# Reduce the dimensions of the dataset using UMAP
filtered_seurat_harmony <- RunUMAP(filtered_seurat_harmony, 
                                   dims = 1:pcs, reduction = "harmony", 
                                   assay = "RNA",
                                   seed.use = 29,
                                   reduction.name = "UMAP_50")

# Find the nearest neighbors for each cell based on the UMAP dimensions
filtered_seurat_harmony <- FindNeighbors(filtered_seurat_harmony, 
                                         reduction = "harmony", 
                                         dims = 1:pcs)

# Cluster cells based on the nearest neighbors found earlier 
filtered_seurat_harmony <- FindClusters(filtered_seurat_harmony, 
                                        resolution = 1.2, 
                                        algorithm = 2, 
                                        random.seed = 29)

# Visualize the clusters on the reduced UMAP dimensions
DimPlot(filtered_seurat_harmony, group.by = c("RNA_snn_res.1.2"), 
        pt.size = 0.1, 
        ncol = 1, 
        label = TRUE, 
        label.size = 6,
        reduction = "UMAP_50",
        cols = alpha(my_cols, 0.7)) + NoLegend()

#############################################################################
#The clustree() function will plot a dendrogram of hierarchical clustering and 
#an accompanying plot showing the optimal number of clusters. 
#The plot shows the cophenetic correlation coefficient (CCC) as a function of 
#the number of clusters, with a red line indicating the optimal number of clusters. 
#The cophenetic correlation coefficient measures how well the dendrogram represents 
#the pairwise distances between data points.

clustree(filtered_seurat_harmony, prefix = "RNA_snn_res.")
#############################################################################

# Set the sample/tissue identity for the entire dataset
Idents(filtered_seurat_harmony) <- "Tissue"

# Create a subset of the data with cells only from the "Thymus" tissue
filtered_seurat_Harmony_thymus <- subset(filtered_seurat_harmony, idents = "Thymus")

# Create a subset of the data with cells only from the "PBMC" tissue
filtered_seurat_Harmony_PBMC <- subset(filtered_seurat_harmony, idents = "PBMC")

# Set a new sample/cell identity for the entire dataset for a downstream analysis
Idents(filtered_seurat_harmony) <- "Test_cell_ident_for_cLISI"

# Create a subset of the data with cells only from specific cell types
filtered_seurat_Harmony_cLisi <- subset(filtered_seurat_harmony, 
                                        idents = c("DN", "DP", "CD8AA", "CD4", "MAIT", "Treg", "B"))

# Obtain the UMAP coordinates for each subset of the data
umap_coords_thymus <- filtered_seurat_Harmony_thymus@reductions$UMAP_50@cell.embeddings
umap_coords_PBMC <- filtered_seurat_Harmony_PBMC@reductions$UMAP_50@cell.embeddings
umap_coords_cLISI <- filtered_seurat_Harmony_cLisi@reductions$UMAP_50@cell.embeddings

# Create a list of metadata variables that we will use in downstream analyses
list_of_metadata = c("Donor", "Batch", "Method", "Tissue", "group.ident", "seurat_clusters", "Test_cell_ident_for_cLISI")

# Compute the LISI score for each subset of the data using the UMAP coordinates and metadata
lisi_after_results_thymus <- compute_lisi(umap_coords_thymus, filtered_seurat_Harmony_thymus@meta.data, list_of_metadata)
lisi_after_results_PBMC <- compute_lisi(umap_coords_PBMC, filtered_seurat_Harmony_PBMC@meta.data, list_of_metadata)
clisi_after_results <- compute_lisi(umap_coords_cLISI, filtered_seurat_Harmony_cLisi@meta.data, list_of_metadata)

# Print a summary of the LISI scores for each subset of the data
summary(lisi_after_results_thymus)
summary(lisi_after_results_PBMC)
summary(clisi_after_results)

# Count the number of cells in each Method category in the PBMC dataset
table(filtered_seurat_Harmony_PBMC@meta.data$Method)

# Combine the batch information from before and after the Harmony correction for the thymus dataset
batch_lisi_before_results_thymus <- lisi_before_results_thymus %>% 
  select(Batch) %>% 
  mutate(Before = Batch, .keep = "none", cell.id = rownames(lisi_before_results_thymus))

batch_lisi_after_results_thymus <- lisi_after_results_thymus %>% 
  select(Batch) %>% 
  mutate(After = Batch, .keep = "none", cell.id = rownames(lisi_after_results_thymus))

batch_lisi <- full_join(batch_lisi_before_results_thymus, batch_lisi_after_results_thymus, by = "cell.id")
rownames(batch_lisi) <- batch_lisi$cell.id

batch_lisi <- batch_lisi %>% pivot_longer(-cell.id)

batch_lisi$cell.id <- NULL

# Create a plot of the distribution of batch-corrected LISI scores for the thymus dataset
pdf(paste0(fig_dir, "Batch_Lisi_thymus.pdf"), width=4, height=2)  
p1 <- ggplot(batch_lisi, aes(x = value, y= name, fill = name)) + 
  geom_density_ridges(alpha = 0.6, scale = 1, show.legend = FALSE, bandwidth = 0.1, size = 0.2) +
  labs(x = "iLISI", y = "Density", title = "Batch mixing") + 
  scale_y_discrete(expand = c(-1, 0)) +    
  coord_cartesian(clip = "off") + 
  xlim(0, 7) +
  theme_classic() + NoLegend()
print(p1)
dev.off()

donor_lisi_before_results_thymus <- lisi_before_results_thymus %>% select(Donor) %>% 
  mutate(Before = Donor, .keep = "none", cell.id = rownames(lisi_before_results_thymus))
donor_lisi_after_results_thymus <- lisi_after_results_thymus %>% select(Donor) %>% 
  mutate(After = Donor, .keep = "none", cell.id = rownames(lisi_after_results_thymus))
donor_lisi <- full_join(donor_lisi_before_results_thymus, donor_lisi_after_results_thymus, by = "cell.id")
rownames(donor_lisi) <- donor_lisi$cell.id

donor_lisi <- donor_lisi %>% pivot_longer(-cell.id)

donor_lisi$cell.id <- NULL

pdf(paste0(fig_dir, "Donor_Lisi_thymus.pdf"), width=4, height=2)  
p1 <- ggplot(donor_lisi, aes(x = value, y= name, fill = name)) + 
  geom_density_ridges(alpha = 0.6, scale = 1, show.legend = FALSE, bandwidth = 0.1, size = 0.2) +
  labs(x = "iLISI", y = "Density", title = "Donor mixing") + 
  scale_y_discrete(expand = c(-1, 0)) + 
  coord_cartesian(clip = "off") + 
  xlim(0, 7) +
  theme_classic() + NoLegend()
print(p1)
dev.off()

method_lisi_before_results_thymus <- lisi_before_results_thymus %>% select(Method) %>% 
  mutate(Before = Method, .keep = "none", cell.id = rownames(lisi_before_results_thymus))
method_lisi_after_results_thymus <- lisi_after_results_thymus %>% select(Method) %>% 
  mutate(After = Method, .keep = "none", cell.id = rownames(lisi_after_results_thymus))
method_lisi <- full_join(method_lisi_before_results_thymus, method_lisi_after_results_thymus, by = "cell.id")
rownames(method_lisi) <- method_lisi$cell.id
method_lisi <- method_lisi %>% pivot_longer(-cell.id)
method_lisi$cell.id <- NULL

pdf(paste0(fig_dir, "Method_Lisi_thymus.pdf"), width=4, height=2)  
p1 <- ggplot(method_lisi, aes(x = value, y= name, fill = name)) + 
  geom_density_ridges(alpha = 0.6, scale = 1, show.legend = FALSE, bandwidth = 0.1, size = 0.2) +
  labs(x = "iLISI", y = "Density", title = "Method mixing") + 
  scale_y_discrete(expand = c(-1, 0)) + 
  coord_cartesian(clip = "off") +
  xlim(0, 4) +
  theme_classic() + NoLegend()
print(p1)
dev.off()

batch_lisi_before_results_PBMC <- lisi_before_results_PBMC %>% select(Batch) %>% 
  mutate(Before = Batch, .keep = "none", cell.id = rownames(lisi_before_results_PBMC))
batch_lisi_after_results_PBMC <- lisi_after_results_PBMC %>% select(Batch) %>% 
  mutate(After = Batch, .keep = "none", cell.id = rownames(lisi_after_results_PBMC))
batch_lisi <- full_join(batch_lisi_before_results_PBMC, batch_lisi_after_results_PBMC, by = "cell.id")
rownames(batch_lisi) <- batch_lisi$cell.id

batch_lisi <- batch_lisi %>% pivot_longer(-cell.id)

batch_lisi$cell.id <- NULL

pdf(paste0(fig_dir, "Batch_Lisi_PBMC.pdf"), width=4, height=2)  
p1 <- ggplot(batch_lisi, aes(x = value, y= name, fill = name)) + 
  geom_density_ridges(alpha = 0.6, scale = 1, show.legend = FALSE, bandwidth = 0.1, size = 0.2) +
  labs(x = "iLISI", y = "Density", title = "Batch mixing") + 
  scale_y_discrete(expand = c(-1, 0)) +     
  coord_cartesian(clip = "off") + 
  xlim(0, 7) +
  theme_classic() + NoLegend()
print(p1)
dev.off()

donor_lisi_before_results_PBMC <- lisi_before_results_PBMC %>% select(Donor) %>% 
  mutate(Before = Donor, .keep = "none", cell.id = rownames(lisi_before_results_PBMC))
donor_lisi_after_results_PBMC <- lisi_after_results_PBMC %>% select(Donor) %>% 
  mutate(After = Donor, .keep = "none", cell.id = rownames(lisi_after_results_PBMC))
batch_lisi <- full_join(donor_lisi_before_results_PBMC, donor_lisi_after_results_PBMC, by = "cell.id")
rownames(batch_lisi) <- batch_lisi$cell.id

batch_lisi <- batch_lisi %>% pivot_longer(-cell.id)

batch_lisi$cell.id <- NULL

pdf(paste0(fig_dir, "Donor_Lisi_PBMC.pdf"), width=4, height=2)  
p1 <- ggplot(batch_lisi, aes(x = value, y= name, fill = name)) + 
  geom_density_ridges(alpha = 0.6, scale = 1, show.legend = FALSE, bandwidth = 0.1, size = 0.2) +
  labs(x = "iLISI", y = "Density", title = "Donor mixing") + 
  scale_y_discrete(expand = c(-1, 0)) +    
  coord_cartesian(clip = "off") +
  xlim(0, 7) +
  theme_classic() + NoLegend()
print(p1)
dev.off()

method_lisi_before_results_PBMC <- lisi_before_results_PBMC %>% select(Method) %>% 
  mutate(Before = Method, .keep = "none", cell.id = rownames(lisi_before_results_PBMC))
method_lisi_after_results_PBMC <- lisi_after_results_PBMC %>% select(Method) %>% 
  mutate(After = Method, .keep = "none", cell.id = rownames(lisi_after_results_PBMC))
batch_lisi <- full_join(method_lisi_before_results_PBMC, method_lisi_after_results_PBMC, by = "cell.id")
rownames(batch_lisi) <- batch_lisi$cell.id
batch_lisi <- batch_lisi %>% pivot_longer(-cell.id)
batch_lisi$cell.id <- NULL

pdf(paste0(fig_dir, "Method_Lisi_PBMC.pdf"), width=4, height=2)  
p1 <- ggplot(batch_lisi, aes(x = value, y= name, fill = name)) + 
  geom_density_ridges(alpha = 0.6, scale = 1, show.legend = FALSE, bandwidth = 0.1, size = 0.2) +
  labs(x = "iLISI", y = "Density", title = "Method mixing") + 
  scale_y_discrete(expand = c(-1, 0)) +
  coord_cartesian(clip = "off") +
  xlim(0, 4) +
  theme_classic() + NoLegend()
print(p1)
dev.off()

clisi_before <- clisi_before_results %>% select(Test_cell_ident_for_cLISI) %>% 
  mutate(Before = Test_cell_ident_for_cLISI, .keep = "none", cell.id = rownames(clisi_before_results))
clisi_after <- clisi_after_results %>% select(Test_cell_ident_for_cLISI) %>% 
  mutate(After = Test_cell_ident_for_cLISI, .keep = "none", cell.id = rownames(clisi_after_results))
clisi_lisi <- full_join(clisi_before, clisi_after, by = "cell.id")
rownames(clisi_lisi) <- clisi_lisi$cell.id
clisi_lisi <- clisi_lisi %>% pivot_longer(-cell.id)
clisi_lisi$cell.id <- NULL

pdf(paste0(fig_dir, "cLisi.pdf"), width=4, height=2)  
p1 <- ggplot(clisi_lisi, aes(x = value, y= name, fill = name)) + 
  geom_density_ridges(alpha = 0.6, scale = 1, show.legend = FALSE, bandwidth = 0.1, size = 0.2) +
  labs(x = "cLISI", y = "Density", title = "cLisi") + 
  scale_y_discrete(expand = c(-1, 1)) + 
  coord_cartesian(clip = "off", xlim=c(0, 4)) + 
  xlim(0, 4) +
  scale_x_continuous(breaks = seq(0, 4, 1)) +
  theme_classic() + NoLegend()
print(p1)
dev.off()

# plot 1: Method
pdf(paste0(fig_dir, "Method_after.pdf"), width=7, height=7)
p1 <- DimPlot(filtered_seurat_harmony, group.by = "Method", pt.size = 0.1, ncol = 1, label = FALSE, label.size = 0, reduction = "UMAP_50", cols = alpha(my_cols, 0.7)) + NoLegend()
print(p1)
dev.off()

# plot 2: Donor
pdf(paste0(fig_dir, "Donor_after.pdf"), width=7, height=7)
p2 <- DimPlot(filtered_seurat_harmony, group.by = "Donor", pt.size = 0.1, ncol = 1, label = FALSE, label.size = 0, reduction = "UMAP_50", cols = alpha(my_cols, 0.7)) + NoLegend()
print(p2)
dev.off()

# plot 3: Batch
pdf(paste0(fig_dir, "Batch_after.pdf"), width=7, height=7)
p3 <- DimPlot(filtered_seurat_harmony, group.by = "Batch", pt.size = 0.1, ncol = 1, label = FALSE, label.size = 0, reduction = "UMAP_50", cols = alpha(my_cols, 0.7)) + NoLegend()
print(p3)
dev.off()

# plot 4: Tissue
pdf(paste0(fig_dir, "Tissue_after.pdf"), width=7, height=7)
p4 <- DimPlot(filtered_seurat_harmony, group.by = "Tissue", pt.size = 0.1, ncol = 1, label = FALSE, label.size = 0, reduction = "UMAP_50", cols = alpha(c("#a40000", "#72bcd5"), 0.3)) + NoLegend()
print(p4)
dev.off()

# plot 5: Clusters
pdf(paste0(fig_dir, "Clusters_after.pdf"), width=7, height=7)
p5 <- DimPlot(filtered_seurat_harmony, group.by = "seurat_clusters", 
              pt.size = 0.1, ncol = 1, 
              label = TRUE, label.size = 6, 
              reduction = "UMAP_50", cols = alpha(my_cols, 0.7)) + NoLegend()
print(p5)
dev.off()

# Count number of cells in filtered_seurat_harmony object
sum(table(filtered_seurat_harmony@meta.data$orig.ident))

# Create new_clusters column in filtered_seurat_harmony object
filtered_seurat_harmony$new_clusters <- case_when(
  filtered_seurat_harmony$RNA_snn_res.1.2 == '13' ~ '0',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '16' ~ '1',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '19' ~ '1',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '17' ~ '2',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '3' ~ '3',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '14' ~ '4',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '12' ~ '5',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '9' ~ '6',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '8' ~ '7',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '11' ~ '8',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '5' ~ '9',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '4' ~ '10',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '6' ~ '11',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '20' ~ '12',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '7' ~ '13',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '2' ~ '14',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '0' ~ '15',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '1' ~ '16',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '15' ~ '17',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '10' ~ '18',
  filtered_seurat_harmony$RNA_snn_res.1.2 == '18' ~ '19'
)

# Convert new_clusters to a factor with levels 0 to 19
filtered_seurat_harmony$new_clusters <- as.factor(filtered_seurat_harmony$new_clusters)
filtered_seurat_harmony$new_clusters <- factor(filtered_seurat_harmony$new_clusters, levels = c(
  "0", "1", "2", "3", "4", "5", "6", "7", "8",
  "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"))

colors_clusters <- c("0" = "#f4c40f", "1" = "#b75347", "2" = "#d8443c", "3" = "#e09351", "4" = "#2b9b81", 
                     "5" = "#421401", "6" = "#92c051", "7" = "#9f5691", "8" = "#17154f", "9" = "#74c8c3", 
                     "10" = "#5a97c1", "11" = "gold", "12" = "#a40000", "13" = "#72bcd5", "14" = "grey50",
                     "15" = "orange", "16" = "blueviolet", "17" = "#0a2e57", "18" = "red", "19" = "#bdbdbd")

# Set the device to produce a PDF file with a given name and size
pdf(paste0(fig_dir, "Renamed_clusters.pdf"), width=7, height=7)

# Create a dim plot of the filtered Seurat object using the "new_clusters" grouping 
# variable, with specified point size, number of columns, label display, label size, 
# point color and repulsion options, and no legend
p <- DimPlot(filtered_seurat_harmony, group.by = c("new_clusters"), 
             pt.size = 0.1, ncol = 1, 
             label = TRUE, label.size = 6, 
             reduction = "UMAP_50", cols = alpha(colors_clusters, 0.7)) + NoLegend()

# Print the dim plot
print(p)

# Close the PDF device
dev.off()

# Set the Idents of the Seurat object to "new_clusters"
Idents(filtered_seurat_harmony) <- "new_clusters"

# Subset the Seurat object to remove cells with "19" as their identity
filtered_seurat_harmony <- subset(filtered_seurat_harmony, idents = "19", invert = TRUE)

# Convert the "new_clusters" column to a factor with specific levels
filtered_seurat_harmony$new_clusters <- factor(filtered_seurat_harmony$new_clusters, levels = c(
  "0", "1", "2", "3", "4", "5", "6", "7", "8",
  "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))

# Assign the "new_clusters" grouping variable to the Seurat object Idents
Idents(filtered_seurat_harmony) <- "new_clusters"

# Set the device to produce a PDF file with a given name and size
pdf(paste0(fig_dir, "Splitted_by_group.ident.pdf"), width=70, height=7)

# Create a dim plot of the filtered Seurat object using the "group.ident" splitting 
# variable, with specified point size, number of columns, label display, label size, 
# point color and repulsion options, and no legend
p <- DimPlot(filtered_seurat_harmony, split.by = "group.ident",
             pt.size = 0.1, 
             ncol = 12, 
             label = FALSE, repel = FALSE, label.size = 6, na.value = "grey95", reduction = "UMAP_50",
             cols = alpha(colors_clusters, 0.7)) + NoLegend()

# Print the dim plot
print(p)

# Close the PDF device
dev.off()

# Assign cell.ident column to "Tissue" and create subsets for PBMC and Thymus cells
Idents(filtered_seurat_harmony) <- "Tissue"
PBMC_cells <- subset(filtered_seurat_harmony, idents = "PBMC")
Thymus_cells <- subset(filtered_seurat_harmony, idents = "Thymus")

# Create dataframes of the number of cells in each cluster for PBMC and Thymus cells
A <- as.data.frame(table(PBMC_cells@meta.data$new_clusters))
B <- as.data.frame(table(Thymus_cells@meta.data$new_clusters))

# Rename columns and join the dataframes
A <- A %>% dplyr::rename(clusters = Var1, number_pbmc = Freq)
B <- B %>% dplyr::rename(clusters = Var1, number_thymus = Freq)
tissue_per_cluster <- inner_join(A, B, by = "clusters") %>% pivot_longer(-clusters)

# Define breaks for y-axis and create theme
brks <- c(0, 0.25, 0.5, 0.75, 1)
My_Theme = theme(
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  axis.title.y = element_text(size = 18))

# Create barplot of the frequency of each cluster in PBMC and Thymus cells
pdf(paste0(fig_dir, "Cells_per_tisue.pdf"), width=12, height=7)
p6 <- ggplot(tissue_per_cluster, aes(x = clusters, y = value, fill=name, na.rm = TRUE)) +
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous(breaks = brks, labels = scales::percent(brks), expand = c(0,0)) +
  labs(x = "Clusters", y = "Frequency (%)", fill = "Tissue") +
  scale_fill_discrete(name = "Tissue", labels = c("PBMC","Thymus")) +
  scale_fill_manual(values=alpha(c("#a40000", "#72bcd5"), 0.9)) +
  theme_classic(base_size = 14) + My_Theme
print(p6)
dev.off()

##################################################################################################################
# Define a function "bootstrap_myclusters" that takes in a dataset "x", a function "FUN", 
# a vector of cluster IDs "clusters", a boolean "transposed", number of cells "n.cells", number of iterations "iterations", 
# and additional parameters "..." as input.

bootstrap_myclusters <- function(x, FUN, clusters=NULL, transposed=FALSE, n.cells=5000, iterations=50, ...) {
  # If no cluster IDs are provided, generate cluster IDs using the function "FUN" and assign them to "clusters".
  if (is.null(clusters)) {
    clusters <- FUN(x, ...)
  }
  # Create a vector of unique cluster IDs and assign them to "cluster.ids".
  cluster.ids <- as.character(sort(unique(clusters)))
  # Create a matrix of zeros with dimensions based on the length of "cluster.ids" and assign it to "output".
  output <- matrix(0, length(cluster.ids), length(cluster.ids))
  # Set the lower triangular portion of "output" to missing values.
  output[lower.tri(output)] <- NA_real_
  # Set the row and column names of "output" to "cluster.ids".
  dimnames(output) <- list(cluster.ids, cluster.ids)
  
  # Repeat the following "iterations" number of times.
  for (i in seq_len(iterations)) {
    # If "transposed" is false, randomly select "n.cells" number of columns from "x" and assign them to "chosen".
    # If "transposed" is true, randomly select "n.cells" number of rows from "x" and assign them to "chosen".
    if (transposed) {
      chosen <- sample(nrow(x), n.cells)
      resampled <- x[chosen,]
    } else {
      chosen <- sample(ncol(x), n.cells)
      resampled <- x[,chosen]
    }
    # Generate cluster IDs for the resampled data using "FUN" and assign them to "reclusters".
    reclusters <- FUN(resampled, ...)
    # Create a contingency table "tab" with "clusters[chosen]" as rows and "reclusters" as columns.
    tab <- table(clusters[chosen], reclusters)
    # For each pair of cluster IDs in "cluster.ids", calculate the spread between their respective rows in "tab" and add the result to "output".
    for (j1 in seq_along(cluster.ids)) {
      spread1 <- tab[cluster.ids[j1],]
      spread1 <- spread1/sum(spread1)
      for (j2 in seq_len(j1)) {
        spread2 <- tab[cluster.ids[j2],]
        spread2 <- spread2/sum(spread2)
        output[j2,j1] <- output[j2,j1] + sum(spread1 * spread2)/iterations
      }
    }
  }
  # Return the "output" matrix.
  output
}

# Define a function "myknn_FUN" that takes in a dataset "x" and returns cluster IDs based on dimension reduction and clustering methods.

myknn_FUN <- function(x) {
  # Perform dimension reduction on the dataset using the "FindNeighbors" function with parameters "verbose = F", "reduction = "harmony"", and "dims=1:15".
  g <- FindNeighbors(x, verbose = F, reduction = "harmony", dims=1:15) 
  # Cluster the reduced dataset using the "FindClusters" function with parameters "verbose = F" and "resolution = 1".
  g <- FindClusters(g, verbose = F, resolution = 1)
  # Convert the cluster IDs to numeric values and return them
  as.numeric(g$new_clusters)}

#This is the cluster or CellType information, you already have stored in Seurat object
originals<- filtered_seurat_harmony$new_clusters
#You can choose iterations of your choice but its memory intensive
coassign <- bootstrap_myclusters(filtered_seurat_harmony, clusters = originals, FUN = myknn_FUN, 
                                 n.cells = ncol(filtered_seurat_harmony)-1, iterations = 50) 

#Plot heatmap of co-assignment probabilities
pheatmap(coassign, cluster_row=F, cluster_col=F, main= "Coassignment probabilities", angle_col = 45,
         color=rev(viridis::magma(100)))

#############################################################################################################

# Define a list of genes representing the Egress signature
Egress <- list(c("KLF2", "CORO1A", "CCR7", "CXCR4", "CXCR6", "FOXO1", "CXCR3", "S1PR1", "S1PR4",
                 "S100A4", "S100A6", "EMP3"))

# Compute a module score for each cell in the Seurat object based on the expression of the Egress signature genes
filtered_seurat_harmony <- AddModuleScore(object = filtered_seurat_harmony, features = Egress,
                                          assay = "RNA", name = 'Egress_score')

# Generate a feature plot that shows the distribution of the Egress score across cells, split by tissue type
pdf(paste0(fig_dir, "Egress_signature.pdf"), width=13, height=8)
p <- SCpubr::do_FeaturePlot(sample = filtered_seurat_harmony, 
                            features = "Egress_score1", split.by = "Tissue",
                            plot.title = "",
                            reduction = "UMAP_50",
                            viridis_color_map = "inferno")

# Print the plot to the PDF device
print(p)

# Close the PDF device
dev.off()

##########
# Define a vector of genes
genes_teichmann <- c("PTCRA", "RAG1", "CD1C", "AQP3", "CD8A", "PDCD1", "GNG4", "TRDC", "TRGC2", 
                     "EGR1", "EGR3", "NR4A1", "IKZF4", "FOXP3", "CTLA4", "STAT1", "IFI6", "CD4", "CD40LG",
                     "CD8B", "CCR9", "SATB1", "CCR7", "SELL", "FOS", "JUN", "JUNB",
                     "KLRB1", "ZBTB16", "GZMK", "NKG7", "EOMES", "IFNG", "TBX21", "GNLY", "GZMB", "RORC", "CCR6")

# Assign a new cluster identity to the Seurat object
Idents(filtered_seurat_harmony) <- "new_clusters"

# Scale the expression data of the Seurat object using the defined vector of genes
filtered_seurat_Harmony <- ScaleData(filtered_seurat_harmony, features = genes_teichmann)

# Create a dot plot for each cluster using the defined vector of genes
# The dot plot displays the gene expression levels for each gene in each cluster
pdf(paste0(fig_dir, "Teichmann_features_per_clusters.pdf"), width=7, height=7)
DotPlot_colors <- c("#FFDAB9", "#a40000")
p <- DotPlot(filtered_seurat_harmony, features = genes_teichmann, dot.scale = 8, cols = DotPlot_colors,
             col.min = -1, col.max = 2, dot.min = 0) + coord_flip()
print(p)
dev.off()

# Define the lists of genes for each module
effectorness_genes <- c("HOPX", "GZMB", "GZMK", "ZEB2", "NKG7", "GNLY", "TBX21", "EOMES", "TYROBP", "PRF1", 
                        "CCL4", "CCL5", "KLRB1", "GZMH", "GZMA", "KLRD1", "CST7", "KLF6", "CXCR4")

naive_genes <- c("SATB1", "TCF7", "LEF1", "CCR7", "SELL", "MYC", "EIF3E", "SOX4", "ID3", "BACH2")

type_3_genes <- c("ZBTB16", "MAF", "BLK", "RORA", "RORC", "AHR", "CCR6", "IL7R", "IL18R1", "IL23R", "IL1R1", 
                  "SCART1", "S100A4", "S100A6")

# Calculate the module scores and add them to the Seurat object
filtered_seurat_Harmony <- AddModuleScore(object = filtered_seurat_Harmony, features = effectorness_genes, 
                                          assay = "RNA", name = "Effectorness")
filtered_seurat_Harmony <- AddModuleScore(object = filtered_seurat_Harmony, features = naive_genes, 
                                          assay = "RNA", name = "Naiveness")
filtered_seurat_Harmony <- AddModuleScore(object = filtered_seurat_Harmony, features = type_3_genes, 
                                          assay = "RNA", name = "Type_3_score")

# Create a PDF plot of the module scores
pdf(file = "Signatures.pdf", width = 25, height = 7)
p2 <- FeaturePlot(object = filtered_seurat_Harmony, reduction = "UMAP_50", 
                  features = c("Effectorness", "Naiveness", "Type_3_score"), 
                  cols = alpha(c("#1E466E", "gold", "#a40000"), 1), pt.size = 0.1, order = TRUE, ncol = 3)
print(p2)
dev.off()

###
# Assign a cluster identity to each cell based on the new_clusters assay
Idents(filtered_seurat_harmony) <- "new_clusters"

# Find marker genes for each cluster using the 'wilcox' test and filter for positive markers
clusters.markers <- FindAllMarkers(filtered_seurat_harmony, test.use = 'wilcox', 
                                   logfc.threshold = 0.4, min.pct = 0.3, only.pos = TRUE)

# Select the top 5 marker genes with the highest average log fold change for each cluster
top5 <- clusters.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Keep only the unique gene names from the top 5 genes
top5_distinct <- base::unique(top5$gene)

# Scale the data for the top 5 marker genes and create a DotPlot for visualization
filtered_seurat_harmony <- ScaleData(filtered_seurat_harmony, features = top5_distinct)
DotPlot_colors <- c("lightsteelblue1", "#1E466E")

pdf(file = "Cluster_signatures.pdf", width = 20, height = 7)
p <- DotPlot(filtered_seurat_harmony, features = top5_distinct, dot.scale = 8, cols = DotPlot_colors,
        col.min = -1, col.max = 2, dot.min = 0) + RotatedAxis()
print(p)
dev.off()

write.csv(clusters.markers, "/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/clusters.LogNorm.markers.csv")

FeaturePlot(filtered_seurat_harmony, features = c("ZBTB16"), split.by = "group.ident", 
            pt.size = 0.1, 
            ncol = 4, 
            label = FALSE, repel = FALSE, label.size = 6,
            order = T,
            cols = alpha(c("#1E466E", "gold", "#a40000"), 1)) + NoLegend()

### Add TCR metadata
TCRS_All_cells <- read.csv("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/TCRS_All_cells.csv", row.names = 1)

filtered_seurat_harmony <- AddMetaData(filtered_seurat_harmony, TCRS_All_cells, col.name = NULL)
filtered_seurat_harmony@meta.data

###################### SAVE INTEGRATED SEURAT OBJECT #######################
Idents(filtered_seurat_harmony) <- "orig.ident"
remove <- c("CD1a_1_Thymus", "CD1c_1_Thymus")
filtered_seurat_harmony <- subset(filtered_seurat_harmony, idents=remove, invert = TRUE)

saveRDS(filtered_seurat_harmony, "/Volumes/Samsung_T5/Human_MAIT_NKT/Data/seurat_filtered_harmony_02_15_23.RDS")
###########################################################################





