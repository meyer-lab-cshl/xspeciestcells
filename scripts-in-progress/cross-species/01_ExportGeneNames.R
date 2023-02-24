###
# Purpose: Export list of mouse, human, pig genes
# Date: Aug 17th 2022
# Author: Salomé Carcy
###

library(Seurat)
library(Matrix)
library(data.table)


#### MOUSE DATA ####
path.data <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/mouse_data"
mouse1 <- Read10X(file.path(path.data, "B6_1"))
mouse2 <- Read10X(file.path(path.data, "B6_2"))

# Create Seurat Object
seur.ms1 <- CreateSeuratObject(mouse1, project="B6_Thymus_NKT_1")
seur.ms2 <- CreateSeuratObject(mouse2, project="B6_Thymus_NKT_2")
seur.mouse <- merge(seur.ms1, y = seur.ms2, add.cell.ids = c('B61', 'B62'), project = "NKT")

# Export list of genes
genes.ms <- rownames(seur.mouse) # 31,053 genes
fwrite(list(genes.ms), file = "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/01_ListGenes/mouse_genes.txt")



#### HUMAN DATA ####
path.data <- "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/human_data"
human5 <- read.csv(file.path(path.data, "CUThy13_220225_SampleTag05_hs_NKT_RSEC_MolsPerCell.csv"), sep=",", header=T, skip=7) # 404 cells
human8 <- read.csv(file.path(path.data, "CUTHY11BDRscRNA_seq_091621_SampleTag08_hs_NKT_RSEC_MolsPerCell.csv"), sep=",", header=T, skip=8) # 1913 cells
human12 <- read.csv(file.path(path.data, "CUTHY12BDRscRNA_seq_211101_SampleTag12_hs_NKT_RSEC_MolsPerCell.csv"), sep=",", header=T, skip=8) # 344 cells

# Invert the dataframes (cells as columns and genes as rows) and transform to dgCMatrix
convert_to_matrix <- function(df){
  rownames(df) <- df$Cell_Index
  df$Cell_Index <- NULL
  df <- t(df)
  df <- Matrix(df, sparse=T)
  return(df)
}
mat.hu5  <- convert_to_matrix(human5) # 404 cells
mat.hu8  <- convert_to_matrix(human8) # 1913 cells
mat.hu12 <- convert_to_matrix(human12) # 344 cells

# Create Seurat Object
seur.h5  <- CreateSeuratObject(mat.hu5, project="Hu_Thymus_NKT_5")
seur.h8  <- CreateSeuratObject(mat.hu8, project="Hu_Thymus_NKT_8")
seur.h12 <- CreateSeuratObject(mat.hu12, project="Hu_Thymus_NKT_12")
seur.human <- merge(seur.h5, y=c(seur.h8, seur.h12), add.cell.ids=c('Hu5', 'Hu8', 'Hu12'), project="HuNKT")

# Export list of genes
genes.hu <- rownames(seur.human) # 28,479 genes
fwrite(list(genes.hu), file = "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/01_ListGenes/human_genes.txt")



#### PIG DATA ####
pig <- Read10X("~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/raw_data/pig_data/GSE192520_RAW/iNKT")

# Create Seurat Object
seur.pig <- CreateSeuratObject(counts = pig, project = "pig_Thymus_NKT", min.cells = 3, min.features = 200)

# Export list of genes
genes.ss <- rownames(seur.pig) # 13,922 genes
fwrite(list(genes.ss), file = "~/Projects/20220809_Thymic-iNKT-CrossSpecies/data/01_ListGenes/pig_genes.txt")

