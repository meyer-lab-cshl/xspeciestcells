# ---------------------------------
#         RUN SHINY APP
# ---------------------------------

# Import seurat objects
library(SingleCellExperiment)
library(Seurat)
library(ShinyCell)


# Create shiny configuration
filtered_seurat_MNN <- FindVariableFeatures(filtered_seurat_MNN) # return 2,000 HVGs
scConf_mnn = createConfig(filtered_seurat_MNN)
makeShinyFiles(filtered_seurat_MNN, scConf_mnn, gex.assay="RNA",
               shiny.prefix="fastMNN", shiny.dir="shinyAppHumaniNKT")

scConf_harm = createConfig(filtered_seurat_Harmony)
makeShinyFiles(filtered_seurat_Harmony, scConf_harm, gex.assay="SCT",
               shiny.prefix="Harmony", shiny.dir="shinyAppHumaniNKT")


# Make the shiny app
makeShinyCodesMulti(
  shiny.title = "scRNAseq human iNKT data",
  shiny.footnotes = "",
  shiny.prefix  = c("fastMNN", "Harmony"),
  shiny.headers = c("LogNormalize + fastMNN", "SCTransform + Harmony"), 
  shiny.dir = "shinyAppHumaniNKT/")
