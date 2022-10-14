# ---------------------------------
#         RUN SHINY APP
# ---------------------------------

# Import seurat objects
library(SingleCellExperiment)
library(Seurat)
library(ShinyCell)


# Create shiny configuration
scConf = createConfig(seur.harmony)
makeShinyApp(seur.harmony, scConf, gene.mapping = TRUE,
             shiny.title = "ShinyCell Mouse iNKT+Thymocyte data") 

filtered_seurat_MNN <- FindVariableFeatures(filtered_seurat_MNN) # return 2,000 HVGs
scConf_mnn = createConfig(filtered_seurat_MNN)
makeShinyFiles(filtered_seurat_MNN, scConf_mnn, gex.assay="RNA",
               shiny.prefix="fastMNN", shiny.dir="shinyAppHumaniNKT")
