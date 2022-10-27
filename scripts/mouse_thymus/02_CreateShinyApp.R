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
