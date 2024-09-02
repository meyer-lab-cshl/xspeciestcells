# Purpose:
# Author: Salomé Carcy
# Date:




# **************
# 1. IMPORT ####
# **************

# Import librairies
library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)
library(ShinyCell)

# Import data
seur.human <- readRDS("./data_geo/seurat_objects/seurat_human_integrated_object_23_12_01.rds")
seur.ms <- readRDS("./data_geo/seurat_objects/seurat_mouse_integrated_object_23_12_02.rds")



# *****************
# 2. SHINY APP ####
# *****************

#___________________________
## 2.1. Human page ####

# Keep only the integrated umap embedding for shiny app
seur.human@reductions$pca <- NULL
seur.human@reductions$initial_umap <- NULL
seur.human@reductions$harmony <- NULL
seur.human[["umap"]] <- CreateDimReducObject(embeddings = seur.human@reductions$umap_integrated@cell.embeddings,
                                             loadings   = seur.human@reductions$umap_integrated@feature.loadings,
                                             projected  = seur.human@reductions$umap_integrated@feature.loadings.projected,
                                             key="UMAP_",
                                             assay="RNA")
seur.human@reductions$umap_integrated <- NULL


# Make donor a factor
seur.human@meta.data$donor_id <- factor(seur.human@meta.data$donor_id, levels=1:13)
seur.human@meta.data$GEP_with_max_usage <- factor(seur.human@meta.data$GEP_with_max_usage, levels=paste0("GEP", 1:11))

# Create shiny config human
scConf_hu = createConfig(seur.human, maxLevels = 60)
scConf_hu = modColours(scConf_hu, meta.to.mod = "tcell_lineage", new.colours = c("#74c476", "#df65b0", "#08519c", "#9e9ac8", "#9ecae1"))
scConf_hu = modColours(scConf_hu, meta.to.mod = "tcell_lineage_tissue",
                       new.colours = c("#74c476", "#c7e9c0",
                                       "#df65b0", "#d4b9da",
                                       "#08519c", "#4292c6",
                                       "#9e9ac8", "#dadaeb",
                                       "#9ecae1", "#deebf7")
)
scConf_hu = modColours(scConf_hu, meta.to.mod = "tissue", new.colours = c("#b2182b", "#9ecae1"))
scConf_hu = modColours(scConf_hu, meta.to.mod = "clusters_integrated_data", new.colours = as.vector(cols_integrated))
scConf_hu = modColours(scConf_hu, meta.to.mod = "clusters_per_lineage", new.colours = as.vector(c(cols_pbmc_cd4, #0 to 5
                                                                                                  cols_thym_cd4, #0 to 6
                                                                                                  cols_pbmc_cd8, #0 to 4
                                                                                                  cols_thym_cd8, #0 to 5
                                                                                                  cols_pbmc_gdt, #0 to 4
                                                                                                  cols_thym_gdt, #0 to 7
                                                                                                  cols_pbmc_nkt, #0 to 3
                                                                                                  cols_thym_nkt, #0 to 6
                                                                                                  cols_pbmc_mait, #0 to 3
                                                                                                  cols_thym_mait #0 to 6
)))
scConf_hu = modColours(scConf_hu, meta.to.mod = "GEP_with_max_usage", new.colours = as.vector(cols_GEPs))

# set default metadata to display
scConf_hu = modDefault(scConf_hu, "clusters_integrated_data", "tissue")

# create data files
makeShinyFiles(obj=seur.human,
               scConf=scConf_hu,
               gex.assay = "RNA",
               gex.slot = "data",
               gene.mapping = FALSE,
               shiny.prefix = "sc_hu",
               shiny.dir = "shinyAppMulti/",
               default.gene1 = "CD4",
               default.gene2 = "CD8A",
               default.multigene = c("CCR6","RORC","GZMB","GNLY","TBX21",
                                     "IFNG","EOMES","NKG7","GZMK","ZBTB16",
                                     "KLRB1", "JUNB", "JUN", "FOS", "SELL",
                                     "CCR7", "SATB1", "CCR9", "CD8B", "CD40LG",
                                     "CD4", "IFI6", "STAT1", "CTLA4", "FOXP3", "IKZF4",
                                     "NR4A1", "EGR3", "EGR1", "TRGC2", "TRDC", "GNG4",
                                     "PDCD1", "CD8A", "AQP3", "CD1C", "RAG1", "PTCRA"),
               default.dimred = c("UMAP_1", "UMAP_2"))
## /end ####


#___________________________
## 2.2. Mouse page ####

# Keep only the integrated umap for the shiny app
seur.ms@reductions$mnn <- NULL

# Create shiny config mouse
scConf_ms = createConfig(seur.ms)
scConf_ms = modColours(scConf_ms, meta.to.mod = "tcell_lineage", new.colours = c("#08519c", "#9e9ac8", "#9ecae1"))
scConf_ms = modColours(scConf_ms,
                       meta.to.mod = "clusters_integrated_data",
                       new.colours = c("#f4c40f","#b75347", "#d8443c", "#e09351",
                                       "#2b9b81", "#421401", "#92c051", "#9f5691",
                                       "#17154f", "#74c8c3", "#5a97c1", "gold", "#a40000")
)
scConf_ms = modColours(scConf_ms,
                       meta.to.mod = "clusters_annotation",
                       new.colours = c("#f4c40f","#b75347", "#d8443c", "#e09351", "#421401",
                                       "#92c051", "#17154f", "#9f5691", "#5a97c1", "#a40000")
)
scConf_ms = modColours(scConf_ms,
                       meta.to.mod = "inkt_mait_tcr",
                       new.colours = c("grey80", "#9ecae1", "#9e9ac8")
)

# set default metadata to display
scConf_ms = modDefault(scConf_ms, "clusters_integrated_data", "study")

# create data files
makeShinyFiles(obj=seur.ms,
               scConf=scConf_ms,
               gex.assay = "RNA",
               gex.slot = "data",
               gene.mapping = FALSE,
               shiny.prefix = "sc_ms",
               shiny.dir = "shinyAppMulti/",
               default.gene1 = "Cd4",
               default.gene2 = "Cd8a",
               default.multigene = c("Ccl5", "Fosb", "Rorc", "Gzmb", "Nkg7", "Zbtb16",
                                     "Tesc", "Lef1", "Tox", "Gzma", "S1pr1", "Klf2"),
               default.dimred = c("UMAP_1", "UMAP_2"))
## /end ####


#___________________________
## 2.3. Shiny app ####

# Code for shiny app
citation = list(
  author  = "L.Loh, S.Carcy, et al.",
  title   = "",
  journal = "",
  volume  = "",
  page    = "",
  year    = "", 
  doi     = "",
  link    = "")
makeShinyCodesMulti(
  shiny.title = "ShinyCell Human & Mouse innate T cell development",
  shiny.footnotes = citation,
  shiny.prefix = c("sc_hu", "sc_ms"),
  shiny.headers = c("Human data", "Mouse data"), 
  shiny.dir = "shinyAppMulti/")

