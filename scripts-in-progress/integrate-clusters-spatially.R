## libraries ####
#################
library(Seurat)
#library(SeuratData)
library(ggplot2)
#library(patchwork)
library(dplyr)
library(harmony)
library(Matrix)
library(RColorBrewer)
library(pals)
library(NMF)

#################
## functions ####
#################

analyzeSPOTS <- function(slice, spe, modls=NULL, sce=NULL, cell_annot=NULL,
                         markers=NULL, variable_features=NULL) { 
  tmp <- subset(spe, subset=orig.ident==slice)
  if (is.null(modls)) {
    res.tmp <- SPOTlight(
      x = sce,
      y = tmp,
      groups = cell_annot,
      mgs = markers,
      hvg = variable_features,
      weight_id = "avg_log2FC",
      group_id = "cell_annot",
      gene_id = "gene")
  } else {
    res.tmp <- runDeconvolution(
      x = tmp,
      mod = modls[["mod"]],
      ref = modls[["topic"]])
  }
}

annotateSPOTS <- function(slice, mat, spe, cell_type) { 
  tmp <- subset(spe, subset=orig.ident==slice)
  spot_annotation <- tibble(spotid=names(tmp$Location), location=tmp$Location)
  
  cell_location <- mat %>%
    as.data.frame %>%
    rownames_to_column("spotid") %>%
    left_join(spot_annotation) %>%
    select(spotid, location, everything()) %>%
    pivot_longer(starts_with(cell_type),
                 names_to = "cluster", values_to="proportion") %>%
    #group_by(location, cluster) %>%
    #summarize(total_proportion = sum(proportion)) %>%
    #group_by(location) %>%
    #mutate(total_proportion = total_proportion/sum(total_proportion),
    mutate(celltype=cell_type,
           slice=slice)
}
############
## data ####
############
colvalues <- colorRampPalette(brewer.pal(12,"Paired"))(18)

thymus <- readRDS("data/seurat_filtered_harmony_02_15_23_thymus_subset_annotations.RDS")
markers.subsets <- readr::read_csv("results/top20markers_thymic_populations.csv")
hvgs.subsets <- readr::read_csv("results/hvgs_thymic_populations.csv")

# From: Heimli M, Flåm ST, Hjorthaug HS, Trinh D, Frisk M, Dumont K-A, Ribarska
# T, Tekpli X, Saare M and Lie BA (2023) Multimodal human thymic pro!ling
# reveals trajectories and cellular milieu for T agonist selection.Front.
# Immunol doi: 10.3389/fimmu.2022.1092028

tspe <- readRDS('data/annotated_spatial.rds')

# format spatial object 
# match up rownames: substitute dash with period to reflect naming in thymus
sct_names <- rownames(tspe@assays$SCT@data)
spatial_names <- rownames(tspe@assays$Spatial@data)
integrated_names <- rownames(tspe@assays$integrated@data)

new_sct <- gsub("-", ".", sct_names)
new_spatial <- gsub("-", ".", spatial_names)
new_integrated <- gsub("-", ".", integrated_names)

rownames(tspe@assays$Spatial@data) <- new_spatial
rownames(tspe@assays$Spatial@counts) <- new_spatial
rownames(tspe@assays$SCT@data) <- new_sct
rownames(tspe@assays$SCT@counts) <- new_sct
rownames(tspe@assays$integrated@data) <- new_integrated

# Rename Spatial assay to RNA and set as default as is needed for input SPOTlight
tspe_rename <- RenameAssays(object = tspe, Spatial = 'RNA')
DefaultAssay(object = tspe_rename) <- "RNA"

################
## analysis ####
################

## downsample cells ####
# in cases when we have more transcriptionally similar cell
# identities we need to increase our N to capture the biological heterogeneity
# between them; go with 100 as suggested in vignette?
set.seed(20)
n_cells <- 10

idx <- split(seq(ncol(thymus)), thymus$subset_annotations)
cells.keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells) n_cells <- n
  sample(i, n_cells)
})
thymus.keep <- thymus[, unlist(cells.keep)]

## Train non-negative matrix factorization model ####
#slices <- c(S1_A1, S1_B1.1, S1_C1.2, S1_D1.3, S2_B1.5, S2_C1.6, S2_D1.7)
slices <- c("S1_A1", "S1_B1", "S1_C1", "S1_D1", "S2_B1", "S2_C1", "S2_D1")

# Train NMF on the thymus single cell data;
# learned model applied to first slice, but can then be used for deconvolution
# of other slices as well 
tspe.tmp <- subset(tspe_rename, subset=orig.ident=="S2_A1")
res.modls <- trainNMF(
  x = thymus.keep,
  y = tspe.tmp,
  groups = as.character(thymus.keep$subset_annotations),
  mgs = markers.subsets,
  hvg = unique(hvgs.subsets$Features),
  weight_id = "avg_log2FC",
  group_id = "cell_annot",
  gene_id = "gene")


# Visualise learned topic profiles
## how specific each topic signature is for each cell identity.
plotTopicProfiles(
  x = res.modls$mod,
  y = as.character(thymus.keep$subset_annotations),
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1) +
  theme(aspect.ratio = 1)

## how consistent is signature for all the cells from the same cell identity.
plotTopicProfiles(
  x = res.modls$mod,
  y = as.character(thymus.keep$subset_annotations),
  facet = TRUE,
  min_prop = 0.01,
  ncol = 10) 

## sanity check learned topics
sign <- basis(modls)
colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
head(sign)


#
#res.tmp <- lapply(slices, analyzeSPOTS, spe=tspe_rename, sce=thymus.keep,
#                  cell_annot=as.character(thymus.keep$subset_annotations),
#                  markers=markers.subsets,
#                  variable_features = unique(subset_hvgs$Features))


## Deconvolute thymus slices by learned ####
# for one slice
res.tmp <- runDeconvolution(
  x = tspe.tmp,
  mod = res.modls[["mod"]],
  ref = res.modls[["topic"]])

# figure out how to run across slices
#res.deconv <- lapply(slices, analyzeSPOTS, spe=tspe_rename, modls=res.modls)

## Visualise deconvolution results ####
thymus.mat <- res.tmp$mat

## correlation of cell types found in same spot ####
plotCorrelationMatrix(thymus.mat)
plotInteractions(thymus.mat, which = "heatmap", metric = "prop")
plotInteractions(thymus.mat, which = "network")

## pie scatter on histology image ####
ct <- colnames(thymus.mat)
mat_nn <- thymus.mat
mat_nn[mat_nn < 0.1] <- 0

#  'paletteMartin' from the 'colorBlindness' package)
paletteMartin <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
           "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
           "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")
           
pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct

p_scatter <- plotSpatialScatterpie(
  x = tspe.tmp,
  y = thymus.mat,
  slice = "S2_A1.4",
  cell_types = colnames(thymus.mat),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.4) +
  scale_fill_manual(
    values = pal,
    breaks = names(pal)) +
  coord_flip() +
  scale_x_reverse() +
  theme(aspect.ratio = 1)

p_histo <- SpatialDimPlot(tspe,
               images="S2_A1.4", label = FALSE, 
               pt.size.factor = 1) 

p_scatter / p_histo


## Annotate across slices and clusters: WIP ####
res.annotated <- res.tmp %>% bind_rows

res.annotated %>%
  group_by(location) %>%
  mutate(spots_per_region = n()) %>%
  group_by(location, cluster, slice) %>%
  summarize(proportion_total = sum(proportion, na.rm = TRUE),
            spots_per_region=unique(spots_per_region)) %>%
  mutate(total_proportion=proportion_total/spots_per_region) %>%
  
  
  ggplot( aes(x=cluster, y=total_proportion, color=location)) +
  #ggplot( aes(x=location, y=total_proportion, color=cluster)) +
  geom_boxplot() +
  #scale_color_manual(values=colvalues) +
  cowplot::theme_minimal_hgrid() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

## other visialisations: WIP ####
p1 <- DimPlot(tspe, reduction = "umap", label = FALSE)

slices.sub <- c("S2_A1.4", "S1_B1.1", "S1_C1.2", "S1_D1.3",
                "S2_B1.5", "S2_C1.6", "S2_D1.7")

spatialplots <- lapply(slices.sub, function(slice) {
  SpatialDimPlot(tspe,
                 images=slice, label = FALSE, 
                 pt.size.factor = 1) +
    theme(legend.position="none")
})

cowplot::plot_grid(plotlist = spatialplots)

p2 <- SpatialDimPlot(tspe,
                     images=c("S2_A1.4", "S1_B1.1", "S1_C1.2", "S1_D1.3",
                              "S2_B1.5", "S2_C1.6"), label = FALSE, 
                     pt.size.factor = 1) +
  theme(legend.position="none")

SpatialDimPlot(tspe,
               images="S1_A1.7", label = FALSE, 
               pt.size.factor = 1) 

p2 / p3





