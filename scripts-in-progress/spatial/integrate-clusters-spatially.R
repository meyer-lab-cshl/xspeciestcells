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
library(ggh4x)

#################
## functions ####
#################


extractSlice <- function(slice, spe) {
  tmp <- subset(spe, subset=orig.ident==slice)
}

analyzeSPOTS <- function(slice, spe.list, modls=NULL, sce=NULL, cell_annot=NULL,
                         markers=NULL, variable_features=NULL) { 
  cat(slice)
  tmp <- spe.list[[names(spe.list) == slice]]
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

visualiseTopics <- function(slice, res.slices) {
  mod <- res.slices[[which(names(res.slices) == slice)]]$NMF
  ## how specific each topic signature is for each cell identity.
  p_topics <- plotTopicProfiles(
    x = mod,
    y = as.character(thymus.keep$subset_annotations),
    facet = FALSE,
    min_prop = 0.01,
    ncol = 1) +
    theme(aspect.ratio = 1)
  ggsave(paste0("results/plots/thymus.batchC.", slice, ".topics.pdf"),
         plot=p_topics, width=8, height=8)
  
  ## how consistent is signature for all the cells from the same cell identity.
  p_topics_cell <- plotTopicProfiles(
    x = mod,
    y = as.character(thymus.keep$subset_annotations),
    facet = TRUE,
    min_prop = 0.01,
    ncol = 10) 
  ggsave(paste0("results/plots/thymus.batchC.", slice, ".topics.percell.pdf"),
         plot=p_topics_cell, width=15, height=20)
}

visualiseSpatial <- function(slice, slices.mapping, res.slices, tspe.slices,
                             palette) {
  cat(slice)
  ## extract deconvolution results ####
  slice.version <- slices.mapping[names(slices.mapping) == slice]
  thymus.mat <- res.slices[[which(names(res.slices) == slice)]]$mat
  thymus.slice <- tspe.slices[[which(names(tspe.slices) == slice)]]
  
  ## correlation of cell types found in same spot ####
  p_corr <- plotCorrelationMatrix(thymus.mat, insig="blank")
  ggsave(plot=p_corr,
         filename = paste0('results/plots/corr_matrix_', slice, '.pdf'),
         width = 7, height=7)
  
  ## pie scatter on histology image ####
  ct <- colnames(thymus.mat)
  mat_nn <- thymus.mat
  mat_nn[mat_nn < 0.1] <- 0
  
  pal <- colorRampPalette(paletteMartin)(length(ct))
  names(pal) <- ct
  
  p_scatter <- plotSpatialScatterpie(
    x = thymus.slice,
    y = thymus.mat,
    slice = slice.version,
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
  
  p_histo <- SpatialDimPlot(thymus.slice,
                            images=slice.version, label = FALSE, 
                            pt.size.factor = 1) 
  
  p_combined <- p_scatter / p_histo
  ggsave(plot=p_combined,
         filename = paste0('results/plots/spatial_plot_', slice, '.pdf'),
         width = 12, height=12)
}

annotateSPOTS <- function(slice, res.slices, tspe.slices,
                          cluster_startswith="thy") { 
  ## extract deconvolution results ####
  slice.version <- slices.mapping[names(slices.mapping) == slice]
  thymus.mat <- res.slices[[which(names(res.slices) == slice)]]$mat
  thymus.slice <- tspe.slices[[which(names(tspe.slices) == slice)]]
  
  spot_annotation <- tibble(spotid=names(thymus.slice$Location),
                            location=thymus.slice$Location)
  
  cell_location <- thymus.mat %>%
    as.data.frame %>%
    rownames_to_column("spotid") %>%
    left_join(spot_annotation) %>%
    select(spotid, location, everything()) %>%
    pivot_longer(starts_with(cluster_startswith),
                 names_to = "cluster", values_to="proportion") %>%
    mutate(slice=slice)
}


############
## data ####
############
colvalues <- colorRampPalette(brewer.pal(12,"Paired"))(18)
#  'paletteMartin' from the 'colorBlindness' package)
paletteMartin <- c("#000000", "#004949", "#009292", "#ff6db6", "#ffb6db",
                            "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
                            "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")
                            
colhisto <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d')


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

# Available image data
slices <- c("S2_A1", "S1_B1", "S1_C1", "S1_D1", "S2_B1", "S2_C1", "S2_D1")
slices.mapping <- c("S2_A1.4", "S1_B1.1", "S1_C1.2", "S1_D1.3", "S2_B1.5",
                    "S2_C1.6", "S2_D1.7")
names(slices.mapping) <- slices

################
## analysis ####
################

## Overview spatial transcriptomics data ####
p_dim <- DimPlot(tspe, reduction = "umap", label = FALSE) +
  scale_color_manual(values=colhisto)
legend_dim <- cowplot::get_legend(p_dim)

p_dim <- p_dim  +
  theme(legend.position="none",
        axis.text =element_text(size=8),
        axis.title =element_text(size=10)) 

spatialplots <- lapply(slices.mapping, function(slice) {
  SpatialDimPlot(tspe,
                 images=slice, label = FALSE, 
                 pt.size.factor = 1.3) +
    labs(title=slice) +
    theme(legend.position="none",
          title=element_text(size=8)) +
    scale_fill_manual(values=colhisto)
})

p_spatial <- cowplot::plot_grid(plotlist = append(list(p_dim=p_dim,
                                                       legend=legend_dim),
                                                  spatialplots),
                                nrow=3)

ggsave(plot=p_spatial, 'results/plots/histology_plots.pdf', width = 9, height=9)

## downsample cells ####
# in cases when we have more transcriptionally similar cell
# identities we need to increase our N to capture the biological heterogeneity
# between them; prioritize picking from one batch
with(thymus@meta.data, table(Batch, subset_annotations))
# -> batch C best?
thymus.batchC <-  subset(thymus, subset=Batch=="C") 
print(thymus.batchC) # 7,379 cells

#go with 100 as suggested in vignette?
set.seed(20)
n_cells <- 100

idx <- split(seq(ncol(thymus.batchC)), thymus.batchC$subset_annotations)
cells.keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells) n_cells <- n
  sample(i, n_cells)
})
thymus.keep <- thymus.batchC[, unlist(cells.keep)]

## Train non-negative matrix factorization model ####

# Train NMF on the thymus single cell data;
# learned model applied to first slice, but can then be used for deconvolution
# of other slices as well 

tspe.slices <- lapply(slices, extractSlice, tspe_rename)
names(tspe.slices) <- slices

res.decon <- parallel::mclapply(slices, analyzeSPOTS, spe.list=tspe.slices,
                                sce=thymus.keep,
                                cell_annot=as.character(thymus.keep$subset_annotations),
                                markers=markers.subsets,
                                variable_features = unique(subset_hvgs$Features),
                                mc.cores=8)
names(res.decon) <- slices
#saveRDS(res.decon, 'results/NMF.modls.thymus.batchC.allslices.rds')


## Visualise learned topic profiles ####
res.topics <- parallel::mclapply(slices, visualiseTopics, res.slices=res.decon,
                              mc.cores=8)


## sanity check learned topics
#sign <- basis(modls)
#colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
#head(sign)

## Visualise spatial distribution ####
slices_tmp <- slices[-4]
res.spatial <- parallel::mclapply(slices_tmp, visualiseSpatial,
                                 slices.mapping=slices.mapping,
                                 res.slices=res.decon, tspe.slices=tspe.slices,
                                 mc.cores=8)

## Annotate across slices and clusters ####
res.annotated <-  parallel::mclapply(slices, annotateSPOTS,
                                     res.slices=res.decon,
                                     tspe.slices=tspe.slices) %>%
  bind_rows

## proportion of clusters/cell states across regions ####
summary.annotated <- res.annotated %>%
  group_by(location, cluster, slice) %>%
  summarize(proportion = mean(proportion, na.rm = TRUE), .groups ='drop') %>%
  group_by(cluster, slice) %>%
  mutate(proportion_cluster_across_regions=proportion/sum(proportion),
         sum_pp = sum(proportion_cluster_across_regions),
         celltype=gsub("_.*", "", cluster)) %>%
  ungroup() %>%
  mutate(celltype=factor(celltype,
                         levels=c("thyCD4", "thyCD8", "thyGDT", "thyNKT",
                                  "thyMAIT"))) %>%
  arrange(celltype) %>%
  mutate(cluster=forcats::fct_inorder(cluster))

# boxplot
p_bp <- ggplot(summary.annotated,
               aes(x=location, y=proportion_cluster_across_regions,
                   color=location)) +
    facet_wrap(~cluster, nrow=5, scales='fixed')  +
    geom_boxplot() +
    scale_color_manual(values=colhisto) +
    cowplot::theme_minimal_hgrid() +
    theme(axis.text.x = element_blank())

# proportions plot
design <- rbind(c(1:7), c(8:13,NA), c(14:19,NA), c(20:26), c(27:31, NA, NA))
p_bar <- ggplot(summary.annotated,
                aes(x=slice, y=proportion_cluster_across_regions,
                    fill=location)) +
  geom_bar( stat='identity') +
  scale_fill_manual(values=colhisto) +
  facet_manual(~cluster, design) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
ggsave(plot=p_bar,
       filename = 'results/plots/proportion_cluster_across_regions.pdf',
       width = 15, height=12)

## average correlation of cell types found in same spot ####

slices.corr <- lapply(slices_tmp, function(slice) {
  mat <- res.decon[[which(names(res.decon) == slice)]]$mat
  mat <- mat[, colSums(mat) > 0]
  corr <- cor(mat)
})
slices.mean.corr <- Reduce(f='+', slices.corr)/length(slices.corr)

p.mat <- ggcorrplot::cor_pmat(x = slices.mean.corr, conf_int = 0.95,
                              method = 'pearson')
p_corr <- ggcorrplot::ggcorrplot(corr = slices.mean.corr, p.mat = p.mat,
                       hc.order = TRUE, 
                       insig = 'blank', lab = FALSE) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.text.x = element_text(angle = 60, vjust = 1), 
        axis.text = element_text(vjust = 0.5))

## visualise cell type specific proportions ####

colvalues <- colorRampPalette(brewer.pal(12,"Paired"))(18)

thymus.S2_A1 <- tspe.slices[[which(names(tspe.slices) == "S2_A1")]]
thymus.mat.S2_A1 <- res.decon[[which(names(res.decon) == "S2_A1")]]$mat

x <- thymus.S2_A1
y <- thymus.mat.S2_A1[, grepl("thyNKT", colnames(thymus.mat.S2_A1))]
y <- t(apply(y, 1, function(x) {
  if(sum(x) == 0) return(x)
  return(x/sum(x))
}
))

x <- SPOTlight:::.extract_coord(x = thymus.S2_A1, slice = "S2_A1.4", img = TRUE)
x <- SPOTlight:::.x_cnames(x)

df <- merge(x, y, by = 0, all = TRUE)
df$location <- thymus.S2_A1$Location
cell_types <- colnames(y)

p_nkt <- ggplot() + coord_fixed()
p_nkt <- p_nkt + scatterpie::geom_scatterpie(data = df,
                                aes(x = coord_x, 
                                    y = abs(coord_y)),
                                    #color=location),
                                cols = cell_types,
                                pie_scale = 0.4,
                                color = NA) + 
  scale_fill_manual(values=colvalues) +
  scale_color_manual(values=colhisto) +
  coord_flip() +
  scale_x_reverse() +
  scale_size_continuous(range=c(0.01, 0.5)) +
  theme_void() +
  theme(legend.key.size = unit(0.5, "lines"),
        aspect.ratio = 1)

p_histo <- SpatialDimPlot(thymus.S2_A1,
                          images="S2_A1.4", label = FALSE, 
                          pt.size.factor = 1.3)  +
  scale_fill_manual(values=colhisto)

p_combined <- p_nkt / p_histo
ggsave(plot=p_combined,
       filename = paste0('results/plots/spatial_plot_', slice, '.pdf'),
       width = 12, height=12)


## WIP ####
df <- merge(x, y, by = 0, all = TRUE) %>%
  as_tibble %>%
  select(-Row.names) %>%
  pivot_longer(starts_with("thyNKT"), names_to="cluster",
               values_to="proportion")

ggplot(data = df, aes(x = coord_x, y = abs(coord_y - ymax))) +
  geom_point(aes(size=proportion, color =cluster)) +
  facet_grid(~cluster) +
  theme_void() +
  theme(legend.key.size = unit(0.5, "lines")) +
  coord_flip() +
  scale_y_reverse() +
  scale_x_reverse() +
  scale_size_continuous(range=c(0.01, 0.5)) +
  theme(aspect.ratio = 1)