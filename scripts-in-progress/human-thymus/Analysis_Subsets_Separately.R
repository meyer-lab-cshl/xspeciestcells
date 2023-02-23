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
library(symphony)
library(SCpubr)
library(ggthemes)

#########################################################################################################################
filtered_seurat_Harmony <- readRDS("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/filtered_seurat_Harmony.07.22.22.RDS")

Idents(filtered_seurat_Harmony) <- "orig.ident"
remove <- c("CD1a_1_Thymus", "CD1c_1_Thymus")
filtered_seurat_Harmony <- subset(filtered_seurat_Harmony, idents=remove, invert = TRUE)
########################################################################################################################


######################################
## Reclustering of Thymic NKT cells ##
######################################
Idents(filtered_seurat_Harmony) <- "group.ident"
NKT_cells_Total <- subset(filtered_seurat_Harmony, ident = c("NKT_Thymus", "NKT_PBMC"))
Idents(NKT_cells_Total) <- "group.ident"
NKT_Thymus_samples <- subset(NKT_cells_Total, ident = c("NKT_Thymus"))
NKT_Thymus_samples <- NormalizeData(NKT_Thymus_samples)
NKT_Thymus_samples <- SCTransform(NKT_Thymus_samples, vars.to.regress = c('percent.mt'), 
                                  verbose=TRUE, method = "glmGamPoi")
NKT_Thymus_samples <- RunHarmony(NKT_Thymus_samples, group.by.vars = c("Method", "Batch"), assay.use = "SCT", 
                                 max.iter.harmony = 20)
NKT_Thymus_samples <- RunPCA(NKT_Thymus_samples)
NKT_Thymus_samples <- RunUMAP(NKT_Thymus_samples, dims = 1:15)
NKT_Thymus_samples <- FindNeighbors(NKT_Thymus_samples, dims = 1:15)
NKT_Thymus_samples <- FindClusters(NKT_Thymus_samples, resolution = 0.7)

DimPlot(NKT_Thymus_samples)
FeaturePlot(NKT_Thymus_samples, features = "RAG1")

NKT_Thymus_samples$new_clusters_NKT <- case_when(
  NKT_Thymus_samples$seurat_clusters == '3' ~ '0',
  NKT_Thymus_samples$seurat_clusters == '2' ~ '1',
  NKT_Thymus_samples$seurat_clusters == '1' ~ '2',
  NKT_Thymus_samples$seurat_clusters == '6' ~ '3',
  NKT_Thymus_samples$seurat_clusters == '5' ~ '4',
  NKT_Thymus_samples$seurat_clusters == '0' ~ '5',
  NKT_Thymus_samples$seurat_clusters == '4' ~ '6'
)

NKT_Thymus_samples$new_clusters_NKT  <- as.factor(NKT_Thymus_samples$new_clusters_NKT)
Idents(NKT_Thymus_samples) <- "new_clusters_NKT"

colors_clusters <- c("0" = "#f4c40f", "1" = "#b75347", "2" = "#d8443c", "3" = "#e09351", "4" = "#2b9b81", 
                     "5" = "#421401", "6" = "#92c051", "7" = "#9f5691", "8" = "#17154f", "9" = "#74c8c3", 
                     "10" = "#5a97c1", "11" = "gold", "12" = "#a40000", "13" = "#72bcd5", "14" = "grey50",
                     "15" = "orange", "16" = "blueviolet", "17" = "#0a2e57", "18" = "olivedrab2")

colors_clusters_NKT <- c("0" = "#d8443c", "1" = "#74c8c3", "2" = "#5a97c1", "3" = "#9f5691", "4" = "gold", 
                     "5" = "#a40000", "6" = "grey50")

DimPlot(NKT_Thymus_samples, cols = colors_clusters_NKT, pt.size = 1)

SCpubr::do_DimPlot(sample = NKT_Thymus_samples, 
                   label = FALSE, 
                   label.color = "black", colors.use = colors_clusters_NKT,
                   legend.position = "right")

NKT_Thymus_samples <- NormalizeData(NKT_Thymus_samples)
SCpubr::do_FeaturePlot(sample = NKT_Thymus_samples, 
                       features = "CCR7",
                       plot.title = "",
                       reduction = "umap",
                       viridis_color_map = "inferno")

##DEGs#################################################################################################################
NKT.clusters.markers <- FindAllMarkers(NKT_Thymus_samples, test.use = 'wilcox', 
                                   logfc.threshold = 0.3, min.pct = 0.3, only.pos = TRUE)
top5 <- NKT.clusters.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5_distinct <- base::unique(top5$gene)

write.csv(NKT.clusters.markers, "/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/clusters.markers.csv")

NKT_Thymus_samples <- ScaleData(NKT_Thymus_samples, features = top5_distinct)
#DotPlot_colors <- met.brewer("Hiroshige", type = "discrete", n = 2)
#DotPlot_colors <- c("#72bcd5", "#1E466E")
DotPlot_colors <- c("lightsteelblue1", "#1E466E")

fig_dir <- "/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/"

pdf(paste0(fig_dir, "Features_per_NKT_clusters.pdf"), width=8, height=7)
p3 <- DotPlot(NKT_Thymus_samples, features = top5_distinct, dot.scale = 8, cols = DotPlot_colors,
              col.min = -1, col.max = 2, dot.min = 0) + RotatedAxis()
print(p3)
dev.off()

###Distribution of cells per Donor####################################################################################
Idents(NKT_Thymus_samples) <- "new_clusters_NKT"
NKT_cell_numbers_thymus <- as.data.frame(NKT_Thymus_samples@meta.data) %>%
  dplyr::select(Donor, new_clusters_NKT) %>% group_by(new_clusters_NKT, Donor) %>% summarize(N = n()) %>% ungroup()
NKT_cell_numbers_thymus <- NKT_cell_numbers_thymus %>% complete(new_clusters_NKT, fill=list(N=0))

My_Theme = theme(
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  axis.title.y = element_text(size = 18))

pdf(paste0(fig_dir, "Nb_NKT_Thymus_per_donor.pdf"), width=7, height=7)
p <- ggplot(NKT_cell_numbers_thymus, aes(x = new_clusters_NKT, y = N, fill = factor(Donor))) +
  geom_bar(stat="identity") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 500)) +
  labs(x = "Clusters", y = "Nb of cells", fill = "Donors") +
  scale_fill_manual(values = c("#ffeda0", "#feb24c", "#f03b20"),
                    labels = levels(NKT_cell_numbers_thymus$Donor), 
                    drop = FALSE, limits = levels(NKT_cell_numbers_thymus$Donor)) +
  theme_classic(base_size = 14) + My_Theme
print(p)
dev.off()

NKT_cell_numbers_thymus_donor2 <- NKT_cell_numbers_thymus %>% dplyr::filter(Donor == 2)

NKT_cell_numbers_thymus_donor2_test <- NKT_cell_numbers_thymus_donor2 %>% 
  mutate(freq = N/sum(NKT_cell_numbers_thymus_donor2$N)*100) %>% ungroup()

NKT_cell_numbers_thymus_donor3 <- NKT_cell_numbers_thymus %>% dplyr::filter(Donor == 3)

NKT_cell_numbers_thymus_donor3_test <- NKT_cell_numbers_thymus_donor3 %>% 
  mutate(freq = N/sum(NKT_cell_numbers_thymus_donor3$N)*100) %>% ungroup()

NKT_cell_numbers_thymus_donor4 <- NKT_cell_numbers_thymus %>% dplyr::filter(Donor == 4)

NKT_cell_numbers_thymus_donor4_test <- NKT_cell_numbers_thymus_donor4 %>% 
  mutate(freq = N/sum(NKT_cell_numbers_thymus_donor4$N)*100) %>% ungroup()

distribution_donors_thymus <- rbind(NKT_cell_numbers_thymus_donor2_test, NKT_cell_numbers_thymus_donor3_test, 
                                    NKT_cell_numbers_thymus_donor4_test)

distribution_donors_thymus$Donor <- factor(distribution_donors_thymus$Donor)

distribution_donors_thymus <- distribution_donors_thymus %>% complete(new_clusters_NKT, fill=list(N=0, freq=0))

pdf(paste0(fig_dir, "Distribution_per_cluster_NKT_Thymus.pdf"), width=7, height=7)
p <- ggplot(data= subset(distribution_donors_thymus, !is.na(Donor)), aes(x = Donor, y = freq, 
                                                                         fill = forcats::fct_rev(factor(new_clusters_NKT)))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = colors_clusters_NKT, 
                    labels = levels(distribution_donors_thymus$new_clusters_NKT), 
                    drop = FALSE, limits = levels(distribution_donors_thymus$new_clusters_NKT),
                    guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Donors", y = "%", fill = "Clusters") + theme_classic() + My_Theme 
print(p)
dev.off()

####VISION###############################################################################################################
signatures <- c("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/MouseNKTsubsets.signatures.gmt",
                "/Volumes/Samsung_T5/Human_MAIT_NKT/Data/c5.go.bp.v2022.1.Hs.symbols.gmt")

vision.obj.seurat.Harmony <- Vision(NKT_Thymus_samples, signatures = signatures, 
                                    projection_methods = NULL, meta = NKT_Thymus_samples@meta.data)

vision.obj.seurat.Harmony <- analyze(vision.obj.seurat.Harmony)

viewResults(vision.obj.seurat.Harmony)

# Display autocorrelation coefficients, p-values for signatures
head(getSignatureAutocorrelation(vision.obj.seurat.Harmony), n = 10)

correlations <- getSignatureAutocorrelation(vision.obj.seurat.Harmony)
correlations <- rownames_to_column(correlations, var = "pathway")
correlations %>% dplyr::filter(str_detect(pathway, "Stage"))

umap <- vision.obj.seurat.Harmony@Projections$Seurat_umap

corr_response <- correlations %>% dplyr::filter(str_detect(pathway, "^GOBP_RESPONSE_TO_PHORBOL_13_ACETATE_12_MYRISTATE$"))
sigScores <- getSignatureScores(vision.obj.seurat.Harmony)[, "GOBP_RESPONSE_TO_PHORBOL_13_ACETATE_12_MYRISTATE"]
ggplot() + aes(x=umap[, 1], y=umap[, 2], color = sigScores) + 
  geom_point(size = 1) +
  xlim(-10, 12) + ylim(-7.5, 5.5) +
  scale_colour_viridis(option = "magma", limits = c(0, 1), oob = scales::squish) +
  ggtitle(glue("GOBP_RESPONSE_TO_PHORBOL_13_ACETATE_12_MYRISTATE\n C' = {round(corr_response[2], digit = 2)}, FDR = {round(corr_response[4], digit = 2)}")) +
  labs(x = "", y = "") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 0.2)) + NoLegend()

### Monocle3 trajectory####################################################################################################
cds <- SeuratWrappers::as.cell_data_set(NKT_Thymus_samples)

#to get cell metadata
colData(cds)
#to gene metadata
fData(cds)
fData(cds)$gene_short_name <- rownames(fData(cds))
fData(cds)

#to get counts
counts(cds)

## assign partitions
recreate_partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate_partition) <- cds@colData@rownames
recreate_partition <- as.factor(recreate_partition)

cds@clusters$UMAP$partitions <- recreate_partition

## assign the cluster info
Idents(NKT_Thymus_samples) <- "new_clusters_NKT"
list.clusters <- NKT_Thymus_samples@active.ident
cds@clusters$UMAP$clusters <- list.clusters

##assign UMAP coordinates
cds@int_colData@listData$reducedDims$UMAP <- NKT_Thymus_samples@reductions$umap@cell.embeddings

#plot
clusters.before.trajectory <- plot_cells(cds,
                                         color_cells_by = "cluster",
                                         label_cell_groups = TRUE,
                                         group_label_size = 5,
                                         label_branch_points = FALSE,
                                         label_roots = FALSE,
                                         label_leaves = FALSE) +
  theme(legend.position = "right")

## learn trajectory
cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by = "cluster",
           label_cell_groups = FALSE,
           group_label_size = 5,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)

##order cells per pseudotime
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[,clusters(cds) == 0]))
cds <- order_cells(cds, root_cells = NULL)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           group_label_size = 5,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           cell_size = 1, show_trajectory_graph = TRUE) + 
  viridis::scale_color_viridis(option = "D") + xlim(-10, 12) + ylim(-7.5, 5.5)

# cells ordered by monocle3 pseudotime
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(new_clusters_NKT, monocle3_pseudotime, median), 
                        fill = new_clusters_NKT)) +
  geom_boxplot() + 
  scale_fill_manual(values = colors_clusters_NKT) +
  ylab("Clusters") + xlab("Pseudotime") + theme_clean()

######TCR analysis########################################################################################################
Idents(NKT_Thymus_samples) <- "TRAV10_TRAJ18"
DimPlot(NKT_Thymus_samples, 
        pt.size = 2, 
        ncol = 1, 
        label = FALSE, label.size = 0, order = T, cols = c("grey10", "red")) +
  NoLegend() + xlim(-10, 12) + ylim(-7.5, 5.5) + ggtitle("TRAV10-TRAJ18 expressing cells")

TCR_data_NKT_cells.df <- as.data.frame(NKT_Thymus_samples@meta.data)
TCR_data_NKT_cells.df <- TCR_data_NKT_cells.df %>%
  dplyr::select(-nCount_RNA, -nFeature_RNA, -percent.mt, -cell.ident, -Sex, -Age_in_weeks, -Donor, -Batch, -Method, -nCount_SCT,
         -nFeature_SCT, -SCT_snn_res.0.7, -seurat_clusters, -orig.ident, -Tissue, -Test_cell_ident_for_cLISI,- SCT_snn_res.0.9,
         -SCT_snn_res.1, -group.ident, -TCRa_g_chain, -TCRa_g_clonotype, -TCRb_d_chain, -TCRb_d_clonotype, -TRAV10_TRAJ18, -TRAV1_TRAJ33, -TRAV1, -new_clusters) %>%
  na.omit() %>%
  dplyr::rename(TRAV = TCR_Alpha_Gamma_V_gene_Dominant, TRAJ = TCR_Alpha_Gamma_J_gene_Dominant, CDR3a = TCR_Alpha_Gamma_CDR3_Translation_Dominant,
                TRBV = TCR_Beta_Delta_V_gene_Dominant, TRBD = TCR_Beta_Delta_D_gene_Dominant, TRBJ = TCR_Beta_Delta_J_gene_Dominant,
                CDR3b = TCR_Beta_Delta_CDR3_Translation_Dominant) %>%
  mutate(CDR3b = paste0("C", CDR3b)) %>%
  mutate(CDR3a = paste0("C", CDR3a))

Total_TRAV_TRAJ_distribution <- TCR_data_NKT_cells.df %>%
  dplyr::select(-new_clusters_NKT, -CDR3a, -TRBV, -TRBD, -TRBJ, -CDR3b) %>%
  dplyr::filter(TRAV == "TRAV10*01") %>%
  mutate(TRAV = str_remove(TRAV, pattern = "\\*[^.]*$"),
         TRAJ = str_remove(TRAJ, pattern = "\\*[^.]*$")) %>%
  group_by(TRAV, TRAJ) %>%
  summarise(n = n()) %>% ungroup() %>% as.data.frame()

#determine how many TRAV and TRAJ are in the samples
V_alpha_genes <- unique(Total_TRAV_TRAJ_distribution$TRAV) %>% sort()
J_alpha_genes <- unique(Total_TRAV_TRAJ_distribution$TRAJ) %>% sort()
Valpha_color_vector <- c(V_alpha_genes, J_alpha_genes)  
Valpha_color_vector <- unique(Valpha_color_vector)

#create a named vector where each TRAV and TRAJ has a dedicated color
#["Austria",["#a40000","#16317d","#007e2f","#ffcd12","#721b3e","#b86092","#00b7a7"]]
Valpha_colors_vector <- c("TRAV10" = "#a40000", "TRAJ18" = "#16317d", "TRAJ24" = "#007e2f", "TRAJ37" = "#ffcd12")

par(cex = 1.5, mar = c(0, 0, 0, 0))
chordWrap <- function(Total_TRAV_TRAJ_distribution){
  circos.clear()
  chordDiagram(Total_TRAV_TRAJ_distribution, 
               annotationTrack = c("grid"), 
               annotationTrackHeight = c(0.1, 0.1),
               grid.col = Valpha_colors_vector,
               row.col = c("#A40000"),
               order = c(rep("TRAV10", 3), rep("TRAJ18", 1), rep("TRAJ24", 1), rep("TRAJ37", 1)),
               transparency = 0.3,
               link.lwd = 0.3, 
               link.lty = 5, 
               link.border = 1,
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(Total_TRAV_TRAJ_distribution))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA) 
}

chordWrap(Total_TRAV_TRAJ_distribution)

Total_TRAV_TRAJ <- TCR_data_NKT_cells.df %>%
  dplyr::select(-new_clusters_NKT, -TRBV, -TRBD, -TRBJ, -CDR3b) %>%
  dplyr::filter(TRAV == "TRAV10*01") %>%
  mutate(TRAV = str_remove(TRAV, pattern = "\\*[^.]*$"),
         TRAJ = str_remove(TRAJ, pattern = "\\*[^.]*$")) %>%
  mutate(CDR3_Length = str_length(CDR3a)) %>%
  dplyr::filter(TRAV == "TRAV10", TRAJ == "TRAJ18", CDR3a %in% str_remove(CDR3a, pattern = "\\*[^*]"))

pdf(paste0(fig_dir, "CDR3 size TRAV10-TRAJ18.pdf"), width=7, height=7)
p <- Total_TRAV_TRAJ %>% group_by(CDR3_Length) %>% summarise(N = n()) %>%
  ggplot(aes(x = CDR3_Length, y = N)) + geom_bar(stat="identity", width = 0.5, fill = "#08306b") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 300)) +
  scale_x_continuous(limits = c(10, 18), breaks = c(seq(from = 10, to = 18, by = 1))) +
  labs(x = "CDR3 Length (aa)", y = "Nb of Sequences") +
  theme_classic(base_size = 14) + My_Theme
print(p)
dev.off()

pdf(paste0(fig_dir, "CDR3 logo TRAV10-TRAJ18.pdf"), width = 10, height = 3)
p <- Total_TRAV_TRAJ %>%
  dplyr::filter(CDR3_Length == "14") %>% 
  dplyr::select(CDR3a) %>% ggseqlogo(seq_type = 'aa', method = 'bits') +
  theme(axis.text=element_text(size=28), axis.title=element_text(size=28, face="bold")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 28, b = 0, l = 0))) +
  theme(legend.position="none")
print(p)
dev.off()

Total_TRBV_TRBJ_distribution <- TCR_data_NKT_cells.df %>%
  dplyr::filter(TRAV == "TRAV10*01") %>%
  dplyr::select(-new_clusters_NKT, -TRAV, -TRAJ, -TRBD, -CDR3a, -CDR3b) %>%
  mutate(TRBV = str_remove(TRBV, pattern = "\\-[^.]*$"),
         TRBV = str_remove(TRBV, pattern = "\\*[^.]*$"),
         TRBJ = str_remove(TRBJ, pattern = "\\*[^.]*$")) %>%
  group_by(TRBV, TRBJ) %>%
  summarise(n = n()) %>% ungroup() %>% as.data.frame()

#determine how many TRAV and TRAJ are in the samples
V_beta_genes <- unique(Total_TRBV_TRBJ_distribution$TRBV) %>% sort()
J_beta_genes <- unique(Total_TRBV_TRBJ_distribution$TRBJ) %>% sort()
Vbeta_color_vector <- c(V_beta_genes, J_beta_genes)  
Vbeta_color_vector <- unique(Vbeta_color_vector)

#create a named vector where each TRAV and TRAJ has a dedicated color
Vbeta_colors <- c("TRBV4" = "#EAF1F4", "TRBV5" = "#93D4DB", "TRBV6" = "#99CCA7", "TRBV7" = "#57AA65", "TRBV12" = "#A7A365", 
                  "TRBV19" = "#C2947E", "TRBV21" = "#785158", "TRBV23" = "#757575", "TRBV25" = "#3D7D8F", "TRBV27" = "#956E73", 
                  "TRBJ1-1" = "#EB924A", "TRBJ1-2" = "#3EACCC", "TRBJ1-3" = "#61C2DA", "TRBJ1-4" = "#F4CBA0", 
                  "TRBJ1-5" = "#A7A7C7", 
                  "TRBJ1-6" = "#E8A1CF", 
                  "TRBJ2-1" = "#E7D9DB", "TRBJ2-2" = "#736B9D", "TRBJ2-3" = "#DEE9E9", "TRBJ2-4" = "#EF4F55", 
                  "TRBJ2-5" = "#76C76C", "TRBJ2-6" = "#1E2223", 
                  "TRBJ2-7" = "#E7D9DB")

par(cex = 1.4, mar = c(0, 0, 0, 0))
chordWrap <- function(Total_TRBV_TRBJ_distribution){
  circos.clear()
  chordDiagram(Total_TRBV_TRBJ_distribution, 
               annotationTrack = c("grid"), 
               annotationTrackHeight = c(0.05, 0.05),
               grid.col = Vbeta_colors, 
               order = c(rep("TRBV4", 2), "TRBV5", rep("TRBV6", 1), "TRBV7", "TRBV12",
                         rep("TRBV19", 3), rep("TRBV21", 1), "TRBV23", rep("TRBV25", 13), "TRBV27", rep("TRBJ1-1", 1),
                         rep("TRBJ1-2", 2), rep("TRBJ1-3", 2), "TRBJ1-4", rep("TRBJ1-5", 1), "TRBJ1-6", 
                         rep("TRBJ2-1", 1),
                         "TRBJ2-2", rep("TRBJ2-3", 4), "TRBJ2-4", rep("TRBJ2-5", 3), "TRBJ2-6", rep("TRBJ2-7", 4)),
               col = Vbeta_colors, transparency = 0.1, 
               link.lwd = 1, link.lty = 1, link.border = 1,
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(Total_TRBV_TRBJ_distribution))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.4))
  }, bg.border = NA) 
}

chordWrap(Total_TRBV_TRBJ_distribution)


###cell identity of each NKT clusters
C0 <- WhichCells(NKT_Thymus_samples, idents = "0")
C0.df <- data.frame(cell_ID = C0, Cluster_NKT = "NKT_C0")
C1 <- WhichCells(NKT_Thymus_samples, idents = "1")
C1.df <- data.frame(cell_ID = C1, Cluster_NKT = "NKT_C1")
C2 <- WhichCells(NKT_Thymus_samples, idents = "2")
C2.df <- data.frame(cell_ID = C2, Cluster_NKT = "NKT_C2")
C3 <- WhichCells(NKT_Thymus_samples, idents = "3")
C3.df <- data.frame(cell_ID = C3, Cluster_NKT = "NKT_C3")
C4 <- WhichCells(NKT_Thymus_samples, idents = "4")
C4.df <- data.frame(cell_ID = C4, Cluster_NKT = "NKT_C4")
C5 <- WhichCells(NKT_Thymus_samples, idents = "5")
C5.df <- data.frame(cell_ID = C5, Cluster_NKT = "NKT_C5")
C6 <- WhichCells(NKT_Thymus_samples, idents = "6")
C6.df <- data.frame(cell_ID = C6, Cluster_NKT = "NKT_C6")

NKT_clustering <- rbind(C0.df, C1.df, C2.df, C3.df, C4.df, C5.df, C6.df)

rownames(NKT_clustering) <- NKT_clustering$cell_ID
NKT_clustering$cell_ID <- NULL

filtered_seurat_Harmony <- AddMetaData(filtered_seurat_Harmony, NKT_clustering, col.name = "NKT_clustering")

filtered_seurat_Harmony$NKT_clustering <- as.factor(filtered_seurat_Harmony$NKT_clustering)

filtered_seurat_Harmony$NKT_clustering <- factor(filtered_seurat_Harmony$NKT_clustering, levels = c(
  "C0", "C1", "C2", "C3", "C4", "C5", "C6"))

Idents(filtered_seurat_Harmony) <- "NKT_clustering"

colors_clusters <- c("C0" = "#f4c40f", "C1" = "#b75347", "C2" = "#92c051", "C3" = "#0a2e57", "C4" = "#421401", 
                     "C5" = "#2b9b81", "C6" = "#d8443c")

DimPlot(filtered_seurat_Harmony, order = TRUE, pt.size = 1.2, cols = colors_clusters, na.value = "grey95")

###################################################################

saveRDS(filtered_seurat_Harmony, "/Volumes/Samsung_T5/Human_MAIT_NKT/Data/Thymic_NKT_102722.RDS")

