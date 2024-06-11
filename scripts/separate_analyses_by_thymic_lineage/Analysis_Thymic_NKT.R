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
library(ggrepel)
library(scico)

My_Theme = theme(
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  axis.title.y = element_text(size = 18))

colors_clusters <- c("0" = "#f4c40f", "1" = "#b75347", "2" = "#d8443c", "3" = "#e09351", "4" = "#2b9b81", 
                     "5" = "#421401", "6" = "#92c051", "7" = "#9f5691", "8" = "#17154f", "9" = "#74c8c3", 
                     "10" = "#5a97c1", "11" = "gold", "12" = "#a40000", "13" = "#72bcd5", "14" = "grey50",
                     "15" = "orange", "16" = "blueviolet", "17" = "#0a2e57", "18" = "#bdbdbd")

#################
## functions ####
#################

# Preprocessing function
preprocess <- function(seurobj, ndim=20, res=0.8, colvalues, celltype, 
                       approxPCA=FALSE){
  print(seurobj)
  cat("Cells:", unique(seurobj$group.ident), "\n")
  
  print("Normalize")
  seurobj <- NormalizeData(seurobj)
  print("HVG")
  seurobj <- FindVariableFeatures(seurobj, selection.method="vst",
                                  nfeatures=2000)
  print("Scale")
  seurobj <- ScaleData(seurobj)
  
  print("PCA")
  seurobj <- RunPCA(seurobj, npcs=50, verbose=FALSE, approx=approxPCA)
  
  print("Harmony")
  seurobj <- RunHarmony(seurobj, reduction = "pca", max.iter.harmony=30,
                        group.by.vars = c("Method", "Batch"))
  print(ElbowPlot(seurobj, ndims=50, reduction="harmony"))
  
  print("UMAP")
  seurobj <- RunUMAP(seurobj, reduction = "harmony", dims = 1:ndim,
                     n.neighbors=50, reduction.key = "UMAP50_")
  
  print("Cluster")
  seurobj <- FindNeighbors(seurobj, reduction = "harmony", dims = 1:ndim)
  seurobj <- FindClusters(seurobj, resolution = res, random.seed = 0)
  
  # Visualization sanity check
  print(DimPlot(seurobj, reduction = "umap", group.by="seurat_clusters",
                label = TRUE, repel = TRUE) +
          #scale_color_manual(values=colvalues) +
          labs(title=celltype))
  
  return(seurobj)
}

#########DATA ###########
all_cells <- readRDS("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/seurat_filtered_harmony_02_15_23.RDS")
########################

########Analysis
# filter for thymus only:
thymus <- subset(all_cells, subset=Tissue=="Thymus")
print(thymus) # 37,369 cells

# Make raw object
thymus.filt <- CreateSeuratObject(counts=thymus@assays$RNA@counts,
                                  meta.data = thymus@meta.data)

print(thymus.filt) # 37,369 cells

thymus.nkt <- preprocess(seurobj = subset(thymus.filt, subset=group.ident=="NKT_Thymus"),
                         celltype = "NKTthy", colvalues=colvalues, ndim=10,
                         res=0.4)

DimPlot(thymus.nkt, group.by="RNA_snn_res.0.4", label=T, repel=T) +
  theme(legend.position="none") + 
  # scale_color_manual(values=colvalues) +
  labs(title="Thymic NKT")

thymus.nkt$new_clusters_NKT <- case_when(
  thymus.nkt$RNA_snn_res.0.4 == '3' ~ '0',
  thymus.nkt$RNA_snn_res.0.4 == '2' ~ '1',
  thymus.nkt$RNA_snn_res.0.4 == '0' ~ '2',
  thymus.nkt$RNA_snn_res.0.4 == '5' ~ '3',
  thymus.nkt$RNA_snn_res.0.4 == '6' ~ '4',
  thymus.nkt$RNA_snn_res.0.4 == '1' ~ '5',
  thymus.nkt$RNA_snn_res.0.4 == '4' ~ '6'
)

thymus.nkt$new_clusters_NKT  <- as.factor(thymus.nkt$new_clusters_NKT)
Idents(thymus.nkt) <- "new_clusters_NKT"

colors_clusters <- c("0" = "#f4c40f", "1" = "#b75347", "2" = "#d8443c", "3" = "#e09351", "4" = "#2b9b81", 
                     "5" = "#421401", "6" = "#92c051", "7" = "#9f5691", "8" = "#17154f", "9" = "#74c8c3", 
                     "10" = "#5a97c1", "11" = "gold", "12" = "#a40000", "13" = "#72bcd5", "14" = "grey50",
                     "15" = "orange", "16" = "blueviolet", "17" = "#0a2e57", "18" = "olivedrab2")

colors_clusters_NKT <- c("0" = "#d8443c", "1" = "#e09351", "2" = "gold", "3" = "#74c8c3", "4" = "#5a97c1",
                         "5" = "#a40000", "6" = "#72bcd5", "7" = "grey50")

DimPlot(thymus.nkt, cols = colors_clusters_NKT, pt.size = 1, label = T)

plot <- SCpubr::do_DimPlot(sample = thymus.nkt, 
                        label = F, repel = T, font.size = 20,
                        label.color = "white", colors.use = colors_clusters_NKT,
                        border.color = "black", border.size = 2,
                        reduction = "umap", pt.size = 1.2)

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/NKT_analysis_clusters_SCPubR.pdf",
       plot, device = "pdf", width=7, height=7, limitsize = F)

thymus.nkt <- NormalizeData(thymus.nkt)

Cd8aa_score <- list(c("GNG4", "CD8A", "NUCB2", "LEF1", "PRKCH", "PTPN7", "SH2D1A", "ELOVL5", "TOX",
                      "AQP3", "BCL6", "ITGA4", "MYB", "RTKN2", "IKZF2", "DUSP2", "MINDY2", "HIVEP3"))

thymus.nkt <- AddModuleScore(object = thymus.nkt, features = Cd8aa_score,
                                     assay = "RNA", name = 'Cd8aa_score')

plot <- SCpubr::do_FeaturePlot(sample = thymus.nkt, 
                            features = "Cd8aa_score1", order = T,
                            plot.title = "", legend.position = "none", min.cutoff = c(-0.50), max.cutoff = c(1),
                            border.color = "black", border.size = 2,
                            reduction = "umap", pt.size = 1.2) &
scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0.1, end = 0.85, direction = -1)

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/CD8aa_signature_NKT.pdf",
       plot, device = "pdf", width=8, height=7, limitsize = F)

Egress <- list(c("KLF2", "CORO1A", "CCR7", "CXCR4", "CXCR6", "FOXO1", "CXCR3", "S1PR1", "S1PR4",
                 "S100A4", "S100A6", "EMP3"))

thymus.nkt <- AddModuleScore(object = thymus.nkt, features = Egress,
                             assay = "RNA", name = 'Egress_score')

plot <- SCpubr::do_FeaturePlot(sample = thymus.nkt, 
                               features = "Egress_score1", order = T,
                               plot.title = "", legend.position = "none", min.cutoff = c(-0.50), max.cutoff = c(0.5),
                               border.color = "black", border.size = 2,
                               reduction = "umap", pt.size = 1.2) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0.1, end = 0.85, direction = -1)

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/Egress_signature_NKT.pdf",
       plot, device = "pdf", width=8, height=7, limitsize = F)

###Distribution of cells per Donor####################################################################################
Idents(thymus.nkt) <- "new_clusters_NKT"

# Extract and summarize NKT cell numbers by donor and new_clusters_NKT
NKT_cell_numbers_thymus <- thymus.nkt@meta.data %>%
  as.data.frame()
NKT_cell_numbers_thymus <-  NKT_cell_numbers_thymus %>% select(Donor, new_clusters_NKT) %>%
  dplyr::count(new_clusters_NKT, Donor) %>%
  rename(N = n) %>%
  ungroup()

# Compute frequency of NKT cells for each donor and new_clusters_NKT
NKT_cell_numbers_thymus_test <- NKT_cell_numbers_thymus %>%
  group_by(Donor) %>%
  mutate(freq = N/sum(N)*100) %>%
  ungroup()

# Combine frequency data for all donors and new_clusters_NKT
distribution_donors_thymus <- NKT_cell_numbers_thymus_test %>%
  complete(new_clusters_NKT, Donor, fill=list(N=0, freq=0)) %>%
  mutate(Donor = factor(Donor))

plot <- ggplot(data= subset(distribution_donors_thymus, !is.na(Donor)), aes(x = Donor, y = freq, 
                                                                         fill = forcats::fct_rev(factor(new_clusters_NKT)))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = colors_clusters_NKT, 
                    labels = levels(distribution_donors_thymus$new_clusters_NKT), 
                    drop = FALSE, limits = levels(distribution_donors_thymus$new_clusters_NKT),
                    guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Donors", y = "%", fill = "Clusters") + theme_classic() + My_Theme 

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/Distribution_Thymic_NKT_per_donor.pdf",
       plot, device = "pdf", width=7, height=7, limitsize = F)

##DEGs#################################################################################################################
Idents(thymus.nkt) <- "new_clusters_NKT"
NKT.clusters.markers <- FindAllMarkers(thymus.nkt, test.use = 'wilcox', 
                                       logfc.threshold = 0.4, min.pct = 0.3, only.pos = TRUE)
top5 <- NKT.clusters.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5_distinct <- base::unique(top5$gene)

write.csv(NKT.clusters.markers, "/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/clusters.markers.NKT.csv")

NKT_Thymus_samples <- ScaleData(thymus.nkt, features = top5_distinct)
DotPlot_colors <- c("lightsteelblue1", "#1E466E")

plot <- DotPlot(NKT_Thymus_samples, features = top5_distinct, dot.scale = 8, cols = DotPlot_colors,
              col.min = -1, col.max = 2, dot.min = 0) + RotatedAxis()
ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/Features_per_NKT_clusters.pdf",
       plot, device = "pdf", width=10, height=5, limitsize = F)

##################

SCpubr::do_FeaturePlot(sample = thymus.nkt, 
                       features = c("ZBTB16"), order = T,
                       plot.title = "", legend.position = "none", min.cutoff = c(0), 
                       max.cutoff = c(3), border.color = "grey", border.size = 1.5,
                       reduction = "umap", pt.size = 10, raster = T, raster.dpi = 2048,
                       viridis_color_map = "inferno")

VanGogh_colors <- met.brewer("VanGogh2", type = "discrete", n = 3)

plot <- SCpubr::do_FeaturePlot(sample = thymus.nkt, 
                            features = c("ZBTB16", "CCR9", "CCR7", "CD4", "FOS", "KLRB1"), order = T,
                            plot.title = "", legend.position = "none", min.cutoff = c(0, 0, 0, 0, 0, 0), max.cutoff = c(2, 2, 2, 2, 3, 2),
                            border.color = "black", border.size = 1.2,
                            reduction = "umap", pt.size = 0.6) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0.1, end = 0.85, direction = -1)

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/Genes1_NKT.pdf", 
       plot, device = "pdf", width=10, height=18)

#&
#scale_colour_gradientn(colours = alpha(rev(VanGogh_colors), 0.7))
scico_palette_names()

plot <- SCpubr::do_FeaturePlot(sample = thymus.nkt, 
                       features = c("EOMES", "GZMK", "CD8A", "CCR6", "FOS", "KLRB1"), order = T,
                       plot.title = "", legend.position = "none", min.cutoff = c(0, 0, 0, 0, 0, 0), max.cutoff = c(2, 2, 2, 2, 3, 2),
                       border.color = "black", border.size = 2,
                       reduction = "umap", pt.size = 12, raster = T, raster.dpi = 2048) &
  scale_colour_scico(palette = "lajolla", alpha = 0.8, begin = 0.1, end = 0.9)


plot <- SCpubr::do_FeaturePlot(sample = thymus.nkt, 
                       features = c("EOMES", "GZMK", "CD8A", "CCR6", "FOS", "KLRB1"), order = T,
                       plot.title = "", legend.position = "top", min.cutoff = c(0, 0, 0, 0, 0, 0), max.cutoff = c(2, 2, 2, 2, 3, 2),
                       border.color = "black", border.size = 1.2,
                       reduction = "umap", pt.size = 0.6) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0.1, end = 0.85, direction = -1)

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/Genes2_NKT.pdf", 
       plot, device = "pdf", width=10, height=18)




### Monocle3 trajectory####################################################################################################
cds <- SeuratWrappers::as.cell_data_set(thymus.nkt)

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
Idents(thymus.nkt) <- "new_clusters_NKT"
list.clusters <- thymus.nkt@active.ident
cds@clusters$UMAP$clusters <- list.clusters

##assign UMAP coordinates
cds@int_colData@listData$reducedDims$UMAP <- thymus.nkt@reductions$umap@cell.embeddings

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
plot <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           group_label_size = 5,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           cell_size = 1, show_trajectory_graph = TRUE) + 
  viridis::scale_color_viridis(option = "D") + xlim(-12, 10) + ylim(-7.5, 7.5)

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/Trajectory_NKT.pdf", 
       plot, device = "pdf", width=7, height=7)

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

TCR_data_NKT_cells.df <- as.data.frame(thymus.nkt@meta.data)
TCR_data_NKT_cells.df <- TCR_data_NKT_cells.df %>%
  dplyr::select(-nCount_RNA, -nFeature_RNA, -percent.mt, -cell.ident, -Sex, -Age_in_weeks, -Donor, -Batch, -Method,
                -seurat_clusters, -orig.ident, -Tissue, -Test_cell_ident_for_cLISI,-RNA_snn_res.0.7, -RNA_snn_res.1.2, -Egress_score1,
                -Effectorness1, -Naiveness1, -Type_3_score1, -new_clusters_id, -RNA_snn_res.0.4, -Cd8aa_score1,
                -group.ident, -TCRa_g_chain, -TCRa_g_clonotype, -TCRb_d_chain, -TCRb_d_clonotype, -TRAV10_TRAJ18, -TRAV1_TRAJ33, -TRAV1, -new_clusters) %>%
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

table(TCR_data_NKT_cells.df$TRAV)
Trav10_cells <- TCR_data_NKT_cells.df %>% dplyr::filter(TRAV == "TRAV10*01") %>% row.names()
length(Trav10_cells)

plot <- SCpubr::do_DimPlot(sample = NKT_Thymus_samples, 
                        label = FALSE,
                        cells.highlight = Trav10_cells,
                        plot.title = "", legend.position = "right",
                        na.value = "grey90", 
                        sizes.highlight = 2, colors.use = "#a40000")


ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/TRAV10_cells.pdf", 
       plot, device = "pdf", width=8, height=7)

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
               link.lwd = 0.1, link.lty = 0.1, link.border = 1,
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(Total_TRBV_TRBJ_distribution))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.4))
  }, bg.border = NA) 
}

chordWrap(Total_TRBV_TRBJ_distribution)

########################
TCR_data_NKT_cells_TCRb_diversity_thymus.df <- TCR_data_NKT_cells.df %>%
  dplyr::filter(TRAV == "TRAV10*01") %>%
  mutate(TRBV = str_remove(TRBV, pattern = "\\-[^.]*$"),
         TRBV = str_remove(TRBV, pattern = "\\*[^.]*$"),
         TRBJ = str_remove(TRBJ, pattern = "\\*[^.]*$")) %>%
  dplyr::select(-TRAV, -TRAJ, -TRBD, -CDR3a)

Vbeta_distribution_per_cluster_thymus <- TCR_data_NKT_cells_TCRb_diversity_thymus.df %>%
  group_by(TRBV, new_clusters_NKT) %>% 
  summarize(N = n(), .groups = "drop") %>% complete(new_clusters_NKT, fill=list(N=0))

distribution_Vb_per_cluster_thymus <- Vbeta_distribution_per_cluster_thymus %>%
  group_split(new_clusters_NKT) %>%
  map_dfr(~ .x %>% mutate(freq = N/sum(N)*100))

distribution_Vb_per_cluster_thymus$new_clusters_NKT <- factor(distribution_Vb_per_cluster_thymus$new_clusters_NKT)
distribution_Vb_per_cluster_thymus <- distribution_Vb_per_cluster_thymus %>% 
  complete(new_clusters_NKT, fill=list(N=0, freq=0))

distribution_Vb_per_cluster_thymus$TRBV <- factor(distribution_Vb_per_cluster_thymus$TRBV, levels = c("TRBV4",
                                                                                                      "TRBV5",
                                                                                                      "TRBV6",
                                                                                                      "TRBV7", "TRBV12",
                                                                                                      "TRBV19", "TRBV21",
                                                                                                      "TRBV23", "TRBV25",
                                                                                                      "TRBV27"))

Vbeta_colors <- c("TRBV4" = "#EAF1F4", "TRBV5" = "#93D4DB", "TRBV6" = "#99CCA7", "TRBV7" = "#57AA65", "TRBV12" = "#A7A365", 
                  "TRBV19" = "#C2947E", "TRBV21" = "#785158", "TRBV23" = "#757575", "TRBV25" = "#3D7D8F", "TRBV27" = "#956E73")

plot <- ggplot(data= subset(distribution_Vb_per_cluster_thymus), aes(x = new_clusters_NKT, y = freq, 
                                                                  fill = forcats::fct_rev(factor(TRBV)))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = Vbeta_colors, 
                    labels = levels(distribution_Vb_per_cluster_thymus$TRBV), 
                    drop = FALSE, limits = levels(distribution_Vb_per_cluster_thymus$TRBV),
                    guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Clusters", y = "%", fill = "Clusters") + theme_classic() + My_Theme

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/Distribution_Vb_per_cluster_NKT_Thymus.pdf", 
       plot, device = "pdf", width=7, height=7)


##############
TCR_data_NKT_cells_TCRb_diversity_thymus_naive.df <- TCR_data_NKT_cells_TCRb_diversity_thymus.df %>%
  dplyr::filter(new_clusters_NKT %in% c("1", "2", "3")) %>%
  mutate(TCRb_d_clonotype = paste(TRBV, TRBJ, CDR3b, sep = "_")) %>%
  dplyr::select(-new_clusters_NKT) %>% group_by(TCRb_d_clonotype) %>%
  summarise(n = n()) %>% arrange(desc(n)) %>% pivot_wider(names_from = TCRb_d_clonotype, values_from = n)

shannon_iteration <- function(){
  TCR_data_NKT_cells_TCRb_diversity_thymus_naive.df %>% 
    rrarefy(sample = 23) %>% 
    diversity(index = "shannon")
}

mean_shannon_thymus_naive <- mean(replicate(100, shannon_iteration()))

TCR_data_NKT_cells_TCRb_diversity_thymus_C5.df <- TCR_data_NKT_cells_TCRb_diversity_thymus.df %>%
  dplyr::filter(new_clusters_NKT %in% c("5")) %>%
  mutate(TCRb_d_clonotype = paste(TRBV, TRBJ, CDR3b, sep = "_")) %>%
  dplyr::select(-new_clusters_NKT) %>% group_by(TCRb_d_clonotype) %>%
  summarise(n = n()) %>% arrange(desc(n)) %>% pivot_wider(names_from = TCRb_d_clonotype, values_from = n)

shannon_iteration <- function(){
  TCR_data_NKT_cells_TCRb_diversity_thymus_C5.df %>% 
    rrarefy(sample = 23) %>% 
    diversity(index = "shannon")
}

mean_shannon_thymus_C5 <- mean(replicate(100, shannon_iteration()))

TCR_data_NKT_cells_TCRb_diversity_thymus_C6.df <- TCR_data_NKT_cells_TCRb_diversity_thymus.df %>%
  dplyr::filter(new_clusters_NKT %in% c("6")) %>%
  mutate(TCRb_d_clonotype = paste(TRBV, TRBJ, CDR3b, sep = "_")) %>%
  dplyr::select(-new_clusters_NKT) %>% group_by(TCRb_d_clonotype) %>%
  summarise(n = n()) %>% arrange(desc(n)) %>% pivot_wider(names_from = TCRb_d_clonotype, values_from = n)

shannon_iteration <- function(){
  TCR_data_NKT_cells_TCRb_diversity_thymus_C6.df %>% 
    rrarefy(sample = 23) %>% 
    diversity(index = "shannon")
}

mean_shannon_thymus_C6 <- mean(replicate(100, shannon_iteration()))

shannon_diversity_NKT <- c(mean_shannon_thymus_naive, mean_shannon_thymus_C5, mean_shannon_thymus_C6)

samples <- c("Clusters 1,2,3", "Cluster 5", "Cluster 6")
df <- data.frame(samples, shannon_diversity_NKT)
df$samples <- factor(df$samples, levels = c("Clusters 1,2,3", "Cluster 5", "Cluster 6"))

plot <- ggplot(df, aes(x = samples, y = shannon_diversity_NKT)) + 
  geom_bar(stat="identity", width = 0.2, fill = "#08306b") +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(3, 3.2)) +
  labs(x = "Samples", y = "Shannon Index") +
  theme_classic(base_size = 14) + My_Theme

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/Shannon_per_cluster_NKT_Thymus.pdf", 
       plot, device = "pdf", width=7, height=7)

####VISION###############################################################################################################
signatures <- c("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/MouseNKTsubsets.signatures.gmt")

vision.obj.seurat.Harmony <- Vision(thymus.nkt, signatures = signatures, 
                                    projection_methods = NULL, meta = thymus.nkt@meta.data)

vision.obj.seurat.Harmony <- analyze(vision.obj.seurat.Harmony)

viewResults(vision.obj.seurat.Harmony)

# Display autocorrelation coefficients, p-values for signatures
head(getSignatureAutocorrelation(vision.obj.seurat.Harmony), n = 10)

correlations <- getSignatureAutocorrelation(vision.obj.seurat.Harmony)
correlations <- rownames_to_column(correlations, var = "pathway")
correlations %>% dplyr::filter(str_detect(pathway, "Stage"))

umap <- vision.obj.seurat.Harmony@Projections$Seurat_umap

corr_response <- correlations %>% dplyr::filter(str_detect(pathway, "^Stage0_signature$"))
sigScores <- getSignatureScores(vision.obj.seurat.Harmony)[, "Stage0_signature"]
plot <- ggplot() + aes(x=umap[, 1], y=umap[, 2], color = sigScores) + 
  geom_point(size = 2) +
  xlim(-12, 11) + ylim(-8, 6) +
  scale_colour_viridis(option = "magma", limits = c(0, 1), oob = scales::squish) +
  ggtitle(glue("Stage0_signature\n C' = {round(corr_response[2], digit = 2)}, FDR = {round(corr_response[4], digit = 2)}")) +
  labs(x = "", y = "") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 0.2)) + NoLegend()

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/Stage0_signature.pdf", 
       plot, device = "pdf", width=7, height=7)

corr_response <- correlations %>% dplyr::filter(str_detect(pathway, "^NKT2_signature$"))
sigScores <- getSignatureScores(vision.obj.seurat.Harmony)[, "NKT2_signature"]
plot <- ggplot() + aes(x=umap[, 1], y=umap[, 2], color = sigScores) + 
  geom_point(size = 2) +
  xlim(-12, 11) + ylim(-8, 6) +
  scale_colour_viridis(option = "magma", limits = c(0, 1), oob = scales::squish) +
  ggtitle(glue("NKT2_signature\n C' = {round(corr_response[2], digit = 2)}, FDR = {round(corr_response[4], digit = 2)}")) +
  labs(x = "", y = "") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 0.2)) + NoLegend()

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/NKT2_signature.pdf", 
       plot, device = "pdf", width=7, height=7)


corr_response <- correlations %>% dplyr::filter(str_detect(pathway, "^NKTp_signature$"))
sigScores <- getSignatureScores(vision.obj.seurat.Harmony)[, "NKTp_signature"]
plot <- ggplot() + aes(x=umap[, 1], y=umap[, 2], color = sigScores) + 
  geom_point(size = 2) +
  xlim(-12, 11) + ylim(-8, 6) +
  scale_colour_viridis(option = "magma", limits = c(0, 1), oob = scales::squish) +
  ggtitle(glue("NKTp_signature\n C' = {round(corr_response[2], digit = 2)}, FDR = {round(corr_response[4], digit = 2)}")) +
  labs(x = "", y = "") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 0.2)) + NoLegend()

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/NKTp_signature.pdf", 
       plot, device = "pdf", width=7, height=7)


corr_response <- correlations %>% dplyr::filter(str_detect(pathway, "^NKT17_signature$"))
sigScores <- getSignatureScores(vision.obj.seurat.Harmony)[, "NKT17_signature"]
plot <- ggplot() + aes(x=umap[, 1], y=umap[, 2], color = sigScores) + 
  geom_point(size = 2) +
  xlim(-12, 11) + ylim(-8, 6) +
  scale_colour_viridis(option = "magma", limits = c(0, 1), oob = scales::squish) +
  ggtitle(glue("NKT17_signature\n C' = {round(corr_response[2], digit = 2)}, FDR = {round(corr_response[4], digit = 2)}")) +
  labs(x = "", y = "") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 0.2)) + NoLegend()

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/NKT17_signature.pdf", 
       plot, device = "pdf", width=7, height=7)

corr_response <- correlations %>% dplyr::filter(str_detect(pathway, "^NKT1c1_signature$"))
sigScores <- getSignatureScores(vision.obj.seurat.Harmony)[, "NKT1c1_signature"]
plot <- ggplot() + aes(x=umap[, 1], y=umap[, 2], color = sigScores) + 
  geom_point(size = 2) +
  xlim(-12, 11) + ylim(-8, 6) +
  scale_colour_viridis(option = "magma", limits = c(0, 1), oob = scales::squish) +
  ggtitle(glue("NKT1c1_signature\n C' = {round(corr_response[2], digit = 2)}, FDR = {round(corr_response[4], digit = 2)}")) +
  labs(x = "", y = "") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 0.2)) + NoLegend()

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/NKT1c1_signature.pdf", 
       plot, device = "pdf", width=7, height=7)

corr_response <- correlations %>% dplyr::filter(str_detect(pathway, "^NKT1c2_signature$"))
sigScores <- getSignatureScores(vision.obj.seurat.Harmony)[, "NKT1c2_signature"]
plot <- ggplot() + aes(x=umap[, 1], y=umap[, 2], color = sigScores) + 
  geom_point(size = 2) +
  xlim(-12, 11) + ylim(-8, 6) +
  scale_colour_viridis(option = "magma", limits = c(0, 1), oob = scales::squish) +
  ggtitle(glue("NKT1c2_signature\n C' = {round(corr_response[2], digit = 2)}, FDR = {round(corr_response[4], digit = 2)}")) +
  labs(x = "", y = "") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 0.2)) + NoLegend()

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/NKT1c2_signature.pdf", 
       plot, device = "pdf", width=7, height=7)

saveRDS(thymus.nkt, "/Volumes/Samsung_T5/Human_MAIT_NKT/Data/thymus.nkt_02_28_23.RDS")

thymus.nkt <- readRDS("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Data/thymus.nkt_02_28_23.RDS")


######################cNMF
counts <- GetAssayData(object = thymus.nkt, slot = "counts")
str(counts)
# Save it as a tabulated table
write.table(t(as.data.frame(counts)), file = "/Volumes/Samsung_T5/Human_MAIT_NKT/Data/counts_NKT_thymus.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE)

##########Load cNMF results
GEP_scores <- read_delim("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/cnmf_NKT_thymus/NKT_thymus_cNMF/NKT_thymus_cNMF.usages.k_7.dt_0_02.consensus.txt", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE) %>% as.data.frame()
rownames(GEP_scores) <- GEP_scores$...1
GEP_scores$...1 <- NULL

str(GEP_scores)

# Normalize each row to sum of 1
GEP_scores_norm <- t(apply(GEP_scores, 1, function(x) x/sum(x)))

thymus.nkt <- AddMetaData(thymus.nkt, GEP_scores$`1`, col.name = "GEP_Scores_1")
thymus.nkt <- AddMetaData(thymus.nkt, GEP_scores$`2`, col.name = "GEP_Scores_2")
thymus.nkt <- AddMetaData(thymus.nkt, GEP_scores$`3`, col.name = "GEP_Scores_3")
thymus.nkt <- AddMetaData(thymus.nkt, GEP_scores$`4`, col.name = "GEP_Scores_4")
thymus.nkt <- AddMetaData(thymus.nkt, GEP_scores$`5`, col.name = "GEP_Scores_5")
thymus.nkt <- AddMetaData(thymus.nkt, GEP_scores$`6`, col.name = "GEP_Scores_6")
thymus.nkt <- AddMetaData(thymus.nkt, GEP_scores$`7`, col.name = "GEP_Scores_7")

SCpubr::do_FeaturePlot(sample = thymus.nkt, 
                      features = "GEP_Scores_5", order = T,
                      plot.title = "", legend.position = "right",
                      viridis_color_map = "inferno")

Genes_scores <- read_delim("/Volumes/Samsung_T5/Human_MAIT_NKT/Data/cnmf_NKT_thymus/NKT_thymus_cNMF/NKT_thymus_cNMF.gene_spectra_score.k_7.dt_0_02.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE) %>% as.data.frame()

Genes_scores_GEP_1 <- t(Genes_scores[1,2:17205]) %>% as.data.frame()
Genes_scores_GEP_1$V1 <- as.numeric(Genes_scores_GEP_1$V1)
Genes_scores_GEP_1 <-  Genes_scores_GEP_1 %>% mutate(score = Genes_scores_GEP_1$V1, .keep = "used") %>% drop_na() %>%
  arrange(desc(score))
Genes_scores_GEP_1 <- Genes_scores_GEP_1 %>% mutate(genes = rownames(Genes_scores_GEP_1))

Genes_scores_GEP_2 <- t(Genes_scores[2,2:17205]) %>% as.data.frame()
Genes_scores_GEP_2$V1 <- as.numeric(Genes_scores_GEP_2$V1)
Genes_scores_GEP_2 <-  Genes_scores_GEP_2 %>% mutate(score = Genes_scores_GEP_2$V1, .keep = "used") %>% drop_na() %>%
  arrange(desc(score))
Genes_scores_GEP_2 <- Genes_scores_GEP_2 %>% mutate(genes = rownames(Genes_scores_GEP_2))

Genes_scores_GEP_3 <- t(Genes_scores[3,2:17205]) %>% as.data.frame()
Genes_scores_GEP_3$V1 <- as.numeric(Genes_scores_GEP_3$V1)
Genes_scores_GEP_3 <-  Genes_scores_GEP_3 %>% mutate(score = Genes_scores_GEP_3$V1, .keep = "used") %>% drop_na() %>%
  arrange(desc(score))
Genes_scores_GEP_3 <- Genes_scores_GEP_3 %>% mutate(genes = rownames(Genes_scores_GEP_3))

Genes_scores_GEP_4 <- t(Genes_scores[4,2:17205]) %>% as.data.frame()
Genes_scores_GEP_4$V1 <- as.numeric(Genes_scores_GEP_4$V1)
Genes_scores_GEP_4 <-  Genes_scores_GEP_4 %>% mutate(score = Genes_scores_GEP_4$V1, .keep = "used") %>% drop_na() %>%
  arrange(desc(score))
Genes_scores_GEP_4 <- Genes_scores_GEP_4 %>% mutate(genes = rownames(Genes_scores_GEP_4))

Genes_scores_GEP_5 <- t(Genes_scores[5,2:17205]) %>% as.data.frame()
Genes_scores_GEP_5$V1 <- as.numeric(Genes_scores_GEP_5$V1)
Genes_scores_GEP_5 <-  Genes_scores_GEP_5 %>% mutate(score = Genes_scores_GEP_5$V1, .keep = "used") %>% drop_na() %>%
  arrange(desc(score))
Genes_scores_GEP_5 <- Genes_scores_GEP_5 %>% mutate(genes = rownames(Genes_scores_GEP_5))

genes_scores_list <- list(
  Genes_scores_GEP_1 %>% select(genes),
  Genes_scores_GEP_2 %>% select(genes),
  Genes_scores_GEP_3 %>% select(genes),
  Genes_scores_GEP_4 %>% select(genes),
  Genes_scores_GEP_5 %>% select(genes)
)

# Combine the list of data frames into a single data frame
genes_df <- data.frame(genes_scores_list)

# Rename the columns to "GEP_1", "GEP_2", "GEP_3", and "GEP_4"
colnames(genes_df) <- c("GEP_1", "GEP_2", "GEP_3", "GEP_4", "GEP_5")

rownames(genes_df) <- NULL
# Print the resulting data frame
Top_genes <- genes_df[1:20,]

SCpubr::do_FeaturePlot(sample = thymus.nkt, 
                       features = c("CRTAM"), order = T,
                       plot.title = "", legend.position = "none",
                       viridis_color_map = "inferno")

