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

thymus.GD <- preprocess(seurobj = subset(thymus.filt, subset=group.ident=="GD_Thymus"),
                         celltype = "GDthy", colvalues=colvalues, ndim=10,
                         res=0.4)

DimPlot(thymus.GD, group.by="RNA_snn_res.0.4", label=T, repel=T) +
  theme(legend.position="none") + 
  # scale_color_manual(values=colvalues) +
  labs(title="Thymic GD")

thymus.GD$new_clusters_GD <- case_when(
  thymus.GD$RNA_snn_res.0.4 == '7' ~ 'GD_c0',
  thymus.GD$RNA_snn_res.0.4 == '1' ~ 'GD_c1',
  thymus.GD$RNA_snn_res.0.4 == '2' ~ 'GD_c2',
  thymus.GD$RNA_snn_res.0.4 == '4' ~ 'GD_c3',
  thymus.GD$RNA_snn_res.0.4 == '6' ~ 'GD_c4',
  thymus.GD$RNA_snn_res.0.4 == '0' ~ 'GD_c5',
  thymus.GD$RNA_snn_res.0.4 == '3' ~ 'GD_c6',
  thymus.GD$RNA_snn_res.0.4 == '5' ~ 'GD_c7'
)

thymus.GD$new_clusters_GD  <- as.factor(thymus.GD$new_clusters_GD)
Idents(thymus.GD) <- "new_clusters_GD"

colors_clusters <- c("0" = "#f4c40f", "1" = "#b75347", "2" = "#d8443c", "3" = "#e09351", "4" = "#2b9b81", 
                     "5" = "#421401", "6" = "#92c051", "7" = "#9f5691", "8" = "#17154f", "9" = "#74c8c3", 
                     "10" = "#5a97c1", "11" = "gold", "12" = "#a40000", "13" = "#72bcd5", "14" = "grey50",
                     "15" = "orange", "16" = "blueviolet", "17" = "#0a2e57", "18" = "olivedrab2")

colors_clusters_GD <- c("GD_c0" = "#d8443c", "GD_c1" = "#e09351", "GD_c2" = "gold", "GD_c3" = "#9f5691", "GD_c4" = "#72bcd5",
                         "GD_c5" = "blueviolet", "GD_c6" = "olivedrab2", "GD_c7" = "grey50")

DimPlot(thymus.GD, cols = colors_clusters_GD, pt.size = 1, label = T)

plot <- SCpubr::do_DimPlot(sample = thymus.GD, 
                           label = F, repel = T, font.size = 20,
                           label.color = "white", colors.use = colors_clusters_GD,
                           border.color = "black", border.size = 2,
                           reduction = "umap", pt.size = 1.2)


ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/GD_analysis_clusters_SCPubR.pdf",
       plot, device = "pdf", width=8, height=7, limitsize = F)

thymus.GD <- NormalizeData(thymus.GD)

Cd8aa_score <- list(c("GNG4", "CD8A", "NUCB2", "LEF1", "PRKCH", "PTPN7", "SH2D1A", "ELOVL5", "TOX",
                      "AQP3", "BCL6", "ITGA4", "MYB", "RTKN2", "IKZF2", "DUSP2", "MINDY2", "HIVEP3"))

thymus.GD <- AddModuleScore(object = thymus.GD, features = Cd8aa_score,
                             assay = "RNA", name = 'Cd8aa_score')

plot <- SCpubr::do_FeaturePlot(sample = thymus.GD, 
                               features = "Cd8aa_score1", order = T,
                               plot.title = "", legend.position = "none", #min.cutoff = c(-0.50), max.cutoff = c(1),
                               border.color = "black", border.size = 2,
                               reduction = "umap", pt.size = 1.2) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0.1, end = 0.85, direction = -1)

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/CD8aa_signature_GD.pdf",
       plot, device = "pdf", width=8, height=7, limitsize = F)

Egress <- list(c("KLF2", "CORO1A", "CCR7", "CXCR4", "CXCR6", "FOXO1", "CXCR3", "S1PR1", "S1PR4",
                 "S100A4", "S100A6", "EMP3"))

thymus.GD <- AddModuleScore(object = thymus.GD, features = Egress,
                             assay = "RNA", name = 'Egress_score')

plot <- SCpubr::do_FeaturePlot(sample = thymus.GD, 
                               features = "Egress_score1", order = T,
                               plot.title = "", legend.position = "none", #min.cutoff = c(-0.50), max.cutoff = c(1),
                               border.color = "black", border.size = 2,
                               reduction = "umap", pt.size = 1.2) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0.1, end = 0.85, direction = -1)

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/Egress_signature_GD.pdf",
       plot, device = "pdf", width=8, height=7, limitsize = F)


###Distribution of cells per Donor####################################################################################
Idents(thymus.GD) <- "new_clusters_GD"

# Extract and summarize NKT cell numbers by donor and new_clusters_NKT
GD_cell_numbers_thymus <- thymus.GD@meta.data %>%
  as.data.frame()
GD_cell_numbers_thymus %>% 
  select(Donor, new_clusters_GD) %>%
  dplyr::count(new_clusters_GD, Donor, name = "N") %>%
  filter(!is.na(new_clusters_GD) & !is.na(Donor)) %>%
  ungroup()

# Compute frequency of NKT cells for each donor and new_clusters_NKT
GD_cell_numbers_thymus_test <- GD_cell_numbers_thymus %>%
  select(Donor, new_clusters_GD) %>%
  dplyr::count(new_clusters_GD, Donor, name = "N") %>%
  filter(!is.na(new_clusters_GD) & !is.na(Donor)) %>%
  group_by(Donor) %>%
  mutate(freq = N/sum(N)*100) %>%
  ungroup()

# Combine frequency data for all donors and new_clusters_NKT
distribution_donors_thymus <- GD_cell_numbers_thymus_test %>%
  complete(new_clusters_GD, Donor, fill=list(N=0, freq=0)) %>%
  mutate(Donor = factor(Donor))

plot <- ggplot(data= subset(distribution_donors_thymus, !is.na(Donor)), aes(x = Donor, y = freq, 
                                                                            fill = forcats::fct_rev(factor(new_clusters_GD)))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = colors_clusters_GD, 
                    labels = levels(distribution_donors_thymus$new_clusters_GD), 
                    drop = FALSE, limits = levels(distribution_donors_thymus$new_clusters_GD),
                    guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Donors", y = "%", fill = "Clusters") + theme_classic() + My_Theme 

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/Distribution_Thymic_GD_per_donor.pdf",
       plot, device = "pdf", width=7, height=7, limitsize = F)

##DEGs#################################################################################################################
Idents(thymus.GD) <- "new_clusters_GD"
GD.clusters.markers <- FindAllMarkers(thymus.GD, test.use = 'wilcox', 
                                       logfc.threshold = 0.4, min.pct = 0.3, only.pos = TRUE)
top5 <- GD.clusters.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5_distinct <- base::unique(top5$gene)

write.csv(GD.clusters.markers, "/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/clusters.markers.GD.csv")

GD_Thymus_samples <- ScaleData(thymus.GD, features = top5_distinct)
DotPlot_colors <- c("lightsteelblue1", "#1E466E")

plot <- DotPlot(GD_Thymus_samples, features = top5_distinct, dot.scale = 8, cols = DotPlot_colors,
                col.min = -1, col.max = 2, dot.min = 0) + RotatedAxis()
ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/Features_per_GD_clusters.pdf",
       plot, device = "pdf", width=10, height=5, limitsize = F)

##################

SCpubr::do_FeaturePlot(sample = thymus.GD, 
                       features = c("GZMK"), order = T,
                       plot.title = "", legend.position = "none",
                       viridis_color_map = "inferno")

plot <- SCpubr::do_FeaturePlot(sample = thymus.GD,
                               features = c("ZBTB16", "CCR9", "CCR7", "CD4", "FOS", "KLRB1"), order = T,
                              plot.title = "", legend.position = "none", min.cutoff = c(0, 0, 0, 0, 0, 0), max.cutoff = c(2, 2, 2, 2, 3, 2),
                              border.color = "black", border.size = 1.2,
                              reduction = "umap", pt.size = 0.6) &
                                 scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0.1, end = 0.85, direction = -1)

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/Genes1_GD.pdf", 
       plot, device = "pdf", width=10, height=18)

plot <- SCpubr::do_FeaturePlot(sample = thymus.GD, 
                               features = c("EOMES", "GZMK", "CD8A", "CCR6", "MKI67", "CD1C"), order = T,
                               plot.title = "", legend.position = "none", min.cutoff = c(0, 0, 0, 0, 0, 0), max.cutoff = c(2, 2, 2, 2, 3, 2),
                               border.color = "black", border.size = 1.2,
                               reduction = "umap", pt.size = 0.6) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0.1, end = 0.85, direction = -1)

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/Genes2_GD.pdf", 
       plot, device = "pdf", width=10, height=18)

plot <- SCpubr::do_FeaturePlot(sample = thymus.GD, 
                               features = c("CD1C", "TCF7", "GATA3", "SELL", "GNG4", "RAG1"), order = T,
                               plot.title = "", legend.position = "none",
                               viridis_color_map = "inferno")

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/Genes3_GD.pdf", 
       plot, device = "pdf", width=10, height=18)

#############################################################################################
TCR_data_GD_cells.df <- as.data.frame(thymus.GD@meta.data)
TCR_data_GD_cells.df <- TCR_data_GD_cells.df %>%
  dplyr::select(-nCount_RNA, -nFeature_RNA, -percent.mt, -cell.ident, -Sex, -Age_in_weeks, -Donor, -Batch, -Method,
                -RNA_snn_res.0.7, -seurat_clusters, -orig.ident, -Tissue, -Test_cell_ident_for_cLISI,
                -group.ident, -TCRa_g_chain, -TCRa_g_clonotype, -TCRb_d_chain, -TCRb_d_clonotype, 
                -TRAV10_TRAJ18, -TRAV1_TRAJ33, -TRAV1, -new_clusters, -Effectorness1, -Naiveness1, -Type_3_score1, -new_clusters_id,
                -Cd8aa_score1) %>%
  na.omit() %>%
  dplyr::rename(TRGV = TCR_Alpha_Gamma_V_gene_Dominant, TRGJ = TCR_Alpha_Gamma_J_gene_Dominant, CDR3g = TCR_Alpha_Gamma_CDR3_Translation_Dominant,
                TRDV = TCR_Beta_Delta_V_gene_Dominant, TRDD = TCR_Beta_Delta_D_gene_Dominant, TRDJ = TCR_Beta_Delta_J_gene_Dominant,
                CDR3d = TCR_Beta_Delta_CDR3_Translation_Dominant) %>%
  mutate(CDR3d = paste0("C", CDR3d)) %>%
  mutate(CDR3g = paste0("C", CDR3g))

Total_TRDV_TRDJ_distribution <- TCR_data_GD_cells.df %>%
  dplyr::select(-CDR3g, -TRGV, -TRDD, -TRGJ, -CDR3d) %>%
  mutate(TRDV = str_remove(TRDV, pattern = "\\*[^.]*$"),
         TRDJ = str_remove(TRDJ, pattern = "\\*[^.]*$")) %>%
  group_by(TRDV, TRDJ) %>%
  summarise(n = n()) %>% ungroup() %>% as.data.frame()

table(TCR_data_GD_cells.df$TRDV)
TRDV1_cells <- TCR_data_GD_cells.df %>% dplyr::filter(TRDV %in% c("TRDV1*01")) %>% row.names()

TRDV2_cells <- TCR_data_GD_cells.df %>% dplyr::filter(TRDV %in% c("TRDV2*01", "TRDV2*03")) %>% row.names()

TRDV3_cells <- TCR_data_GD_cells.df %>% dplyr::filter(TRDV %in% c("TRDV3*01")) %>% row.names()

plot <- SCpubr::do_DimPlot(sample = thymus.GD, 
                        label = FALSE,
                        cells.highlight = TRDV1_cells,
                        legend.position = "none",
                        na.value = "grey90",
                        sizes.highlight = 4, colors.use = "#a40000")

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/TRDV1_GD.pdf", 
       plot, device = "pdf", width=7, height=7)

plot <- SCpubr::do_DimPlot(sample = thymus.GD, 
                        label = FALSE,
                        cells.highlight = TRDV2_cells,
                        legend.position = "none",
                        na.value = "grey90",
                        sizes.highlight = 4, colors.use = "blue")

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/TRDV2_GD.pdf", 
       plot, device = "pdf", width=7, height=7)


plot <- SCpubr::do_DimPlot(sample = thymus.GD, 
                        label = FALSE,
                        cells.highlight = TRDV3_cells,
                        legend.position = "none",
                        na.value = "grey90",
                        sizes.highlight = 4, colors.use = "#318f49")

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/TRDV3_GD.pdf", 
       plot, device = "pdf", width=7, height=7)

table(TCR_data_GD_cells.df$TRGV)
TRGV9_cells <- TCR_data_GD_cells.df %>% dplyr::filter(TRGV %in% c("TRGV9*01")) %>% row.names()

plot <- SCpubr::do_DimPlot(sample = thymus.GD, 
                        label = FALSE,
                        cells.highlight = TRGV9_cells, 
                        legend.position = "none",
                        na.value = "grey90",
                        sizes.highlight = 4, colors.use = "gold")

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/TRGV9_GD.pdf", 
       plot, device = "pdf", width=7, height=7)


######################
TRDV1_TCRs <- TCR_data_GD_cells.df %>% dplyr::select(-RNA_snn_res.1.2, -Egress_score1) %>%
  filter(str_detect(TRDV, "\\TRDV1")) %>%
  dplyr::mutate(TRGV = str_remove(TRGV, pattern = "\\*[^.]*$"), 
                TRGJ = str_remove(TRGJ, pattern = "\\*[^.]*$"),
                TRDV = str_remove(TRDV, pattern = "\\*[^.]*$"),
                TRDJ = str_remove(TRDJ, pattern = "\\*[^.]*$"))

TRDV1_TRDJ_naive <- TRDV1_TCRs  %>%
  group_by(TRDV, TRDJ) %>%
  summarise(n = n()) %>% ungroup() %>% as.data.frame()

cols <- hcl.colors(20, "Temps")

Vdelta_colors_vector <- c("TRDV1" = "#a40000", "TRDJ1" = "#16317d", "TRDJ2" = "#3BB38C", 
                          "TRDJ3" = "#007e2f", "TRDJ4" = "#ffcd12")
chordWrap <- function(TRDV1_TRDJ_naive){
  circos.clear()
  chordDiagram(TRDV1_TRDJ_naive, annotationTrack = c("grid"), 
               annotationTrackHeight = c(0.1, 0.1),
               grid.col = Vdelta_colors_vector,
               row.col = c("#A40000"),
               order = c(rep("TRDV1", 4), "TRDJ1", "TRDJ2", "TRDJ3", "TRDJ4"),
               transparency = 0.3,
               link.lwd = 0.3, 
               link.lty = 5, 
               link.border = 1,
               col = Vdelta_colors_vector, 
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(TRDV1_TRDJ_naive))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA) 
}

chordWrap(TRDV1_TRDJ_naive)

########################
TRGV_TRDV1_TRDJ_naive <- TRDV1_TCRs  %>%
  group_by(TRGV, TRGJ) %>%
  summarise(n = n()) %>% ungroup() %>% as.data.frame()

cols <- hcl.colors(25, "Temps")

Vgamma_colors_vector <- c("TRGV1" = "#CF597E", "TRGV2" = "#a40000", "TRGV3" = "#16317d", "TRGV4" = "#007e2f", "TRGV5" = "#ffcd12",
                          "TRGV5P" = "#D86279", "TRGV8" = "#1FA990", "TRGV9" = "#E9A26A", "TRGV10" = "#DB6577", "TRGV11" = "#E16C72", 
                          "TRGJ1" = "#93CB83", 
                          "TRGJP" = "#C7D88D", "TRGJP1" = "#EADB94", "TRGJP2" = "#EACE85")

chordWrap <- function(TRGV_TRDV1_TRDJ_naive){
  circos.clear()
  chordDiagram(TRGV_TRDV1_TRDJ_naive, annotationTrack = c("grid"), 
               annotationTrackHeight = c(0.1, 0.1),
               grid.col = Vgamma_colors_vector,
               row.col = c("#A40000"),
               order = c(rep("TRGV1", 2), rep("TRGV2", 3), rep("TRGV3", 4), rep("TRGV4", 3), rep("TRGV5", 3), 
                         rep("TRGV5P", 1), rep("TRGV8", 4), rep("TRGV9", 4), rep("TRGV10", 3), rep("TRGV11", 2),
                         rep("TRGJ1", 8), rep("TRGJP", 3), rep("TRGJP1", 9), rep("TRGJP2", 9)),
               transparency = 0.3,
               link.lwd = 0.3, 
               link.lty = 5, 
               link.border = 1,
               col = Vdelta_colors_vector, 
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(TRGV_TRDV1_TRDJ_naive))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA) 
}

chordWrap(TRGV_TRDV1_TRDJ_naive)


###################################
TRDV2_TCRs <- TCR_data_GD_cells.df %>% dplyr::select(-RNA_snn_res.1.2, -Egress_score1) %>%
  filter(str_detect(TRDV, "\\TRDV2")) %>%
  dplyr::mutate(TRGV = str_remove(TRGV, pattern = "\\*[^.]*$"), 
                TRGJ = str_remove(TRGJ, pattern = "\\*[^.]*$"),
                TRDV = str_remove(TRDV, pattern = "\\*[^.]*$"),
                TRDJ = str_remove(TRDJ, pattern = "\\*[^.]*$"))

TRDV2_TCRs_naive <- TRDV2_TCRs %>% filter(new_clusters_GD %in% c("GD_c1", "GD_c2", "GD_c3", "GD_c5", "GD_c6"))

TRDV2_TRDJ_naive <- TRDV2_TCRs_naive  %>%
  group_by(TRDV, TRDJ) %>%
  summarise(n = n()) %>% ungroup() %>% as.data.frame()

cols <- hcl.colors(20, "Temps")

Vdelta_colors_vector <- c("TRDV2" = "#a40000", "TRDJ1" = "#16317d", "TRDJ3" = "#007e2f", "TRDJ4" = "#ffcd12")
chordWrap <- function(TRDV2_TRDJ_naive){
  circos.clear()
  chordDiagram(TRDV2_TRDJ_naive, annotationTrack = c("grid"), 
               annotationTrackHeight = c(0.1, 0.1),
               grid.col = Vdelta_colors_vector,
               row.col = c("#A40000"),
               order = c(rep("TRDV2", 3), "TRDJ1", "TRDJ3", "TRDJ4"),
               transparency = 0.3,
               link.lwd = 0.3, 
               link.lty = 5, 
               link.border = 1,
               col = Vdelta_colors_vector, 
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(TRDV2_TRDJ_naive))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA) 
}

chordWrap(TRDV2_TRDJ_naive)

############################
TRDV2_TCRs_naive <- TRDV2_TCRs %>% filter(new_clusters_GD %in% c("GD_c1", "GD_c2", "GD_c3", "GD_c5", "GD_c6"))

TRGV_TRDV2_TRDJ_naive <- TRDV2_TCRs_naive  %>%
  group_by(TRGV, TRGJ) %>%
  summarise(n = n()) %>% ungroup() %>% as.data.frame()

cols <- hcl.colors(20, "Temps")

Vgamma_colors_vector <- c("TRGV2" = "#a40000", "TRGV3" = "#16317d", "TRGV4" = "#007e2f", "TRGV5" = "#ffcd12",
                          "TRGV8" = "#1FA990", "TRGV9" = "#E9A26A", "TRGV10" = "#DB6577", "TRGJ1" = "#93CB83", 
                          "TRGJP" = "#C7D88D", "TRGJP1" = "#EADB94", "TRGJP2" = "#EACE85")
chordWrap <- function(TRGV_TRDV2_TRDJ_naive){
  circos.clear()
  chordDiagram(TRGV_TRDV2_TRDJ_naive, annotationTrack = c("grid"), 
               annotationTrackHeight = c(0.1, 0.1),
               grid.col = Vgamma_colors_vector,
               row.col = c("#A40000"),
               order = c(rep("TRGV2", 1), "TRGV3", "TRGV4", "TRGV5", "TRGV8", "TRGV9", "TRGV10",
                         "TRGJ1", "TRGJP", "TRGJP1", "TRGJP2"),
               transparency = 0.3,
               link.lwd = 0.3, 
               link.lty = 5, 
               link.border = 1,
               col = Vdelta_colors_vector, 
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(TRGV_TRDV2_TRDJ_naive))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA) 
}

chordWrap(TRGV_TRDV2_TRDJ_naive)

#################################
TRDV2_TCRs_effector <- TRDV2_TCRs %>% filter(new_clusters_GD %in% c("GD_c7"))

TRDV2_TCRs_effector <- TRDV2_TCRs_effector %>%
  group_by(TRDV, TRDJ) %>%
  summarise(n = n()) %>% ungroup() %>% as.data.frame()

cols <- hcl.colors(20, "Temps")

Vdelta_colors_vector <- c("TRDV2" = "#a40000", "TRDJ1" = "#16317d", "TRDJ2" = "#3BB38C", "TRDJ3" = "#007e2f")
chordWrap <- function(TRDV2_TCRs_effector){
  circos.clear()
  chordDiagram(TRDV2_TCRs_effector, annotationTrack = c("grid"), 
               annotationTrackHeight = c(0.1, 0.1),
               grid.col = Vdelta_colors_vector,
               row.col = c("#A40000"),
               order = c(rep("TRDV2", 3), "TRDJ1", "TRDJ2", "TRDJ3"),
               transparency = 0.3,
               link.lwd = 0.1, 
               link.lty = 5, 
               link.border = 0.1,
               col = Vdelta_colors_vector, 
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(TRDV2_TCRs_effector))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA) 
}

chordWrap(TRDV2_TCRs_effector)

TRDV2_TCRs_effector <- TRDV2_TCRs %>% filter(new_clusters_GD %in% c("GD_c7"))

TRGV_TRDV2_TRDJ_effector <- TRDV2_TCRs_effector  %>%
  group_by(TRGV, TRGJ) %>%
  summarise(n = n()) %>% ungroup() %>% as.data.frame()

cols <- hcl.colors(20, "Temps")

Vgamma_colors_vector <- c("TRGV2" = "#a40000", "TRGV3" = "#16317d", "TRGV4" = "#007e2f", "TRGV5" = "#ffcd12",
                          "TRGV8" = "#1FA990", "TRGV9" = "#E9A26A", "TRGV10" = "#DB6577", "TRGJ1" = "#93CB83", 
                          "TRGJP" = "#C7D88D", "TRGJP1" = "#EADB94", "TRGJP2" = "#EACE85")

chordWrap <- function(TRGV_TRDV2_TRDJ_effector){
  circos.clear()
  chordDiagram(TRGV_TRDV2_TRDJ_effector, annotationTrack = c("grid"), 
               annotationTrackHeight = c(0.1, 0.1),
               grid.col = Vgamma_colors_vector,
               row.col = c("#A40000"),
               order = c(rep("TRGV2", 1), "TRGV3", "TRGV5", "TRGV8", rep("TRGV9", 3), rep("TRGV10", 2),
                         rep("TRGJ1", 5), "TRGJP", rep("TRGJP1", 3)),
               transparency = 0.3,
               link.lwd = 0.1, 
               link.lty = 5, 
               link.border = 0.1,
               col = Vdelta_colors_vector, 
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(TRGV_TRDV2_TRDJ_effector))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA) 
}

chordWrap(TRGV_TRDV2_TRDJ_effector)

##########################
TRDV2_TRGV9.df <- TCR_data_GD_cells.df %>% dplyr::select(-RNA_snn_res.1.2, -Egress_score1, -RNA_snn_res.0.4) %>%
  mutate(TRDV = str_remove(TRDV, pattern = "\\*[^.]*$"),
         TRDJ = str_remove(TRDJ, pattern = "\\*[^.]*$"),
         TRGV = str_remove(TRGV, pattern = "\\*[^.]*$"),
         TRGJ = str_remove(TRGJ, pattern = "\\*[^.]*$")) %>%
  filter(TRDV %in% c("TRDV2") & TRGV == "TRGV9" &
           TRGJ == "TRGJP") %>% 
  mutate(CDR3d_length = str_length(CDR3d),
         CDR3g_length = str_length(CDR3g))

pdf(paste0(fig_dir, "CDR3 size TRGV9_TRGJP_of_TRGV9_PBMC.pdf"), width = 7, height = 7)
p <- TRDV2_TRGV9.df %>% group_by(CDR3_length) %>% summarize(N = n()) %>%
  ggplot(aes(x = CDR3_length, y = N)) + geom_bar(stat="identity", width = 0.5, fill = "#08306b") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 350)) +
  scale_x_continuous(limits = c(10, 32), breaks = c(seq(from = 10, to = 32, by = 1))) +
  labs(x = "CDR3 Length (aa)", y = "Nb of Sequences") +
  theme_classic(base_size = 14) + My_Theme
print(p)
dev.off()

pdf(paste0(fig_dir, "CDR3 logo TRGV9_TRGJP-TRDV2_14aa.pdf"), width = 11, height = 3)
p <- TCR_data_GD_cells.length.df %>%
  dplyr::filter(CDR3_length == "14") %>% 
  select(CDR3g) %>% ggseqlogo(seq_type = 'aa', method = 'bits') +
  theme(axis.text=element_text(size=28), axis.title=element_text(size=28, face="bold")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 28, b = 0, l = 0))) +
  theme(legend.position="none")
print(p)
dev.off()

pdf(paste0(fig_dir, "CDR3 logo TRGV9_TRGJP-TRDV2_15aa.pdf"), width = 11, height = 3)
p <- TCR_data_GD_cells.length.df %>%
  dplyr::filter(CDR3_length == "15") %>% 
  select(CDR3g) %>% ggseqlogo(seq_type = 'aa', method = 'bits') +
  theme(axis.text=element_text(size=28), axis.title=element_text(size=28, face="bold")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 28, b = 0, l = 0))) +
  theme(legend.position="none")
print(p)
dev.off()

pdf(paste0(fig_dir, "CDR3 logo TRGV9_TRGJP-TRDV2_16aa.pdf"), width = 11, height = 3)
p <- TCR_data_GD_cells.length.df %>%
  dplyr::filter(CDR3_length == "16") %>% 
  select(CDR3g) %>% ggseqlogo(seq_type = 'aa', method = 'bits') +
  theme(axis.text=element_text(size=28), axis.title=element_text(size=28, face="bold")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 28, b = 0, l = 0))) +
  theme(legend.position="none")
print(p)
dev.off()



GD2_CDR3_size <- GD2_TCRS %>% filter(TRGV == "TRGV9", !str_detect(CDR3d, pattern = "\\*")) %>% 
  mutate(CDR3_d_length = str_length(CDR3d)) %>% as.data.frame()

ggplot(GD2_CDR3_size) + 
  geom_bar(aes(x=CDR3_d_length, y = ..prop..,)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  labs(x = "CDR3 Delta Length", y = "Frequency (%)") +
  theme_classic() + NoLegend()

GD2_pairing <- GD2_TCRS  %>%
  group_by(TRGV, TRGJ) %>%
  summarise(n = n()) %>% ungroup() %>% as.data.frame()

chordWrap <- function(GD2_pairing){
  circos.clear()
  chordDiagram(GD2_pairing, annotationTrack = c("grid"), annotationTrackHeight = c(0.05, 0.05),
               col = cols, transparency = 0.1, link.lwd = 1, link.lty = 1, link.border = 1,
               order = c("TRGV1", rep("TRGV2", 3), rep("TRGV3", 3), rep("TRGV4", 2), rep("TRGV5", 2),
                         rep("TRGV8", 3), rep("TRGV9", 4), rep("TRGV10", 3), rep("TRGV11", 2),
                         rep("TRGJ1", 9), rep("TRGJP1", 6), rep("TRGJP2", 6), rep("TRGJP", 2)),
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(GD2_pairing))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA) 
}

chordWrap(GD2_pairing)

GD2_TRGV9 <- GD2_TCRS %>% filter(TRGV == "TRGV9", TRGJ == "TRGJP") %>%
  group_by(CDR3g) %>%
  summarise(n = n()) %>% ungroup() %>% as.data.frame()
###################################################################

saveRDS(thymus.GD, "/Volumes/Samsung_T5/Human_MAIT_NKT/Data/thymus.GD.03_09_23.RDS")

thymus.GD <- readRDS("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Data/thymus.GD.03_09_23.RDS")

###################################################################
