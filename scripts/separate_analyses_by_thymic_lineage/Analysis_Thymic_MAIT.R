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
preprocess <- function(seurobj, ndim = 20, res = 0.8, colvalues, celltype, 
                       approxPCA = FALSE){
  print(seurobj)
  cat("Cells:", unique(seurobj$group.ident), "\n")
  
  print("Normalize")
  seurobj <- NormalizeData(seurobj)
  print("HVG")
  seurobj <- FindVariableFeatures(seurobj, selection.method="vst",
                                  nfeatures = 2000)
  print("Scale")
  seurobj <- ScaleData(seurobj)
  
  print("PCA")
  seurobj <- RunPCA(seurobj, npcs = 50, verbose = FALSE, approx=approxPCA)
  
  print("Harmony")
  seurobj <- RunHarmony(seurobj, reduction = "pca", max.iter.harmony = 30,
                        group.by.vars = c("Method", "Batch"))
  print(ElbowPlot(seurobj, ndims = 50, reduction="harmony"))
  
  print("UMAP")
  seurobj <- RunUMAP(seurobj, reduction = "harmony", dims = 1:ndim,
                     n.neighbors = 50, reduction.key = "UMAP50_")
  
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

thymus.MAIT <- preprocess(seurobj = subset(thymus.filt, subset=group.ident=="MAIT_Thymus"),
                         celltype = "MAITthy", colvalues = colvalues, ndim = 10,
                         res = 0.2)

DimPlot(thymus.MAIT, group.by="RNA_snn_res.0.2", label=T, repel=T) +
  theme(legend.position="none") + 
  labs(title="Thymic MAIT")

thymus.MAIT$new_clusters_MAIT <- case_when(
  thymus.MAIT$RNA_snn_res.0.2 == '5' ~ 'MAIT_c0',
  thymus.MAIT$RNA_snn_res.0.2 == '4' ~ 'MAIT_c1',
  thymus.MAIT$RNA_snn_res.0.2 == '0' ~ 'MAIT_c2',
  thymus.MAIT$RNA_snn_res.0.2 == '3' ~ 'MAIT_c3',
  thymus.MAIT$RNA_snn_res.0.2 == '1' ~ 'MAIT_c4',
  thymus.MAIT$RNA_snn_res.0.2 == '6' ~ 'MAIT_c5',
  thymus.MAIT$RNA_snn_res.0.2 == '2' ~ 'MAIT_c6'
)

thymus.MAIT$new_clusters_MAIT  <- as.factor(thymus.MAIT$new_clusters_MAIT)
Idents(thymus.MAIT) <- "new_clusters_MAIT"

colors_clusters <- c("0" = "#f4c40f", "1" = "#b75347", "2" = "#d8443c", "3" = "#e09351", "4" = "#2b9b81", 
                     "5" = "#421401", "6" = "#92c051", "7" = "#9f5691", "8" = "#17154f", "9" = "#74c8c3", 
                     "10" = "#5a97c1", "11" = "gold", "12" = "#a40000", "13" = "#72bcd5", "14" = "grey50",
                     "15" = "orange", "16" = "blueviolet", "17" = "#0a2e57", "18" = "olivedrab2")

colors_clusters_MAIT <- c("MAIT_c0" = "#d8443c", "MAIT_c1" = "#e09351", "MAIT_c2" = "gold", "MAIT_c3" = "#74c8c3", 
                          "MAIT_c4" = "#a40000", "MAIT_c5" = "#5a97c1", "MAIT_c6" = "orange")

DimPlot(thymus.MAIT, cols = colors_clusters_MAIT, pt.size = 1, label = T)

plot <- SCpubr::do_DimPlot(sample = thymus.MAIT, 
                           label = F, repel = T, font.size = 20,
                           label.color = "white", colors.use = colors_clusters_MAIT,
                           border.color = "black", border.size = 2,
                           reduction = "umap", pt.size = 1.2)

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/MAIT_analysis_clusters_SCPubR.pdf",
       plot, device = "pdf", width=7, height=7, limitsize = F)

thymus.MAIT <- NormalizeData(thymus.MAIT)

Cd8aa_score <- list(c("GNG4", "CD8A", "NUCB2", "LEF1", "PRKCH", "PTPN7", "SH2D1A", "ELOVL5", "TOX",
                      "AQP3", "BCL6", "ITGA4", "MYB", "RTKN2", "IKZF2", "DUSP2", "MINDY2", "HIVEP3"))

thymus.MAIT <- AddModuleScore(object = thymus.MAIT, features = Cd8aa_score,
                             assay = "RNA", name = 'Cd8aa_score')

plot <- SCpubr::do_FeaturePlot(sample = thymus.MAIT, 
                               features = "Cd8aa_score1", order = T,
                               plot.title = "", legend.position = "none", min.cutoff = c(-0.50), max.cutoff = c(1),
                               border.color = "black", border.size = 2,
                               reduction = "umap", pt.size = 1.2) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0, end = 1, direction = -1)

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/CD8aa_signature_MAIT.pdf",
       plot, device = "pdf", width=8, height=7, limitsize = F)

Egress <- list(c("KLF2", "CORO1A", "CCR7", "CXCR4", "CXCR6", "FOXO1", "CXCR3", "S1PR1", "S1PR4",
                 "S100A4", "S100A6", "EMP3"))

thymus.MAIT <- AddModuleScore(object = thymus.MAIT, features = Egress,
                             assay = "RNA", name = 'Egress_score')

plot <- SCpubr::do_FeaturePlot(sample = thymus.MAIT, 
                               features = "Egress_score1", order = T,
                               plot.title = "", legend.position = "none", min.cutoff = c(-0.50), max.cutoff = c(0.5),
                               border.color = "black", border.size = 2,
                               reduction = "umap", pt.size = 1.2) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0, end = 1, direction = -1)

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/Egress_signature_MAIT.pdf",
       plot, device = "pdf", width=8, height=7, limitsize = F)

###Distribution of cells per Donor####################################################################################
Idents(thymus.MAIT) <- "new_clusters_MAIT"

# Extract and summarize NKT cell numbers by donor and new_clusters_MAIT
MAIT_cell_numbers_thymus <- thymus.MAIT@meta.data %>%
  as.data.frame()
MAIT_cell_numbers_thymus <-  MAIT_cell_numbers_thymus %>% select(Donor, new_clusters_MAIT) %>%
  dplyr::count(new_clusters_MAIT, Donor) %>%
  rename(N = n) %>%
  ungroup()

# Compute frequency of NKT cells for each donor and new_clusters_MAIT
MAIT_cell_numbers_thymus_test <- MAIT_cell_numbers_thymus %>%
  group_by(Donor) %>%
  mutate(freq = N/sum(N)*100) %>%
  ungroup()

# Combine frequency data for all donors and new_clusters_MAIT
distribution_donors_thymus <- MAIT_cell_numbers_thymus_test %>%
  complete(new_clusters_MAIT, Donor, fill=list(N=0, freq=0)) %>%
  mutate(Donor = factor(Donor))

plot <- ggplot(data= subset(distribution_donors_thymus, !is.na(Donor)), aes(x = Donor, y = freq, 
                                                                            fill = forcats::fct_rev(factor(new_clusters_MAIT)))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = colors_clusters_MAIT, 
                    labels = levels(distribution_donors_thymus$new_clusters_MAIT), 
                    drop = FALSE, limits = levels(distribution_donors_thymus$new_clusters_MAIT),
                    guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Donors", y = "%", fill = "Clusters") + theme_classic() + My_Theme 

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/Distribution_Thymic_MAIT_per_donor.pdf",
       plot, device = "pdf", width=7, height=7, limitsize = F)

##DEGs#################################################################################################################
Idents(thymus.MAIT) <- "new_clusters_MAIT"
MAIT.clusters.markers <- Seurat::FindAllMarkers(thymus.MAIT, test.use = 'wilcox', 
                                       logfc.threshold = 0.4, min.pct = 0.3, only.pos = TRUE)
top5 <- MAIT.clusters.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5_distinct <- base::unique(top5$gene)

write.csv(MAIT.clusters.markers, "/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/clusters.markers.MAIT.csv")

MAIT_Thymus_samples <- ScaleData(thymus.MAIT, features = top5_distinct)
DotPlot_colors <- c("lightsteelblue1", "#1E466E")

plot <- DotPlot(MAIT_Thymus_samples, features = top5_distinct, dot.scale = 8, cols = DotPlot_colors,
                col.min = -1, col.max = 2, dot.min = 0) + RotatedAxis()
ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/Features_per_MAIT_clusters.pdf",
       plot, device = "pdf", width=10, height=5, limitsize = F)

##################PLOTTIÑG Various Genes

plot <- SCpubr::do_FeaturePlot(sample = thymus.MAIT, 
                               features = c("ZBTB16", "CCR9", "CCR7", "CD4", "FOS", "KLRB1"), order = T,
plot.title = "", legend.position = "none", min.cutoff = c(0, 0, 0, 0, 0, 0), max.cutoff = c(2, 2, 2, 2, 3, 2),
border.color = "black", border.size = 1.2,
reduction = "umap", pt.size = 0.6) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0.1, end = 0.85, direction = -1)

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/Genes1_MAIT.pdf", 
       plot, device = "pdf", width=10, height=18)

plot <- SCpubr::do_FeaturePlot(sample = thymus.MAIT, 
                               features = c("EOMES", "GZMK", "CD8A", "RORC", "FOS", "KLRB1"), order = T,
                               plot.title = "", legend.position = "none", min.cutoff = c(0, 0, 0, 0, 0, 0), max.cutoff = c(2, 2, 2, 2, 3, 2),
                               border.color = "black", border.size = 1.2,
                               reduction = "umap", pt.size = 0.6) &
  scale_colour_scico(palette = "lapaz", alpha = 0.8, begin = 0.1, end = 0.85, direction = -1)

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/Genes2_MAIT.pdf", 
       plot, device = "pdf", width=10, height=18)

plot <- SCpubr::do_FeaturePlot(sample = thymus.MAIT, 
                               features = c("CCR6", "CD8A", "SLC4A10", "RORC", "IL23R", "NR4A1"), order = T,
                               plot.title = "", legend.position = "none",
                               viridis_color_map = "inferno")

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/Genes3_MAIT.pdf", 
       plot, device = "pdf", width=10, height=18)

### Monocle3 trajectory####################################################################################################
cds <- SeuratWrappers::as.cell_data_set(thymus.MAIT)

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
Idents(thymus.MAIT) <- "new_clusters_MAIT"
list.clusters <- thymus.MAIT@active.ident
cds@clusters$UMAP$clusters <- list.clusters

##assign UMAP coordinates
cds@int_colData@listData$reducedDims$UMAP <- thymus.MAIT@reductions$umap@cell.embeddings

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
  viridis::scale_color_viridis(option = "D") + xlim(-8, 20) + ylim(-7.5, 7.5)

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/Trajectory_MAIT.pdf", 
       plot, device = "pdf", width=7, height=7)

# cells ordered by monocle3 pseudotime
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(new_clusters_MAIT, monocle3_pseudotime, median), 
                        fill = new_clusters_MAIT)) +
  geom_boxplot() + 
  scale_fill_manual(values = colors_clusters_MAIT) +
  ylab("Clusters") + xlab("Pseudotime") + theme_clean()

######TCR analysis########################################################################################################
TCR_data_MAIT_cells.df <- as.data.frame(thymus.MAIT@meta.data)
TCR_data_MAIT_cells.df <- TCR_data_MAIT_cells.df %>%
  dplyr::select(-nCount_RNA, -nFeature_RNA, -percent.mt, -cell.ident, -Sex, -Age_in_weeks, -Donor, -Batch, -Method,
                -seurat_clusters, -orig.ident, -Tissue, -Test_cell_ident_for_cLISI,-RNA_snn_res.0.7, -RNA_snn_res.1.2, -Egress_score1,
                -Effectorness1, -Naiveness1, -Type_3_score1, -new_clusters_id, -RNA_snn_res.0.2, -Cd8aa_score1,
                -group.ident, -TCRa_g_chain, -TCRa_g_clonotype, -TCRb_d_chain, -TCRb_d_clonotype, -TRAV10_TRAJ18, -TRAV1_TRAJ33, -TRAV1, -new_clusters) %>%
  na.omit() %>%
  dplyr::rename(TRAV = TCR_Alpha_Gamma_V_gene_Dominant, TRAJ = TCR_Alpha_Gamma_J_gene_Dominant, CDR3a = TCR_Alpha_Gamma_CDR3_Translation_Dominant,
                TRBV = TCR_Beta_Delta_V_gene_Dominant, TRBD = TCR_Beta_Delta_D_gene_Dominant, TRBJ = TCR_Beta_Delta_J_gene_Dominant,
                CDR3b = TCR_Beta_Delta_CDR3_Translation_Dominant) %>%
  mutate(CDR3b = paste0("C", CDR3b)) %>%
  mutate(CDR3a = paste0("C", CDR3a))

Total_TRAV_TRAJ_distribution <- TCR_data_MAIT_cells.df %>%
  dplyr::select(-CDR3a, -TRBV, -TRBD, -TRBJ, -CDR3b) %>%
  mutate(TRAV = str_remove(TRAV, pattern = "\\*[^.]*$"),
         TRAJ = str_remove(TRAJ, pattern = "\\*[^.]*$")) %>%
  dplyr::filter(TRAV == "TRAV1-2") %>%
  group_by(TRAV, TRAJ) %>%
  summarise(n = n()) %>% ungroup() %>% as.data.frame()

table(TCR_data_MAIT_cells.df$TRAV)
Trav1_cells <- TCR_data_MAIT_cells.df %>% dplyr::filter(TRAV %in% c("TRAV1-1*01", "TRAV1-1*02",
                                                                    "TRAV1-2*01")) %>% row.names()
length(Trav1_cells)

plot <- SCpubr::do_DimPlot(sample = MAIT_Thymus_samples, 
                           label = FALSE,
                           cells.highlight = Trav1_cells,
                           plot.title = "", legend.position = "right",
                           na.value = "grey90", 
                           sizes.highlight = 2, colors.use = "#a40000")

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/TRAV1_cells.pdf", 
       plot, device = "pdf", width=8, height=7)

#determine how many TRAV and TRAJ are in the samples
V_alpha_genes <- unique(Total_TRAV_TRAJ_distribution$TRAV) %>% sort()
J_alpha_genes <- unique(Total_TRAV_TRAJ_distribution$TRAJ) %>% sort()
Valpha_color_vector <- c(V_alpha_genes, J_alpha_genes)  
Valpha_color_vector <- unique(Valpha_color_vector)

#create a named vector where each TRAV and TRAJ has a dedicated color
#["Austria",["#a40000","#16317d","#007e2f","#ffcd12","#721b3e","#b86092","#00b7a7"]]
Valpha_colors_vector <- c("TRAV1-2" = "#9eb4e0", "TRAJ6" = "#2f357c", "TRAJ8" = "#b0799a", "TRAJ11" = "#e69b00", 
                          "TRAJ12" = "#355828",
                          "TRAJ20" = "#16317d", "TRAJ27" = "#007e2f", "TRAJ30" = "#6c5d9e", "TRAJ31" = "#bf3729",
                          "TRAJ32" = "#e48171", "TRAJ33" = "#f5bb50", "TRAJ35" = "#9d9cd5", "TRAJ37" = "#ffcd12", 
                          "TRAJ40" = "#17154f",
                          "TRAJ41" = "#f6b3b0")

par(cex = 1.5, mar = c(0, 0, 0, 0))
chordWrap <- function(Total_TRAV_TRAJ_distribution){
  circos.clear()
  chordDiagram(Total_TRAV_TRAJ_distribution, 
               annotationTrack = c("grid"), 
               annotationTrackHeight = c(0.1, 0.1),
               grid.col = Valpha_colors_vector,
               row.col = c("#A40000"),
               order = c(rep("TRAV1-2", 14), rep("TRAJ6", 1), rep("TRAJ8", 1), rep("TRAJ11", 1),
                         rep("TRAJ12", 1), rep("TRAJ20", 1), rep("TRAJ27", 1),
                         rep("TRAJ30", 1), rep("TRAJ31", 1), rep("TRAJ32", 1),
                         rep("TRAJ33", 1), rep("TRAJ35", 1), rep("TRAJ37", 1),
                         rep("TRAJ40", 1), rep("TRAJ41", 1)),
               link.lwd = 0.1, 
               link.lty = 0.1, 
               link.border = 1,
               col = Valpha_colors_vector, transparency = 0.1,
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(Total_TRAV_TRAJ_distribution))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA) 
}

chordWrap(Total_TRAV_TRAJ_distribution)

Total_TRAV_TRAJ <- TCR_data_MAIT_cells.df %>%
  dplyr::select(- new_clusters_MAIT, -TRBV, -TRBD, -TRBJ, -CDR3b) %>%
  dplyr::filter(TRAV %in% c("TRAV1-1*01", "TRAV1-1*02",
                            "TRAV1-2*01")) %>%
  mutate(TRAV = str_remove(TRAV, pattern = "\\*[^.]*$"),
         TRAJ = str_remove(TRAJ, pattern = "\\*[^.]*$")) %>%
  mutate(CDR3_Length = str_length(CDR3a)) %>%
  dplyr::filter(TRAV == "TRAV1-2", TRAJ == "TRAJ33", CDR3a %in% str_remove(CDR3a, pattern = "\\*[^*]"))

pdf(paste0(fig_dir, "CDR3 size TRAV1-2-TRAJ33.pdf"), width=7, height=7)
p <- Total_TRAV_TRAJ %>% group_by(CDR3_Length) %>% summarise(N = n()) %>%
  ggplot(aes(x = CDR3_Length, y = N)) + geom_bar(stat="identity", width = 0.5, fill = "#08306b") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 300)) +
  scale_x_continuous(limits = c(8, 14), breaks = c(seq(from = 8, to = 14, by = 1))) +
  labs(x = "CDR3 Length (aa)", y = "Nb of Sequences") +
  theme_classic(base_size = 14) + My_Theme
print(p)
dev.off()

pdf(paste0(fig_dir, "CDR3 logo TRAV1-2-TRAJ33.pdf"), width = 10, height = 3)
p <- Total_TRAV_TRAJ %>%
  dplyr::filter(CDR3_Length == "11") %>% 
  dplyr::select(CDR3a) %>% ggseqlogo(seq_type = 'aa', method = 'bits') +
  theme(axis.text=element_text(size=28), axis.title=element_text(size=28, face="bold")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 28, b = 0, l = 0))) +
  theme(legend.position="none")
print(p)
dev.off()

Total_TRAV_TRAJ <- TCR_data_MAIT_cells.df %>%
  dplyr::select(- new_clusters_MAIT, -TRBV, -TRBD, -TRBJ, -CDR3b) %>%
  dplyr::filter(TRAV %in% c("TRAV1-1*01", "TRAV1-1*02",
                            "TRAV1-2*01")) %>%
  mutate(TRAV = str_remove(TRAV, pattern = "\\*[^.]*$"),
         TRAJ = str_remove(TRAJ, pattern = "\\*[^.]*$")) %>%
  mutate(CDR3_Length = str_length(CDR3a)) %>%
  dplyr::filter(TRAV == "TRAV1-2", TRAJ == "TRAJ20", CDR3a %in% str_remove(CDR3a, pattern = "\\*[^*]"))

pdf(paste0(fig_dir, "CDR3 size TRAV1-2-TRAJ20.pdf"), width=7, height=7)
p <- Total_TRAV_TRAJ %>% group_by(CDR3_Length) %>% summarise(N = n()) %>%
  ggplot(aes(x = CDR3_Length, y = N)) + geom_bar(stat="identity", width = 0.5, fill = "#08306b") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 20)) +
  scale_x_continuous(limits = c(8, 14), breaks = c(seq(from = 8, to = 14, by = 1))) +
  labs(x = "CDR3 Length (aa)", y = "Nb of Sequences") +
  theme_classic(base_size = 14) + My_Theme
print(p)
dev.off()

pdf(paste0(fig_dir, "CDR3 logo TRAV1-2-TRAJ20.pdf"), width = 10, height = 3)
p <- Total_TRAV_TRAJ %>%
  dplyr::filter(CDR3_Length == "11") %>% 
  dplyr::select(CDR3a) %>% ggseqlogo(seq_type = 'aa', method = 'bits') +
  theme(axis.text=element_text(size=28), axis.title=element_text(size=28, face="bold")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 28, b = 0, l = 0))) +
  theme(legend.position="none")
print(p)
dev.off()

Total_TRAV_TRAJ <- TCR_data_MAIT_cells.df %>%
  dplyr::select(- new_clusters_MAIT, -TRBV, -TRBD, -TRBJ, -CDR3b) %>%
  dplyr::filter(TRAV %in% c("TRAV1-1*01", "TRAV1-1*02",
                            "TRAV1-2*01")) %>%
  mutate(TRAV = str_remove(TRAV, pattern = "\\*[^.]*$"),
         TRAJ = str_remove(TRAJ, pattern = "\\*[^.]*$")) %>%
  mutate(CDR3_Length = str_length(CDR3a)) %>%
  dplyr::filter(TRAV == "TRAV1-2", TRAJ == "TRAJ12", CDR3a %in% str_remove(CDR3a, pattern = "\\*[^*]"))

pdf(paste0(fig_dir, "CDR3 size TRAV1-2-TRAJ12.pdf"), width=7, height=7)
p <- Total_TRAV_TRAJ %>% group_by(CDR3_Length) %>% summarise(N = n()) %>%
  ggplot(aes(x = CDR3_Length, y = N)) + geom_bar(stat="identity", width = 0.5, fill = "#08306b") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 30)) +
  scale_x_continuous(limits = c(8, 14), breaks = c(seq(from = 8, to = 14, by = 1))) +
  labs(x = "CDR3 Length (aa)", y = "Nb of Sequences") +
  theme_classic(base_size = 14) + My_Theme
print(p)
dev.off()

pdf(paste0(fig_dir, "CDR3 logo TRAV1-2-TRAJ12.pdf"), width = 10, height = 3)
p <- Total_TRAV_TRAJ %>%
  dplyr::filter(CDR3_Length == "11") %>% 
  dplyr::select(CDR3a) %>% ggseqlogo(seq_type = 'aa', method = 'bits') +
  theme(axis.text=element_text(size=28), axis.title=element_text(size=28, face="bold")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 28, b = 0, l = 0))) +
  theme(legend.position="none")
print(p)
dev.off()


Total_TRBV_TRBJ_distribution <- TCR_data_MAIT_cells.df %>%
  dplyr::filter(TRAV %in% c("TRAV1-1*01", "TRAV1-1*02",
                            "TRAV1-2*01")) %>%
  dplyr::select(-new_clusters_MAIT, -TRAV, -TRAJ, -TRBD, -CDR3a, -CDR3b) %>%
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
Vbeta_colors <- c("TRBV3" = "#3c7c3d" ,"TRBV4" = "#EAF1F4", "TRBV5" = "#93D4DB", "TRBV6" = "#99CCA7", "TRBV7" = "#57AA65", 
                  "TRBV10" = "#2c6b67", "TRBV11" = "#122c43", "TRBV15" = "#192813", 
                  "TRBV19" = "#C2947E", "TRBV20" = "#c2d6a4",
                  "TRBV21" = "#785158", "TRBV23" = "#757575", "TRBV24" = "#9cc184",
                  "TRBV25" = "#3D7D8F", "TRBV27" = "#956E73", "TRBV28" = "#1f5b25", "TRBV29" = "#c3d6ce",
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
               order = c(rep("TRBV3", 1), rep("TRBV4", 12), rep("TRBV5", 2), rep("TRBV6", 13), rep("TRBV7", 2), rep("TRBV10", 3),
                         rep("TRBV11", 1), rep("TRBV15", 5),
                         rep("TRBV19", 5), rep("TRBV20", 12), rep("TRBV21", 1), rep("TRBV23", 2), 
                         rep("TRBV24", 1), rep("TRBV25", 2), rep("TRBV27", 2), rep("TRBV28", 10), rep("TRBV29", 2),
                         rep("TRBJ1-1", 6),
                         rep("TRBJ1-2", 4), rep("TRBJ1-3", 3), rep("TRBJ1-4", 5), rep("TRBJ1-5", 4), rep("TRBJ1-6", 3), 
                         rep("TRBJ2-1", 11),
                         rep("TRBJ2-2", 6), rep("TRBJ2-3", 10), rep("TRBJ2-4", 2), rep("TRBJ2-5", 10), 
                         rep("TRBJ2-6", 6), rep("TRBJ2-7", 6)),
               col = sample(Vbeta_colors), transparency = 0.1, 
               link.lwd = 0.5, link.lty = 0.1, link.border = 1,
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(Total_TRBV_TRBJ_distribution))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.4))
  }, bg.border = NA) 
}

chordWrap(Total_TRBV_TRBJ_distribution)

########################
TCR_data_MAIT_cells_TCRb_diversity_thymus.df <- TCR_data_MAIT_cells.df %>%
  mutate(TRBV = str_remove(TRBV, pattern = "\\-[^.]*$"),
         TRBV = str_remove(TRBV, pattern = "\\*[^.]*$"),
         TRBJ = str_remove(TRBJ, pattern = "\\*[^.]*$")) %>%
  dplyr::filter(TRAV %in% c("TRAV1-1*01", "TRAV1-1*02",
                            "TRAV1-2*01")) %>%
  dplyr::select(-TRAV, -TRAJ, -TRBD, -CDR3a)

TCR_data_MAIT_cells_TCRb_diversity_thymus_naive.df <- TCR_data_MAIT_cells_TCRb_diversity_thymus.df %>%
  dplyr::filter(new_clusters_MAIT %in% c("MAIT_c2", "MAIT_c3", "MAIT_c4")) %>%
  mutate(TCRb_d_clonotype = paste(TRBV, TRBJ, CDR3b, sep = "_")) %>%
  dplyr::select(-new_clusters_MAIT) %>% group_by(TCRb_d_clonotype) %>%
  summarise(n = n()) %>% arrange(desc(n)) %>% pivot_wider(names_from = TCRb_d_clonotype, values_from = n)

shannon_iteration <- function(){
  TCR_data_MAIT_cells_TCRb_diversity_thymus_naive.df %>% 
    rrarefy(sample = 88) %>% 
    diversity(index = "shannon")
}

mean_shannon_thymus_naive <- mean(replicate(100, shannon_iteration()))

TCR_data_MAIT_cells_TCRb_diversity_thymus_C6.df <- TCR_data_MAIT_cells_TCRb_diversity_thymus.df %>%
  dplyr::filter(new_clusters_MAIT %in% c("MAIT_c6")) %>%
  mutate(TCRb_d_clonotype = paste(TRBV, TRBJ, CDR3b, sep = "_")) %>%
  dplyr::select(-new_clusters_MAIT) %>% group_by(TCRb_d_clonotype) %>%
  summarise(n = n()) %>% arrange(desc(n)) %>% pivot_wider(names_from = TCRb_d_clonotype, values_from = n)

shannon_iteration <- function(){
  TCR_data_MAIT_cells_TCRb_diversity_thymus_C6.df %>% 
    rrarefy(sample = 88) %>% 
    diversity(index = "shannon")
}

mean_shannon_thymus_c6 <- mean(replicate(100, shannon_iteration()))

shannon_diversity_NKT <- c(mean_shannon_thymus_naive, mean_shannon_thymus_c6)
samples <- c("Clusters 2/3/4", "Cluster 6")
df <- data.frame(samples, shannon_diversity_NKT)
df$samples <- factor(df$samples, levels = c("Clusters 2/3/4", "Cluster 6"))

plot <- ggplot(df, aes(x = samples, y = shannon_diversity_NKT)) + 
  geom_bar(stat="identity", width = 0.2, fill = "#08306b") +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(3, 5)) +
  labs(x = "Samples", y = "Shannon Index") +
  theme_classic(base_size = 14) + My_Theme

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/Shannon Index MAIT.pdf", 
       plot, device = "pdf", width=7, height=7)


saveRDS(thymus.MAIT, "/Volumes/Samsung_T5/Human_MAIT_NKT/Data/thymus.MAIT_03_02_23.RDS")

thymus.MAIT <- readRDS("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Data/thymus.MAIT_03_02_23.RDS")

