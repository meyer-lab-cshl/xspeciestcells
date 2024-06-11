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
all_cells <- readRDS("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Data/seurat_filtered_harmony_02_15_23.RDS")
########################

########Analysis
# filter for thymus only:
blood <- subset(all_cells, subset=Tissue=="PBMC")
print(blood) # 41238 cells

# Make raw object
blood.filt <- CreateSeuratObject(counts=blood@assays$RNA@counts,
                                 meta.data = blood@meta.data)

print(blood.filt) # 37,369 cells

blood.CD8 <- preprocess(seurobj = subset(blood.filt, subset=group.ident=="CD8_PBMC"),
                        celltype = "CD4blood", colvalues=colvalues, ndim=10,
                        res = 0.5)

DimPlot(blood.CD8, group.by="RNA_snn_res.0.5", label=T, repel=T) +
  theme(legend.position="none") + 
  labs(title="Blood CD8")

SCpubr::do_FeaturePlot(sample = blood.CD8, 
                       features = c("GZMK", "GZMB"), order = T,
                       plot.title = "", legend.position = "right",
                       viridis_color_map = "inferno")

blood.CD8$new_clusters_CD8 <- case_when(
  blood.CD8$RNA_snn_res.0.5 == '0' ~ '0',
  blood.CD8$RNA_snn_res.0.5 == '3' ~ '1',
  blood.CD8$RNA_snn_res.0.5 == '4' ~ '2',
  blood.CD8$RNA_snn_res.0.5 == '2' ~ '3',
  blood.CD8$RNA_snn_res.0.5 == '1' ~ '4',
  blood.CD8$RNA_snn_res.0.5 == '5' ~ '5'
)

blood.CD8$new_clusters_CD8  <- as.factor(blood.CD8$new_clusters_CD8)
Idents(blood.CD8) <- "new_clusters_CD8"

colors <- paletteer_c(`"ggthemes::Orange-Blue Diverging"`, n = 15)
colors_clusters_CD8 <- c("0" = "#DF6D27FF", "1" = "#E9BE99FF", "2" = "grey40", "3" = "grey70", "4" = "#a40000",
                         "5" = "gold")

#2B5C8AFF
#77ACD3FF

plot <- SCpubr::do_DimPlot(sample = blood.CD8, 
                           label = FALSE, repel = T, #font.size = 20,
                           label.color = "white", colors.use = colors_clusters_CD8,
                           reduction = "umap", pt.size = 7, raster = T, raster.dpi = 2048,
                           legend.position = "right")

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/CD8_analysis_clusters_SCPubR_blood.pdf",
       plot, device = "pdf", width=7, height=7, limitsize = F)

GEPs <- read.csv("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Data/genes_per_gep_post_rank_threshold.csv")

GEP_6 <- GEPs %>% dplyr::select(GEP6) %>% na.omit()
GEP_6 <- as.vector(GEP_6$GEP6)
GEP_6_score <- list(GEP_6) 

blood.CD8 <- AddModuleScore(object = blood.CD8, features = GEP_6_score,
                            assay = "RNA", name = 'GEP_6_score')

GEP_1 <- GEPs %>% dplyr::select(GEP1) %>% na.omit()
GEP_1 <- as.vector(GEP_1$GEP1)
GEP_1_score <- list(GEP_1)

blood.CD8 <- AddModuleScore(object = blood.CD8, features = GEP_1_score,
                            assay = "RNA", name = 'GEP_1_score')

GEP_4 <- GEPs %>% dplyr::select(GEP4) %>% na.omit()
GEP_4 <- as.vector(GEP_4$GEP4)
GEP_4_score <- list(GEP_4)

blood.CD8 <- AddModuleScore(object = blood.CD8, features = GEP_4_score,
                            assay = "RNA", name = 'GEP_4_score')

plot <- SCpubr::do_FeaturePlot(sample = blood.CD8, 
                               features = c("GEP_6_score1", "GEP_1_score1", "GEP_4_score1"), order = T,
                               plot.title = "", legend.position = "right", ncol = 4,
                               min.cutoff = c(0, 0, 0, 0), max.cutoff = c(0.08, 0.08, 0.08, 0.08),
                               reduction = "umap", pt.size = 10, raster = T, raster.dpi = 2048,
                               viridis_color_map = "inferno")

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/GEP_signature_CD8_blood.pdf",
       plot, device = "pdf", width=21, height=7, limitsize = F)

plot <- SCpubr::do_ViolinPlot(sample = blood.CD8, 
                              features = c("GEP_6_score1", "GEP_1_score1", "GEP_4_score1"),
                              share.y.lims = TRUE, ncol = 3, colors.use = colors_clusters_CD8)

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/GEP_signature_CD8_blood_Violin.pdf",
       plot, device = "pdf", width=15, height=7, limitsize = F)

###Distribution of cells per Donor####################################################################################
Idents(blood.CD8) <- "new_clusters_CD8"

# Extract and summarize NKT cell numbers by donor and new_clusters_NKT
CD8_cell_numbers_blood <- blood.CD8@meta.data %>%
  as.data.frame()
CD8_cell_numbers_blood <- CD8_cell_numbers_blood %>%
  dplyr::select(Donor, new_clusters_CD8) %>%
  dplyr::count(new_clusters_CD8, Donor) %>%
  dplyr::rename(N = n) %>%
  ungroup()

# Compute frequency of NKT cells for each donor and new_clusters_NKT
CD8_cell_numbers_blood_test <- CD8_cell_numbers_blood %>%
  group_by(Donor) %>%
  mutate(freq = N/sum(N)*100) %>%
  ungroup()

# Combine frequency data for all donors and new_clusters_NKT
distribution_donors_blood <- CD8_cell_numbers_blood_test %>%
  complete(new_clusters_CD8, Donor, fill=list(N=0, freq=0)) %>%
  mutate(Donor = factor(Donor))

plot <- ggplot(data= subset(distribution_donors_blood, !is.na(Donor)), aes(x = Donor, y = freq, 
                                                                           fill = forcats::fct_rev(factor(new_clusters_CD8)))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = colors_clusters_CD8, 
                    labels = levels(distribution_donors_blood$new_clusters_CD8), 
                    drop = FALSE, limits = levels(distribution_donors_blood$new_clusters_CD8),
                    guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Donors", y = "%", fill = "Clusters") + theme_classic() + My_Theme 

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/Distribution_Blood_CD8_per_donor.pdf",
       plot, device = "pdf", width=7, height=7, limitsize = F)

######################
plot <- SCpubr::do_FeaturePlot(sample = blood.CD8, 
                               features = c("CCR7", "LEF1"), order = T,
                               plot.title = "", legend.position = "right",
                               viridis_color_map = "inferno")

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/CCR7_LEF1_Blood_CD8.pdf",
       plot, device = "pdf", width=14, height=7, limitsize = F)

##DEGs#################################################################################################################
Idents(blood.CD8) <- "new_clusters_CD8"
CD8.clusters.markers <- FindAllMarkers(blood.CD8, test.use = 'wilcox', 
                                        logfc.threshold = 0.4, min.pct = 0.3, only.pos = TRUE)
top5 <- CD8.clusters.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5_distinct <- base::unique(top5$gene)

write.csv(CD8.clusters.markers, "/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/clusters.markers.CD8.blood.csv")

CD8_blood_samples <- ScaleData(blood.CD8, features = top5_distinct)
DotPlot_colors <- c("lightsteelblue1", "#1E466E")

plot <- DotPlot(CD8_blood_samples, features = top5_distinct, dot.scale = 8, cols = DotPlot_colors,
                col.min = -1, col.max = 2, dot.min = 0) + RotatedAxis()
ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/Features_per_CD8_clusters.blood.pdf",
       plot, device = "pdf", width=10, height=5, limitsize = F)

##################
saveRDS(blood.CD8, "/Volumes/Samsung_T5 1/Human_MAIT_NKT/Data/blood.CD8_03_16_23.RDS")

blood.CD8 <- readRDS("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Data/blood.CD8_03_16_23.RDS")

##################
TCR_data_CD8_cells.df <- as.data.frame(blood.CD8@meta.data)
TCR_data_CD8_cells.df <- TCR_data_CD8_cells.df %>%
  dplyr::select(-nCount_RNA, -nFeature_RNA, -percent.mt, -cell.ident, -Sex, -Age_in_weeks, -Donor, -Batch, -Method,
                -seurat_clusters, -orig.ident, -Tissue, -Test_cell_ident_for_cLISI,-RNA_snn_res.0.7, -RNA_snn_res.1.2, -Egress_score1,
                -Effectorness1, -Naiveness1, -Type_3_score1, -new_clusters_id,
                -group.ident, -TCRa_g_chain, -TCRa_g_clonotype, -TCRb_d_chain, -TCRb_d_clonotype, -TRAV10_TRAJ18, -TRAV1_TRAJ33, -TRAV1, -new_clusters) %>%
  na.omit() %>%
  dplyr::rename(TRAV = TCR_Alpha_Gamma_V_gene_Dominant, TRAJ = TCR_Alpha_Gamma_J_gene_Dominant, CDR3a = TCR_Alpha_Gamma_CDR3_Translation_Dominant,
                TRBV = TCR_Beta_Delta_V_gene_Dominant, TRBD = TCR_Beta_Delta_D_gene_Dominant, TRBJ = TCR_Beta_Delta_J_gene_Dominant,
                CDR3b = TCR_Beta_Delta_CDR3_Translation_Dominant) %>%
  mutate(CDR3b = paste0("C", CDR3b)) %>%
  mutate(CDR3a = paste0("C", CDR3a))

Total_TRAV_TRAJ_distribution <- TCR_data_CD8_cells.df %>%
  dplyr::select(-CDR3a, -TRBV, -TRBD, -TRBJ, -CDR3b) %>%
  mutate(TRAV = str_remove(TRAV, pattern = "\\*[^.]*$"),
         TRAJ = str_remove(TRAJ, pattern = "\\*[^.]*$")) %>%
  dplyr::filter(TRAV == "TRAV1-2") %>%
  group_by(TRAV, TRAJ) %>%
  summarise(n = n()) %>% ungroup() %>% as.data.frame()

table(TCR_data_CD8_cells.df$TRAV)
Trav1_cells <- TCR_data_CD8_cells.df %>% dplyr::filter(TRAV %in% c("TRAV1-1*01", "TRAV1-1*02",
                                                                    "TRAV1-2*01")) %>% row.names()
length(Trav1_cells)

plot <- SCpubr::do_DimPlot(sample = blood.CD8, 
                           label = FALSE,
                           cells.highlight = Trav1_cells,
                           plot.title = "", legend.position = "right",
                           na.value = "grey90", 
                           sizes.highlight = 2, colors.use = "#a40000")

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/TRAV1_cells.CD8.blood.pdf", 
       plot, device = "pdf", width=8, height=7)

