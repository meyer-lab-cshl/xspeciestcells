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

blood.GD <- preprocess(seurobj = subset(blood.filt, subset=group.ident=="GD_PBMC"),
                         celltype = "GDblood", colvalues=colvalues, ndim=10,
                         res = 0.3)

DimPlot(blood.GD, group.by="RNA_snn_res.0.3", label=T, repel=T) +
  theme(legend.position="none") + 
  labs(title="Blood GD")

SCpubr::do_FeaturePlot(sample = blood.GD, 
                       features = c("GZMK", "GZMB"), order = T,
                       plot.title = "", legend.position = "right",
                       viridis_color_map = "inferno")

blood.GD$new_clusters_GD <- case_when(
  blood.GD$RNA_snn_res.0.3 == '2' ~ '0',
  blood.GD$RNA_snn_res.0.3 == '0' ~ '1',
  blood.GD$RNA_snn_res.0.3 == '3' ~ '2',
  blood.GD$RNA_snn_res.0.3 == '1' ~ '3',
  blood.GD$RNA_snn_res.0.3 == '4' ~ '4',
)

blood.GD$new_clusters_GD  <- as.factor(blood.GD$new_clusters_GD)
Idents(blood.GD) <- "new_clusters_GD"

colors <- paletteer_c(`"ggthemes::Orange-Blue Diverging"`, n = 15)
colors_clusters_GD <- c("0" = "#DF6D27FF", "1" = "grey40", "2" = "#AB6969", "3" = "#a40000", "4" = "gold")

DimPlot(blood.GD, cols = colors_clusters_GD, pt.size = 1, label = T)

plot <- SCpubr::do_DimPlot(sample = blood.GD, 
                           label = TRUE, repel = T, font.size = 20,
                           label.color = "white", colors.use = colors_clusters_GD,
                           legend.position = "right")

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/GD_analysis_clusters_SCPubR_blood.pdf",
       plot, device = "pdf", width=7, height=7, limitsize = F)

blood.GD <- NormalizeData(blood.GD)

GEPs <- read.csv("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Data/cnmf_NKT_thymus/NKT_thymus_cNMF/TopGenes_perGEP_ALL.cells.csv")

GEP_4 <- GEPs %>% select(GEP_4)
GEP_4 <- as.vector(GEP_4$GEP_4)
GEP_4_score <- list(GEP_4)

blood.GD <- AddModuleScore(object = blood.GD, features = GEP_4_score,
                             assay = "RNA", name = 'GEP_4_score')

GEP_7 <- GEPs %>% select(GEP_7)
GEP_7 <- as.vector(GEP_7$GEP_7)
GEP_7_score <- list(GEP_7)

blood.GD <- AddModuleScore(object = blood.GD, features = GEP_7_score,
                             assay = "RNA", name = 'GEP_7_score')

GEP_6 <- GEPs %>% select(GEP_6)
GEP_6 <- as.vector(GEP_6$GEP_6)
GEP_6_score <- list(GEP_6)

blood.GD <- AddModuleScore(object = blood.GD, features = GEP_6_score,
                             assay = "RNA", name = 'GEP_6_score')

plot <- SCpubr::do_FeaturePlot(sample = blood.GD, 
                               features = c("GEP_6_score1","GEP_4_score1", "GEP_7_score1"), order = T,
                               plot.title = "", legend.position = "right", ncol = 3,
                               viridis_color_map = "inferno")

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/GEP_signature_GD_blood.pdf",
       plot, device = "pdf", width=21, height=7, limitsize = F)

###Distribution of cells per Donor####################################################################################
Idents(blood.GD) <- "new_clusters_GD"

# Extract and summarize NKT cell numbers by donor and new_clusters_NKT
GD_cell_numbers_blood <- blood.GD@meta.data %>%
  as.data.frame()
GD_cell_numbers_blood <- GD_cell_numbers_blood %>%
  dplyr::select(Donor, new_clusters_GD) %>%
  dplyr::count(new_clusters_GD, Donor) %>%
  dplyr::rename(N = n) %>%
  ungroup()

# Compute frequency of NKT cells for each donor and new_clusters_NKT
GD_cell_numbers_blood_test <- GD_cell_numbers_blood %>%
  group_by(Donor) %>%
  mutate(freq = N/sum(N)*100) %>%
  ungroup()

# Combine frequency data for all donors and new_clusters_NKT
distribution_donors_blood <- GD_cell_numbers_blood_test %>%
  complete(new_clusters_GD, Donor, fill=list(N=0, freq=0)) %>%
  mutate(Donor = factor(Donor))

plot <- ggplot(data= subset(distribution_donors_blood, !is.na(Donor)), aes(x = Donor, y = freq, 
                                                                           fill = forcats::fct_rev(factor(new_clusters_GD)))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = colors_clusters_GD, 
                    labels = levels(distribution_donors_blood$new_clusters_GD), 
                    drop = FALSE, limits = levels(distribution_donors_blood$new_clusters_GD),
                    guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Donors", y = "%", fill = "Clusters") + theme_classic() + My_Theme 

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/Distribution_Blood_GD_per_donor.pdf",
       plot, device = "pdf", width=7, height=7, limitsize = F)

######################
plot <- SCpubr::do_FeaturePlot(sample = blood.GD, 
                               features = c("CD4", "CD8A"), order = T,
                               plot.title = "", legend.position = "right",
                               viridis_color_map = "inferno")

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/CD4_CD8A_Blood_GD.pdf",
       plot, device = "pdf", width=14, height=7, limitsize = F)

##DEGs#################################################################################################################
Idents(blood.GD) <- "new_clusters_GD"
GD.clusters.markers <- FindAllMarkers(blood.GD, test.use = 'wilcox', 
                                        logfc.threshold = 0.4, min.pct = 0.3, only.pos = TRUE)
top5 <- GD.clusters.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5_distinct <- base::unique(top5$gene)

write.csv(GD.clusters.markers, "/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/clusters.markers.GD.blood.csv")

GD_blood_samples <- ScaleData(blood.GD, features = top5_distinct)
DotPlot_colors <- c("lightsteelblue1", "#1E466E")

plot <- DotPlot(GD_blood_samples, features = top5_distinct, dot.scale = 8, cols = DotPlot_colors,
                col.min = -1, col.max = 2, dot.min = 0) + RotatedAxis()
ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/Features_per_GD_clusters.blood.pdf",
       plot, device = "pdf", width=10, height=5, limitsize = F)

##################
saveRDS(blood.GD, "/Volumes/Samsung_T5 1/Human_MAIT_NKT/Data/blood.GD_03_16_23.RDS")

#################

TCR_data_GD_cells.df <- as.data.frame(blood.GD@meta.data)
TCR_data_GD_cells.df <- TCR_data_GD_cells.df %>%
  dplyr::select(-nCount_RNA, -nFeature_RNA, -percent.mt, -cell.ident, -Sex, -Age_in_weeks, -Donor, -Batch, -Method,
                -RNA_snn_res.0.7, -seurat_clusters, -orig.ident, -Tissue, -Test_cell_ident_for_cLISI,
                -group.ident, -TCRa_g_chain, -TCRa_g_clonotype, -TCRb_d_chain, -TCRb_d_clonotype, 
                -TRAV10_TRAJ18, -TRAV1_TRAJ33, -TRAV1, -new_clusters, -Effectorness1, -Naiveness1, -Type_3_score1, -new_clusters_id) %>%
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

plot <- SCpubr::do_DimPlot(sample = blood.GD, 
                           label = FALSE,
                           cells.highlight = TRDV1_cells,
                           legend.position = "none",
                           na.value = "grey90",
                           sizes.highlight = 3, colors.use = "#a40000")

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/TRDV1_GD.blood.pdf", 
       plot, device = "pdf", width=7, height=7)

plot <- SCpubr::do_DimPlot(sample = blood.GD, 
                           label = FALSE,
                           cells.highlight = TRDV2_cells,
                           legend.position = "none",
                           na.value = "grey90",
                           sizes.highlight = 3, colors.use = "blue")

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/TRDV2_GD.blood.pdf", 
       plot, device = "pdf", width=7, height=7)

plot <- SCpubr::do_DimPlot(sample = blood.GD, 
                           label = FALSE,
                           cells.highlight = TRDV3_cells,
                           legend.position = "none",
                           na.value = "grey90",
                           sizes.highlight = 3, colors.use = "#318f49")

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/TRDV3_GD.blood.pdf", 
       plot, device = "pdf", width=7, height=7)

table(TCR_data_GD_cells.df$TRGV)
TRGV9_cells <- TCR_data_GD_cells.df %>% dplyr::filter(TRGV %in% c("TRGV9*01")) %>% row.names()

plot <- SCpubr::do_DimPlot(sample = blood.GD, 
                           label = FALSE,
                           cells.highlight = TRGV9_cells, 
                           legend.position = "none",
                           na.value = "grey90",
                           sizes.highlight = 3, colors.use = "gold")

ggsave("/Volumes/Samsung_T5/Human_MAIT_NKT/Final_Figures/TRGV9_GD.pdf", 
       plot, device = "pdf", width=7, height=7)
