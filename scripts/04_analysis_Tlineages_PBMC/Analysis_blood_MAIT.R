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

blood.MAIT <- preprocess(seurobj = subset(blood.filt, subset=group.ident=="MAIT_PBMC"),
                        celltype = "MAITblood", colvalues=colvalues, ndim=10,
                        res = 0.4)

DimPlot(blood.MAIT, group.by="RNA_snn_res.0.4", label=T, repel=T) +
  theme(legend.position="none") + 
  labs(title="Blood MAIT")

SCpubr::do_FeaturePlot(sample = blood.MAIT, 
                       features = c("GZMK", "GZMB"), order = T,
                       plot.title = "", legend.position = "right",
                       viridis_color_map = "inferno")

blood.MAIT$new_clusters_MAIT <- case_when(
  blood.MAIT$RNA_snn_res.0.4 == '0' ~ '0',
  blood.MAIT$RNA_snn_res.0.4 == '1' ~ '1',
  blood.MAIT$RNA_snn_res.0.4 == '2' ~ '2',
  blood.MAIT$RNA_snn_res.0.4 == '3' ~ '3'
)

blood.MAIT$new_clusters_MAIT  <- as.factor(blood.MAIT$new_clusters_MAIT)
Idents(blood.MAIT) <- "new_clusters_MAIT"

colors <- paletteer_c(`"ggthemes::Orange-Blue Diverging"`, n = 15)
colors_clusters_MAIT <- c("0" = "grey90", "1" = "grey70", "2" = "grey40", "3" = "#a40000")

DimPlot(blood.MAIT, cols = colors_clusters_MAIT, pt.size = 1, label = T)

plot <- SCpubr::do_DimPlot(sample = blood.MAIT, 
                           label = TRUE, repel = T, font.size = 20,
                           label.color = "white", colors.use = colors_clusters_MAIT,
                           legend.position = "right")

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/MAIT_analysis_clusters_SCPubR_blood.pdf",
       plot, device = "pdf", width=7, height=7, limitsize = F)

blood.MAIT <- NormalizeData(blood.MAIT)

GEPs <- read.csv("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Data/cnmf_NKT_thymus/NKT_thymus_cNMF/TopGenes_perGEP_ALL.cells.csv")

GEP_4 <- GEPs %>% select(GEP_4)
GEP_4 <- as.vector(GEP_4$GEP_4)
GEP_4_score <- list(GEP_4)

blood.MAIT <- AddModuleScore(object = blood.MAIT, features = GEP_4_score,
                            assay = "RNA", name = 'GEP_4_score')

GEP_7 <- GEPs %>% select(GEP_7)
GEP_7 <- as.vector(GEP_7$GEP_7)
GEP_7_score <- list(GEP_7)

blood.MAIT <- AddModuleScore(object = blood.MAIT, features = GEP_7_score,
                            assay = "RNA", name = 'GEP_7_score')

GEP_6 <- GEPs %>% select(GEP_6)
GEP_6 <- as.vector(GEP_6$GEP_6)
GEP_6_score <- list(GEP_6)

blood.MAIT <- AddModuleScore(object = blood.MAIT, features = GEP_6_score,
                            assay = "RNA", name = 'GEP_6_score')

plot <- SCpubr::do_FeaturePlot(sample = blood.MAIT, 
                               features = c("GEP_6_score1","GEP_4_score1", "GEP_7_score1"), order = T,
                               plot.title = "", legend.position = "right", ncol = 3,
                               viridis_color_map = "inferno")

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/GEP_signature_MAIT_blood.pdf",
       plot, device = "pdf", width=21, height=7, limitsize = F)

###Distribution of cells per Donor####################################################################################
Idents(blood.MAIT) <- "new_clusters_MAIT"

# Extract and summarize NKT cell numbers by donor and new_clusters_NKT
MAIT_cell_numbers_blood <- blood.MAIT@meta.data %>%
  as.data.frame()
MAIT_cell_numbers_blood <- MAIT_cell_numbers_blood %>%
  dplyr::select(Donor, new_clusters_MAIT) %>%
  dplyr::count(new_clusters_MAIT, Donor) %>%
  dplyr::rename(N = n) %>%
  ungroup()

# Compute frequency of NKT cells for each donor and new_clusters_NKT
MAIT_cell_numbers_blood_test <- MAIT_cell_numbers_blood %>%
  group_by(Donor) %>%
  mutate(freq = N/sum(N)*100) %>%
  ungroup()

# Combine frequency data for all donors and new_clusters_NKT
distribution_donors_blood <- MAIT_cell_numbers_blood_test %>%
  complete(new_clusters_MAIT, Donor, fill=list(N=0, freq=0)) %>%
  mutate(Donor = factor(Donor))

plot <- ggplot(data= subset(distribution_donors_blood, !is.na(Donor)), aes(x = Donor, y = freq, 
                                                                           fill = forcats::fct_rev(factor(new_clusters_MAIT)))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = colors_clusters_NKT, 
                    labels = levels(distribution_donors_blood$new_clusters_MAIT), 
                    drop = FALSE, limits = levels(distribution_donors_blood$new_clusters_MAIT),
                    guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Donors", y = "%", fill = "Clusters") + theme_classic() + My_Theme 

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/Distribution_Blood_MAIT_per_donor.pdf",
       plot, device = "pdf", width=7, height=7, limitsize = F)

######################
plot <- SCpubr::do_FeaturePlot(sample = blood.nkt, 
                               features = c("CD4", "CD8A"), order = T,
                               plot.title = "", legend.position = "right",
                               viridis_color_map = "inferno")

ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/CD4_CD8A_Blood_NKT.pdf",
       plot, device = "pdf", width=14, height=7, limitsize = F)

##DEGs#################################################################################################################
Idents(blood.MAIT) <- "new_clusters_MAIT"
MAIT.clusters.markers <- FindAllMarkers(blood.MAIT, test.use = 'wilcox', 
                                       logfc.threshold = 0.4, min.pct = 0.3, only.pos = TRUE)
top5 <- MAIT.clusters.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5_distinct <- base::unique(top5$gene)

write.csv(MAIT.clusters.markers, "/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/clusters.markers.MAIT.blood.csv")

MAIT_blood_samples <- ScaleData(blood.MAIT, features = top5_distinct)
DotPlot_colors <- c("lightsteelblue1", "#1E466E")

plot <- DotPlot(MAIT_blood_samples, features = top5_distinct, dot.scale = 8, cols = DotPlot_colors,
                col.min = -1, col.max = 2, dot.min = 0) + RotatedAxis()
ggsave("/Volumes/Samsung_T5 1/Human_MAIT_NKT/Final_Figures/Features_per_MAIT_clusters.blood.pdf",
       plot, device = "pdf", width=10, height=5, limitsize = F)

##################

saveRDS(blood.MAIT, "/Volumes/Samsung_T5 1/Human_MAIT_NKT/Data/blood.MAIT_03_16_23.RDS")
