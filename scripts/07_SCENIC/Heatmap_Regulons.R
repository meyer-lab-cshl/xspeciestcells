# Load necessary libraries
library(ggplot2)
library(tidyverse)
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(cowplot)
library(scico)

# Read regulons data
regulons_path <- "/Volumes/Gapin-lab/For Tonya/new_data_target_gene_thresholds_for_scenic_10242023/auc_matrix/nes2_maskdropouts_full_obj_seurat_filtered_harmony_08_28_23_raw_counts_100_multirun_reg100-target95_regulons_pruned_final_GEP_INTERSECTION_TARGET_PRUNED_auc_mtx.csv"
regulons  <-  read.csv(regulons_path,  row.names  =  1)
regulons$cell_id  <-  rownames(regulons)
head(regulons, 5)

# Read and filter regulons to keep
regulons_to_keep_path <- "/Volumes/Gapin-lab/For Tonya/regulons_to_keep_min5_target.txt"
regulons_to_keep <- read_table(regulons_to_keep_path, col_names = FALSE)$X1

# Load processed Seurat object
seurat_object_path <- "/Volumes/Samsung_T5/Human_MAIT_NKT/Data/seurat_filtered_harmony_08_28_23.RDS"
filtered_seurat_harmony  <-  readRDS(seurat_object_path)

# Import GEP usage data
gep_usage_path <- "/Users/laurentgapin/Desktop/cNMF_non_imputed_all_cells/non_imputed_cNMF_allcells/non_imputed_cNMF_allcells.usages.k_12.dt_0_02.consensus.txt"
gep_usage <- read.table(gep_usage_path, header=TRUE)

# Rename GEP columns based on specified order
colnames(gep_usage) <- paste0("gep", c(2,5,3,1,4,12,6,7,8,10,9,11), "_usage")

# Exclude GEP12 (batch-driven)
gep_usage <- gep_usage[, !colnames(gep_usage) %in% "gep12_usage"]

# Assign cells to the GEP with the highest usage
df <- gep_usage
df$gep_assign <- gsub("_.*", "", colnames(df)[apply(df, 1, which.max)])
table(df$gep_assign)

# Add GEP assignment to Seurat metadata
filtered_seurat_harmony <- AddMetaData(filtered_seurat_harmony, df$gep_assign, col.name = "gep_assign")
filtered_seurat_harmony@meta.data$gep_assign <- factor(filtered_seurat_harmony@meta.data$gep_assign, levels = c("gep1",
                                                                                                                "gep2",
                                                                                                                "gep3",
                                                                                                                "gep4",
                                                                                                                "gep5",
                                                                                                                "gep6",
                                                                                                                "gep7",
                                                                                                                "gep8",
                                                                                                                "gep9",
                                                                                                                "gep10",
                                                                                                                "gep11"))
# Downsample Seurat object to 10,000 cells
downsampled_seurat <- subset(filtered_seurat_harmony, cells = sample(Cells(filtered_seurat_harmony), 10000))

# Extract metadata and join with regulons data
metadata <- downsampled_seurat@meta.data %>% 
  dplyr::select(cell_id, Sex, Donor, Batch, Tissue, new_clusters, new_clusters_id, gep_assign)
metadata_regulons <- left_join(metadata, regulons, by = "cell_id")

# Order and factorize data
metadata_regulons$new_clusters <- factor(metadata_regulons$new_clusters, levels = seq(0,17, 1))
metadata_regulons <- metadata_regulons[order(metadata_regulons$new_clusters, decreasing = FALSE),]
head(metadata_regulons, 3)

# Prepare data for analysis
data <- metadata_regulons %>% select(-Sex, -Donor, -Batch, -Tissue, -new_clusters_id, -new_clusters, -gep_assign)
rownames(data) <- data$cell_id
data$cell_id <- NULL
data <- t(data) # Transpose data
data[1:5, 1:5]

# Filter data based on regulons to keep
select_data <- as.data.frame(data) %>% filter(row.names(data) %in% regulons_to_keep)
dim(select_data)
select_data[1:5, 1:5]

metadata_df  <-  downsampled_seurat@meta.data
metadata_df  <-  metadata_df  %>%  dplyr::select(cell_id,  Sex,  Donor,  Batch,  Tissue,  new_clusters,  new_clusters_id, gep_assign)
metadata_df_regulons  <-  left_join(metadata_df,  regulons  ,  by  =  "cell_id")
metadata_df_regulons$new_clusters  <-  factor(metadata_df_regulons$new_clusters,  levels  =  seq(0,17,  1))
metadata_df_regulons  <-  metadata_df_regulons[order(metadata_df_regulons$new_clusters,  decreasing  =  F),]
dim(metadata_df_regulons)
metadata_df_regulons[1:3,  1:10]

data_2  <-  metadata_df_regulons  %>%  select(-Sex,  -Donor,  -Batch,  -Tissue, -new_clusters_id, -new_clusters, -gep_assign)
data_2[1:3,  1:7]
rownames(data_2) <- data_2$cell_id
data_2$cell_id <- NULL
data_2[1:3,  1:7]
data_2 <- t(data_2)
data_2[1:3,  1:7]

select_data_2  <-  as.data.frame(data_2)  %>%  filter(row.names(data_2)  %in%  regulons_to_keep)
select_data_2[1:5,  1:5]
dim(select_data_2)
select_data_2 <- t(select_data_2)
select_data_2 <- as.data.frame(select_data_2)
select_data_2$new_clusters_id <- metadata_df_regulons$new_clusters_id

# Calculate the 90th percentile for each column by group
percentiles_90 <- select_data_2 %>%
  select_if(is.numeric) %>%
  summarise(across(everything(), ~ quantile(., 0.90, na.rm = TRUE)))

result <- select_data_2 %>%
  group_by(new_clusters_id) %>%
  summarise(across(where(is.numeric), 
                   ~ mean(. > percentiles_90[[cur_column()]], na.rm = TRUE), 
                   .names = "prop_above_{.col}"))

# Filter the percentiles dataframe to keep only columns with at least one value > 10%
cols_to_keep <- result %>%
  select(-new_clusters_id) %>%
  summarise_all(~any(. > 0.20)) %>%
  as.data.frame() %>%
  names(.)[unlist(.)] %>% colnames() %>% gsub("prop_above_|_p90", "", .)

data  <-  metadata_df_regulons  %>%  select(-Sex,  -Donor,  -Batch,  -Tissue, -new_clusters_id, -new_clusters, -gep_assign)
data[1:3,  1:7]
rownames(data) <- data$cell_id
data$cell_id <- NULL
data[1:3,  1:7]
data  <-  t(data) 

select_data_to_plot  <-  as.data.frame(data)  %>%  filter(row.names(data)  %in%  cols_to_keep)
dim(select_data_to_plot)
#check  the  number  of  columns  with  regulons
select_data_to_plot[1:5,  1:5]

# Set up colors and other aesthetics for heatmap
heat_colors <- scico(100, palette = "vik")
colors_clusters <- c("0"  =  "#f4c40f",  "1"  =  "#b75347",  "2"  =  "#d8443c",  "3"  =  "#e09351",  "4"  =  "#2b9b81",  
                     "5"  =  "#421401",  "6"  =  "#92c051",  "7"  =  "#9f5691",  "8"  =  "#17154f",  "9"  =  "#74c8c3",  
                     "10"  =  "#5a97c1",  "11"  =  "gold",  "12"  =  "#a40000",  "13"  =  "#72bcd5",  "14"  =  "grey50",
                     "15"  =  "orange",  "16"  =  "blueviolet",  "17"  =  "#0a2e57") 
names(colors_clusters) <- seq(0, 17)

# Prepare metadata for annotation
meta <- metadata_df_regulons[,1:8]
meta <- meta %>% filter(!new_clusters == 18)
meta$new_clusters <- factor(meta$new_clusters, seq(0, 17, 1))
rownames(meta) <- meta$cell_id

# Prepare GEP assign and tissue colors
heat_colors  <-  scico(100,  alpha  =  NULL,  begin  =  0,  end  =  1,  direction  =  1,  palette  =  "vik",  categorical  =  FALSE)

colors_clusters  <-  c("0"  =  "#f4c40f",  "1"  =  "#b75347",  "2"  =  "#d8443c",  "3"  =  "#e09351",  "4"  =  "#2b9b81",  
                       "5"  =  "#421401",  "6"  =  "#92c051",  "7"  =  "#9f5691",  "8"  =  "#17154f",  "9"  =  "#74c8c3",  
                       "10"  =  "#5a97c1",  "11"  =  "gold",  "12"  =  "#a40000",  "13"  =  "#72bcd5",  "14"  =  "grey50",
                       "15"  =  "orange",  "16"  =  "blueviolet",  "17"  =  "#0a2e57")

colors_tissues <- c("#72bcd5", "#a40000")
names(colors_tissues) <- unique(meta$Tissue)

colors_GEPs <-  c("gep1"  =  "#f5bb50",  "gep2"  =  "#ada43b",  "gep3"  =  "#b0799a",  "gep4"  =  "#f6b3b0",  "gep5"  =  "#bf3729",  
                  "gep6"  =  "#17154f",  "gep7"  =  "#355828",  "gep8"  =  "#e48171",  "gep9"  =  "#2f357c",  "gep10"  =  "#6c5d9e",  
                  "gep11"  =  "#e69b00")


# Reorder rows in the data based on specific order
the_order <- c("PBX1", "ETV6", "HOXA10", "MTF2", "ZNF711", "ETS2", "HES1", "FOXO6", "HOXA3",
               "SOX4", "ILF2", "TFDP2", "YBX1", "NFYB", "MYB", "NONO", "SFPQ", "CHD1", "NR3C1", "HDAC2",
               "HMGA1", "FOXM1", "TFDP1", "E2F1", "MYBL2", "E2F8", "E2F2", "E2F7", "ZNF69",
               "BCL11A", "STAT5A", "KLF5" , "VDR", "IRF8",
               "BACH2", "IRF4", "NFKB1", "REL", "ZFX", "ZNF333", "RFX3", "VEZF1", "EGR2", "EGR3", "EGR1","SMAD3", "STAT6",
               "TCF3", "ETV5", "GATA3", "RAD21", "ELF1", "IKZF2", "BCL11B", "BCL6", "LEF1", "NFATC1", "SREBF2", "HOXB2", "PPARG", 
               "NR2C2", "HMGXB4", "ELK1", "RFX5", 
               "RXRA", "IKZF1", "YY1", "BPTF", "REST", "EP300", "MAX", "FOSL1", "GTF2I",
               "CEBPA", "FOXP3", "NR2F1", "ZNF610", "IRF6", 
               "IRF2" , "ELK4", "IRF7", "STAT2", 
               "FOXO1", "FOXP1", "STAT1", "KLF4",
               "BCLAF1", "KLF9", "ETS1", "FOS",  "JUNB", "JUN", "FOSB",  "JUND", 
               "ELK3", "MBD2", "CREM", "NFE2L2", "NR1D2", "XBP1", "MYBL1", "RORA", "MAF", "CEBPD", "FOSL2", "EOMES",
               "IRF1" , "KLF2",  "KLF6", "RUNX3", "PRDM1", "FLI1", "KLF12", "IKZF3", "ZBTB20", 
               "ELF4", "ZBTB7A", "NFATC2", "TBX21",  "NFATC3", "ZBTB44", "KLF3", 
               "ATF4", "E2F3", "KLF13", "ETV7", "HCLS1", "HES7", "HNF1B", "PURA", "SIN3A",
               "GATA6", "MTF1",
               "IRF5", 
               "VPS72" , "ZBTB25",
               "CLOCK", "ELF2", "GATAD1", "KLF8", "MYC", "NFIC", "SP1", "SP2", "SP3", "TAF6", "ZBTB11", "ZFY", "ZNF460")

select_data_to_plot <- select_data_to_plot[the_order,] 

column_cuts <- which(duplicated(metadata_df_regulons$new_clusters) == FALSE)
column_cuts <- column_cuts[2:length(column_cuts)] - 1
# Values to be removed
values_to_remove <- c(9542, 9361, 7267)

# Remove the specified values
column_cuts <- column_cuts[!column_cuts %in% values_to_remove]

# Create heatmap
heatmap <- pheatmap::pheatmap(select_data_to_plot,
                              color = heat_colors,  
                              breaks = seq(-5, 5, by = 0.1),
                              scale = "row",
                              cellheight = 6,
                              cluster_rows = FALSE,
                              cluster_cols = FALSE,
                              show_colnames = FALSE,
                              gaps_row = c(9, 14, 19, 28, 38, 41, 46, 60, 64, 73, 77, 82, 88, 94, 105, 122),# Specify gaps in rows
                              gaps_col = column_cuts, 
                              annotation_col = meta %>% dplyr::select(c(Tissue, new_clusters, gep_assign)),
                              annotation_colors = list(new_clusters = colors_clusters, gep_assign = colors_GEPs, Tissue = colors_tissues),
                              fontsize_row = 6, width = 10, height = 10,
                              main = "Scenic Run with Full Data \n ± 10Kb \n with AUC scores row scaled")




