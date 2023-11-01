# Purpose:
# Author: Salomé Carcy
# Date:




# **************
# 1. IMPORT ####
# **************

# Import librairies
library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)
library(Seurat)
library(scico)
source("./scripts-final/colors_universal.R")

# Import human data
seur.nkt  <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.nkt.RDS")
seur.mait <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.mait.RDS")
seur.gdt  <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.gd.RDS")
seur.cd4  <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.cd4.RDS")
seur.cd8  <- readRDS("./data/human-thymus/HumanThymus_23_PlotThymicGEPs/seurat_filtered_harmony_02_15_23_thymus.cd8.RDS")

# Plot DimPlot
DimPlot(seur.nkt,  group.by="cell_annot", label=T)
DimPlot(seur.mait, group.by="cell_annot", label=T)
DimPlot(seur.gdt,  group.by="cell_annot", label=T)
DimPlot(seur.cd4,  group.by="cell_annot", label=T)
DimPlot(seur.cd8,  group.by="cell_annot", label=T)

# Cleanup a bit
seur.nkt@meta.data[,33:75]  <- NULL
seur.mait@meta.data[,33:75] <- NULL
seur.gdt@meta.data[,33:75]  <- NULL
seur.cd4@meta.data[,33:75]  <- NULL
seur.cd8@meta.data[,33:75]  <- NULL


# Import mouse data
seur.nkt.ms <- readRDS("./data/cross-species/00_Reproduce_UMAPs/ms_nkt_seurobj.rds")
seur.mait.ms <- readRDS("./data/cross-species/00_Reproduce_UMAPs/ms_mait_seurobj.rds")
ortholog.df <- read.csv("./data/cross-species/03_BiomartTable/big_ass_ortholog_table.csv")
DimPlot(seur.nkt.ms, group.by="cell_type", label=T)
DimPlot(seur.mait.ms, group.by="cell_type", label=T)



# *****************
# 2. FUNCTIONS ####
# *****************

plot_four_genes <- function(seur, genelist, ordercells=F){
  plist <- list()
  # identify which gene has highest max count
  max_per_gene <- lapply(genelist, function(x) max(seur@assays$RNA@data[x,]))
  names(max_per_gene) <- genelist
  gene_with_max_value <- names(which.max(max_per_gene))

  for (gene in genelist){
    print(gene)
    if(gene=="ZBTB16"){ordercells=T} # put PLZF at front
    # get legend once
    if(gene==gene_with_max_value){
      plegend <- ggpubr::get_legend(
        SCpubr::do_FeaturePlot(seur, features=gene, ncol=1, legend.position="right", legend.title="", border.size=1, pt.size=2, order=T) &
          scale_colour_scico(palette = "lapaz", alpha = 0.8, direction = -1)
      )
    }
    # featureplot
    p <- SCpubr::do_FeaturePlot(seur, features=gene, order=ordercells, ncol=1, legend.position="none", border.size=1, pt.size=3) +
      labs(title=gene)+
      theme(plot.background = element_rect(color = "black"),
            plot.title = element_text(hjust = 0.5, vjust=0.1, size=20)) &
      scale_colour_scico(palette = "lapaz", alpha = 0.8, direction = -1)
    plist[[gene]] <- ggrastr::rasterise(p, layers="Point", dpi=300)
  }
  
  
  # combine featureplots
  pcombined <- plot_grid(
    plot_grid(plotlist=plist, ncol=2, scale=1),
    ggpubr::as_ggplot(plegend),
    ncol=2, rel_widths = c(5, .5), scale=c(0.95, 0.5)
    )
  
  return(pcombined)
}

plot_genesignature <- function(seur, genesignature, ordercells=F){
  # get legend
  plegend <- ggpubr::get_legend(
    SCpubr::do_FeaturePlot(seur, features=genesignature, legend.position="right", border.size=1, pt.size=2, order=ordercells)+
      scale_color_viridis_c(limits=c(0,max(seur@meta.data[,genesignature])), option="B")
  )
  # get featureplot
  p <- SCpubr::do_FeaturePlot(seur, features=genesignature, legend.position="none", border.size=1, pt.size=2, order=ordercells)+
    scale_color_viridis_c(limits=c(0,max(seur@meta.data[,genesignature])), option="B")+
    theme(plot.background = element_rect(color = "black", linewidth = 2))
  p <- ggrastr::rasterise(p, layers="Point", dpi=300)
  # combine
  pcombined <- plot_grid(
    p,
    ggpubr::as_ggplot(plegend) + theme(plot.margin=margin(0,50,0,30)),
    ncol=2, rel_widths = c(5, 2), scale=c(0.95, 0.5)
  )
  return(pcombined)
}




# ****************
# 3. ANALYSIS ####
# ****************

#___________________________
## 3.1. Score gene signatures ####
gene_signatures <- list("effector_new"=c("HOPX", "GZMB", "NKG7", "TBX21", "PRF1", "GZMA", "KLRD1", "KLF6",
                                     "CCR6", "RORC", "JUNB", "FOS", "RORA", "FOSB"),
                        "effector"= c("HOPX", "GZMB", "GZMK", "ZEB2", "NKG7", "GNLY", "TBX21", "EOMES", "TYROBP", "PRF1",
                                      "CCL4", "CCL5", "KLRB1", "GZMH", "GZMA", "KLRD1", "CST7", "KLF6", "CXCR4"),
                        "naive"=c("SATB1", "TCF7", "LEF1", "CCR7", "SELL", "MYC", "EIF3E", "SOX4", "ID3", "BACH2"),
                        "egress"=c("KLF2", "CORO1A", "CCR7", "CXCR4", "CXCR6", "FOXO1", "CXCR3", "S1PR1", "S1PR4",
                                   "S100A4", "S100A6", "EMP3"))

seur.nkt   <- AddModuleScore(seur.nkt,  name = names(gene_signatures), features=gene_signatures, seed=1)
seur.mait  <- AddModuleScore(seur.mait, name = names(gene_signatures), features=gene_signatures, seed=1)
seur.gdt   <- AddModuleScore(seur.gdt,  name = names(gene_signatures), features=gene_signatures, seed=1)
seur.cd4   <- AddModuleScore(seur.cd4,  name = names(gene_signatures), features=gene_signatures, seed=1)
seur.cd8   <- AddModuleScore(seur.cd8,  name = names(gene_signatures), features=gene_signatures, seed=1)

colnames(seur.nkt@meta.data)[35:38]  <- names(gene_signatures)
colnames(seur.mait@meta.data)[35:38] <- names(gene_signatures)
colnames(seur.gdt@meta.data)[35:38]  <- names(gene_signatures)
colnames(seur.cd4@meta.data)[35:38]  <- names(gene_signatures)
colnames(seur.cd8@meta.data)[35:38]  <- names(gene_signatures)

# Score same signatures in mouse nkt/mait
gene_signatures_ms <- list("effector_new"=c("Hopx", "Gzmb", "Nkg7", "Tbx21", "Prf1", "Gzma", "Klrd1", "Klf6",
                                        "Ccr6", "Rorc", "Tmem176a", "Tmem176b", "Junb", "Fos", "Rora", "Fosb"),
                           "effector"=c("Hopx", "Gzmb", "Gzmk", "Zeb2", "Nkg7", "Tbx21", "Eomes", "Tyrobp", "Prf1",
                                        "Ccl4", "Ccl5", "Klrb1a", "Gzmg", "Gzma", "Klrd1", "Cst7", "Klf6", "Cxcr4"), # manual curating
                           "naive"=ortholog.df[ortholog.df$hu_symbol %in% gene_signatures$naive, "ms_symbol_data"],
                           "egress"=ortholog.df[ortholog.df$hu_symbol %in% gene_signatures$egress, "ms_symbol_data"])

seur.nkt.ms   <- AddModuleScore(seur.nkt.ms,  name = names(gene_signatures_ms), features=gene_signatures_ms, seed=1)
seur.mait.ms  <- AddModuleScore(seur.mait.ms, name = names(gene_signatures_ms), features=gene_signatures_ms, seed=1)

colnames(seur.nkt.ms@meta.data)[8:10]  <- names(gene_signatures_ms)
colnames(seur.mait.ms@meta.data)[10:12] <- names(gene_signatures_ms)
## /end ####


#___________________________
## 3.2. Plot gene signatures ####

# human
plot_genesignature(seur.nkt, genesignature = "naive", ordercells=F)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_naivesig.jpeg", width=7, height=5)
plot_genesignature(seur.mait, genesignature = "naive", ordercells=F)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_naivesig.jpeg", width=7, height=5)

plot_genesignature(seur.nkt, genesignature = "effector_new", ordercells=F)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_effectorsigSC.jpeg", width=7, height=5)
plot_genesignature(seur.mait, genesignature = "effector_new", ordercells=F)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_effectorsigSC.jpeg", width=7, height=5)
plot_genesignature(seur.gdt, genesignature = "effector_new", ordercells=F)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_gdt_effectorsigSC.jpeg", width=7, height=5)



# mouse
plot_genesignature(seur.nkt.ms, genesignature = "naive")
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_naivesig_mouse.jpeg", width=7, height=5)
plot_genesignature(seur.mait.ms, genesignature = "naive")
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_naivesig_mouse.jpeg", width=7, height=5)

plot_genesignature(seur.nkt.ms, genesignature = "effector")
plot_genesignature(seur.mait.ms, genesignature = "effector")


## /end ####



#___________________________
## 3.3. Plot genes of interest ####

# naive genes
plot_four_genes(seur.nkt, c("CCR9", "CCR7", "SELL", "TCF7"))
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_naivegenes.jpeg", width=6, height=5)
plot_four_genes(seur.mait, c("CCR9", "CCR7", "SELL", "TCF7"))
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_naivegenes.jpeg", width=6, height=5)


# effector genes
plot_four_genes(seur.nkt, c("ZBTB16", "EOMES", "KLRB1", "GZMK"))
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_effectorgenes.jpeg", width=6, height=5)
plot_four_genes(seur.mait, c("ZBTB16", "EOMES", "KLRB1", "GZMK"))
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_effectorgenes.jpeg", width=6, height=5)
plot_four_genes(seur.gdt, c("KLRD1", "EOMES", "KLRB1", "GZMK"), ordercells=T)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_gdt_effectorgenes.jpeg", width=6, height=5)


# CD4 and CD8
FeaturePlot(seur.nkt, features=c("CD4", "CD8A"), blend=T, blend.threshold=0.01, cols=c("lightgrey", "red", "blue"), order=T, pt.size=2)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_cd4cd8.jpeg", width=18, height=5)
FeaturePlot(seur.mait, features=c("CD4", "CD8A"), blend=T, cols=c("lightgrey", "red", "blue"), order=T, pt.size=2)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_cd4cd8.jpeg", width=18, height=5)
FeaturePlot(seur.gdt, features=c("CD4", "CD8A"), blend=T, cols=c("lightgrey", "red", "blue"), order=T, pt.size=2)
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_gdt_cd4cd8.jpeg", width=18, height=5)


plot_four_genes(seur.nkt, c("CD4", "CD8A"))
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_effectorgenes.jpeg", width=6, height=3)
plot_four_genes(seur.mait, c("ZBTB16", "EOMES", "GZMK", "KLRB1"))
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_effectorgenes.jpeg", width=6, height=5)
plot_four_genes(seur.gdt, c("ZBTB16", "EOMES", "GZMK", "KLRB1"))
ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_gdt_effectorgenes.jpeg", width=6, height=5)


# CCR9 and CCR7
FeaturePlot(seur.nkt, features=c("CCR7", "CCR9"), blend=T, blend.threshold=0.01, cols=c("lightgrey", "red", "blue"), order=T, pt.size=2)

## /end ####


#___________________________
## 3.4. Plot lineage composition in clusters 12-17 ####

# Import data
seur.human <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_08_28_23.RDS")
DimPlot(seur.human, group.by="new_clusters", reduction="UMAP_50")

# keep only thymocytes
seur.human <- subset(seur.human, subset=Tissue=="Thymus")
df <- seur.human@meta.data[,c("new_clusters", "cell.ident")]

df %>%
  mutate(new_clusters=paste0("c", new_clusters)) %>%
  group_by(new_clusters, cell.ident) %>%
  summarise(nb_cells=n()) %>%
  ungroup() %>%
  group_by(new_clusters) %>%
  mutate(total_cells_per_cluster=sum(nb_cells),
         freq_lineage_per_cluster=nb_cells*100/total_cells_per_cluster) %>%
  filter(new_clusters %in% paste0("c", 12:17))%>%
  # plot
  ggplot(aes(x=new_clusters, y=nb_cells, fill=cell.ident)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=cols_lineages, name="")+
  labs(x="", y="#cells")+
  theme_cowplot()+
  theme(axis.text = element_text(size=20),
        axis.title=element_text(size=20))

## /end ####


#___________________________
## 3.5. Run Metaneighbor between T lineages ####
library(MetaNeighbor)
library(SummarizedExperiment)
library(gplots)
library(RColorBrewer)

# rename CD8 clusters
# table(seur.cd8@meta.data$cell_annot, useNA="ifany")
seur.cd8@meta.data[,"cell_annot"] <- case_when(
  seur.cd8@meta.data$cell_annot=="thyCD8_DP"     ~ "CD8_c0",
  seur.cd8@meta.data$cell_annot=="thyCD8_cd8aa1" ~ "CD8_c1",
  seur.cd8@meta.data$cell_annot=="thyCD8_cd8aa2" ~ "CD8_c2",
  seur.cd8@meta.data$cell_annot=="thyCD8_ccr9"   ~ "CD8_c3",
  seur.cd8@meta.data$cell_annot=="thyCD8_ccr7"   ~ "CD8_c4",
  seur.cd8@meta.data$cell_annot=="thyCD8_idk"    ~ "CD8_c5",
  .default="NA"
)
# table(seur.cd8@meta.data$cell_annot, useNA="ifany")

# rename CD4 clusters
# table(seur.cd4@meta.data$cell_annot, useNA="ifany")
seur.cd4@meta.data[,"cell_annot"] <- case_when(
  seur.cd4@meta.data$cell_annot=="thyCD4_ISP"      ~ "CD4_c0",
  seur.cd4@meta.data$cell_annot=="thyCD4_DPp"      ~ "CD4_c1",
  seur.cd4@meta.data$cell_annot=="thyCD4_DPq"      ~ "CD4_c2",
  seur.cd4@meta.data$cell_annot=="thyCD4_ccr9"     ~ "CD4_c3",
  seur.cd4@meta.data$cell_annot=="thyCD4_ccr7"     ~ "CD4_c4",
  seur.cd4@meta.data$cell_annot=="thyCD4_Tagonist" ~ "CD4_c5",
  seur.cd4@meta.data$cell_annot=="thyCD4_Treg"     ~ "CD4_c6",
  .default="NA"
)
# table(seur.cd4@meta.data$cell_annot, useNA="ifany")

# get union of HVGs
hvg.all <- unique(c(VariableFeatures(seur.nkt),
                    VariableFeatures(seur.mait),
                    VariableFeatures(seur.gdt),
                    VariableFeatures(seur.cd4),
                    VariableFeatures(seur.cd8)))
length(hvg.all) # 5864 genes

# merge seurat objects
seur.thym <- merge(seur.cd4, y = c(seur.cd8, seur.mait, seur.nkt, seur.gdt))
table(seur.thym@meta.data$cell_annot, useNA="ifany")

# make into summarized experiment for metaneighbor
se <- SummarizedExperiment(assays=seur.thym@assays[["RNA"]]@counts,
                           colData=seur.thym@meta.data[,c("cell.ident", "cell_annot")])

# run metaneighbor
mtn <- MetaNeighborUS(var_genes=hvg.all,
                      dat=se,
                      study_id=seur.thym$cell.ident,
                      cell_type=seur.thym$cell_annot,
                      fast_version=FALSE)

# plot full dendrogram
heatmap.2(mtn,
          # trace
          trace="none",
          # superimpose a density histogram on color key
          density.info="none",
          # color scale
          col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100)),
          breaks=seq(0,1,length=101),
          key.xlab="AUROC",
          # text labels
          cexRow=0.6,
          cexCol=0.6,
          # colRow=col_text,
          # colCol=col_text,
          # margins
          margins=c(7,7))

# Bubble plot
library(ggrepel)
mtn.df <- reshape2::melt(mtn)
mtn.df1 <- mtn.df %>%
  mutate(Var1 = gsub(".*\\|", "", Var1)) %>%
  mutate(Var2 = gsub(".*\\|", "", Var2)) %>%
  as_tibble() %>%
  dplyr::rename(auroc=value) %>%
  filter(Var1 %in% grep("CD4_c|CD8_c", Var1, value=T),
         Var2 %in% grep("NKT_c|MAIT_c", Var2, value=T))

ggplot(mtn.df1, aes(x=factor(Var1, levels=c("CD8_c0", "CD4_c0", "CD4_c1", "CD4_c2", "CD8_c1", # DP/cd8aa
                                           "CD4_c3", "CD8_c3", # ccr9
                                           "CD4_c4", "CD8_c4", # ccr7
                                           "CD8_c5", "CD4_c5", "CD4_c6", 
                                           "CD8_c2")), # effector
             y=factor(Var2, levels=rev(c("NKT_c0", "MAIT_c0", "MAIT_c1", # DP/cd8aa
                                         "NKT_c1", "MAIT_c2", # ccr9
                                         "NKT_c2", "MAIT_c4", "MAIT_c5", # ccr7
                                         "NKT_c3", "NKT_c4", "MAIT_c3", # Treg/agonist
                                         "NKT_c5", "NKT_c6", "MAIT_c6"))))) +
  geom_point(aes(size = auroc, color= auroc))+
  geom_text(data=mtn.df1 %>% filter(auroc>0.8) %>% mutate(across("auroc", \(x) round(x,2))), aes(label=auroc), color="white")+
  scale_size_continuous(limits=c(0,1), breaks=seq(0,1, by=0.2), range = c(1, 15))+
  scale_color_gradient2(low="#2166ac", mid="white", high="#a50f15", midpoint=0.5, limits=c(0,1), name="AUROC", breaks=seq(0,1, by=0.2))+
  labs(x="innate T",y="conventional T", size="AUROC")+
  theme_cowplot()+
  theme(legend.position="bottom", legend.key.width = unit(0.8, 'cm'),
        axis.text = element_text(size=15), axis.text.x = element_text(angle=45, hjust=1), axis.title=element_text(size=20))
ggsave("./data/human-thymus/HumanData_13_CorrelationGEbtwTlineages/metaneighbor/20231031_metaneighbor_unionHVGs.pdf", width=10, height=10)


clusters_naive <- c("CD4_c3", "CD8_c3","NKT_c1", "MAIT_c2","CD4_c4", "CD8_c4","NKT_c2", "MAIT_c4", "MAIT_c5")
mtn.df %>%
  mutate(Var1 = gsub(".*\\|", "", Var1)) %>%
  mutate(Var2 = gsub(".*\\|", "", Var2)) %>%
  as_tibble() %>%
  dplyr::rename(auroc=value) %>%
  mutate(across("auroc", \(x) round(x,2))) %>%
  filter(Var1 %in% clusters_naive,
         Var2 %in% clusters_naive) %>%
  # PLOT
  ggplot(aes(x=factor(Var1, levels=clusters_naive), # effector
             y=factor(Var2, levels=rev(clusters_naive)))) +
  geom_point(aes(size = auroc, color= auroc))+
  # geom_text(aes(label=auroc), color="white")+
  scale_size_continuous(limits=c(0,1), breaks=seq(0,1, by=0.2), range = c(1, 15))+
  scale_color_gradient2(low="#2166ac", mid="white", high="#a50f15", midpoint=0.5, limits=c(0,1), name="AUROC", breaks=seq(0,1, by=0.2))+
  labs(x="",y="", size="AUROC")+
  theme_cowplot()+
  theme(legend.position="bottom", legend.key.width = unit(0.8, 'cm'),
        axis.text = element_text(size=15), axis.text.x = element_text(angle=45, hjust=1), axis.title=element_text(size=20))

## /end ####


# *****************
# 4. TCR USAGE ####
# *****************

#___________________________
## 4.1. GD cells ####

gd.tcr <- seur.gdt@meta.data %>%
  rownames_to_column("cellid") %>%
  as_tibble() %>%
  select(cellid, cell_annot, grep("_gene_Dominant", colnames(seur.gdt@meta.data), value=T)) %>%
  dplyr::rename(TRGV = TCR_Alpha_Gamma_V_gene_Dominant,
                TRGJ = TCR_Alpha_Gamma_J_gene_Dominant,
                TRDV = TCR_Beta_Delta_V_gene_Dominant,
                TRDD = TCR_Beta_Delta_D_gene_Dominant,
                TRDJ = TCR_Beta_Delta_J_gene_Dominant) %>%
  na.omit(!cell_annot) %>%
  mutate(TRGV = str_remove(TRGV, pattern = "\\*[^.]*$"),
         TRGJ = str_remove(TRGJ, pattern = "\\*[^.]*$"),
         TRDV = str_remove(TRDV, pattern = "\\*[^.]*$"),
         TRDD = str_remove(TRDD, pattern = "\\*[^.]*$"),
         TRDJ = str_remove(TRDJ, pattern = "\\*[^.]*$"))

# highlight TRDV1, TRDV2, TRDV3 and TRGV9
plot_grid(
  SCpubr::do_DimPlot(sample = seur.gdt, cells.highlight = pull(gd.tcr[gd.tcr$TRDV=="TRDV1","cellid"]), plot.title = "TRDV1",
                     legend.position = "none", na.value = "grey90", colors.use = "#a40000"),
  SCpubr::do_DimPlot(sample = seur.gdt, cells.highlight = pull(gd.tcr[gd.tcr$TRDV=="TRDV2","cellid"]), plot.title = "TRDV2",
                     legend.position = "none", na.value = "grey90", colors.use = "blue"),
  SCpubr::do_DimPlot(sample = seur.gdt, cells.highlight = pull(gd.tcr[gd.tcr$TRDV=="TRDV3","cellid"]), plot.title = "TRDV3",
                     legend.position = "none", na.value = "grey90", colors.use = "#318f49"),
  SCpubr::do_DimPlot(sample = seur.gdt, cells.highlight = pull(gd.tcr[gd.tcr$TRGV=="TRGV9","cellid"]), plot.title = "TRGV9",
                     legend.position = "none", na.value = "grey90", colors.use = "gold"),
  nrow=2)
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_gdt_highlight_tcrusage.jpeg", width=5, height=5)


# barplot TCR usage per cluster
cols_TRDJ <- c("TRDJ1" = "#16317d", "TRDJ2" = "#3BB38C", "TRDJ3" = "#007e2f", "TRDJ4" = "#ffcd12")
gd.tcr %>%
  mutate(cluster=ifelse(cell_annot=="GDT_c7", "GD_c7", "GD_c0-6")) %>%
  filter(TRDV == "TRDV2") %>%
  summarise(n=n(), .by=c(cluster, TRDV, TRDJ)) %>%
  ggplot(aes(x=cluster, y=n, fill=TRDJ))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=cols_TRDJ)+
  labs(x='',y="#cells", title="TRDV2")+
  theme_cowplot()

cols_TRG <- c("TRGV1" = "#CF597E", "TRGV2" = "#a40000", "TRGV3" = "#16317d", "TRGV4" = "#007e2f", "TRGV5" = "#ffcd12",
              "TRGV5P" = "#D86279", "TRGV8" = "#1FA990", "TRGV9" = "#E9A26A", "TRGV10" = "#DB6577", "TRGV11" = "#E16C72", 
              "TRGJ1" = "#93CB83", 
              "TRGJP" = "#C7D88D", "TRGJP1" = "#EADB94", "TRGJP2" = "#EACE85")
gd.tcr %>%
  mutate(cluster=ifelse(cell_annot=="GDT_c7", "GD_c7", "GD_c0-6")) %>%
  filter(TRDV == "TRDV2") %>%
  summarise(n=n(), .by=c(cluster, TRGV)) %>%
  ggplot(aes(x=cluster, y=n, fill=TRGV))+
  geom_bar(stat="identity", position="fill")+
  scale_fill_manual(values=cols_TRG, name="")+
  labs(x='',y="%cells", title="γV")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=25),
        axis.title.y=element_text(size=25),
        title=element_text(size=35))
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_gdt_tcrusage_gammachain.jpeg", width=5, height=7)

## /end ####


#___________________________
## 4.2. iNKT cells ####

nkt.tcr <- seur.nkt@meta.data %>%
  rownames_to_column("cellid") %>%
  as_tibble() %>%
  select(cellid, cell_annot, grep("_gene_Dominant", colnames(seur.gdt@meta.data), value=T)) %>%
  dplyr::rename(TRAV = TCR_Alpha_Gamma_V_gene_Dominant,
                TRAJ = TCR_Alpha_Gamma_J_gene_Dominant,
                TRBV = TCR_Beta_Delta_V_gene_Dominant,
                TRBD = TCR_Beta_Delta_D_gene_Dominant,
                TRBJ = TCR_Beta_Delta_J_gene_Dominant) %>%
  na.omit(!cell_annot) %>%
  mutate(TRAV = str_remove(TRAV, pattern = "\\*[^.]*$"),
         TRAJ = str_remove(TRAJ, pattern = "\\*[^.]*$"),
         TRBV = str_remove(TRBV, pattern = "\\*[^.]*$"),
         TRBD = str_remove(TRBD, pattern = "\\*[^.]*$"),
         TRBD = case_when(TRBD=="" ~ "NA", .default=TRBD),
         TRBJ = str_remove(TRBJ, pattern = "\\*[^.]*$")) %>%
  filter(TRAV=="TRAV10")

# highlight TRAV10
SCpubr::do_DimPlot(sample = seur.nkt, cells.highlight = pull(nkt.tcr[,"cellid"]), plot.title = "TRAV10",
                   legend.position = "none", na.value = "grey90", colors.use = "#a40000")
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_highlightTRAV.jpeg", width=3, height=3)


# barplot TCR usage per cluster
cols.TRBV <- c("TRBV4-2" = "#EAF1F4", "TRBV5-1" = "#93D4DB", "TRBV6-5" = "#99CCA7", "TRBV7-3" = "#57AA65", "TRBV12-2" = "#A7A365", 
              "TRBV19" = "#C2947E", "TRBV21-1" = "#785158", "TRBV23-1" = "#757575", "TRBV25-1" = "#3D7D8F", "TRBV27" = "#956E73")
cols.TRBD <- c("NA"="grey", "TRBD1"="#9DCD84", "TRBD2"="#E9AA6C")
cols.TRBJ <- c("TRBJ1-1" = "#EB924A", "TRBJ1-2" = "#3EACCC", "TRBJ1-3" = "#61C2DA", "TRBJ1-4" = "#F4CBA0", "TRBJ1-5" = "#A7A7C7", "TRBJ1-6" = "#E8A1CF", 
               "TRBJ2-1" = "#E7D9DB", "TRBJ2-2" = "#736B9D", "TRBJ2-3" = "#DEE9E9", "TRBJ2-4" = "#EF4F55", "TRBJ2-5" = "#76C76C", "TRBJ2-6" = "#1E2223", "TRBJ2-7" = "#E7D9DB")

plot_grid(
  nkt.tcr %>%
    summarise(n=n(), .by=c(cell_annot, TRBV)) %>%
    group_by(cell_annot) %>% filter(sum(n)>10) %>% ungroup() %>%
    ggplot(aes(x=cell_annot, y=n, fill=TRBV))+
    geom_bar(stat="identity", position="fill")+
    scale_fill_manual(values=cols.TRBV, name="")+
    labs(x='',y="%cells", title="βV")+
    theme_cowplot()+
    theme(axis.text.x=element_text(angle=45, hjust=1, size=20),
          axis.title.y=element_text(size=20),
          title=element_text(size=30)),
  nkt.tcr %>%
    summarise(n=n(), .by=c(cell_annot, TRBD)) %>%
    group_by(cell_annot) %>% filter(sum(n)>10) %>% ungroup() %>%
    ggplot(aes(x=cell_annot, y=n, fill=TRBD))+
    geom_bar(stat="identity", position="fill")+
    scale_fill_manual(values=cols.TRBD, name="")+
    labs(x='',y="%cells", title="βD")+
    theme_cowplot()+
    theme(axis.text.x=element_text(angle=45, hjust=1, size=20),
          axis.title.y=element_text(size=20),
          title=element_text(size=30)),
  nkt.tcr %>%
    summarise(n=n(), .by=c(cell_annot, TRBJ)) %>%
    group_by(cell_annot) %>% filter(sum(n)>10) %>% ungroup() %>%
    ggplot(aes(x=cell_annot, y=n, fill=TRBJ))+
    geom_bar(stat="identity", position="fill")+
    scale_fill_manual(values=cols.TRBJ, name="")+
    labs(x='',y="%cells", title="βJ")+
    theme_cowplot()+
    theme(axis.text.x=element_text(angle=45, hjust=1, size=20),
          axis.title.y=element_text(size=20),
          title=element_text(size=30)),
  ncol=3)
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_nkt_tcrusage_betachain.jpeg", width=15, height=5)
## /end ####


#___________________________
## 4.3. MAIT cells ####

mait.tcr <- seur.mait@meta.data %>%
  rownames_to_column("cellid") %>%
  as_tibble() %>%
  select(cellid, cell_annot, grep("_gene_Dominant", colnames(seur.gdt@meta.data), value=T)) %>%
  dplyr::rename(TRAV = TCR_Alpha_Gamma_V_gene_Dominant,
                TRAJ = TCR_Alpha_Gamma_J_gene_Dominant,
                TRBV = TCR_Beta_Delta_V_gene_Dominant,
                TRBD = TCR_Beta_Delta_D_gene_Dominant,
                TRBJ = TCR_Beta_Delta_J_gene_Dominant) %>%
  na.omit(!cell_annot) %>%
  mutate(TRAV = str_remove(TRAV, pattern = "\\*[^.]*$"),
         TRAJ = str_remove(TRAJ, pattern = "\\*[^.]*$"),
         TRBV = str_remove(TRBV, pattern = "\\*[^.]*$"),
         TRBD = str_remove(TRBD, pattern = "\\*[^.]*$"),
         TRBD = case_when(TRBD=="" ~ "NA", .default=TRBD),
         TRBJ = str_remove(TRBJ, pattern = "\\*[^.]*$")) %>%
  filter(TRAV=="TRAV1-2" & TRBV !="TRDV2" & TRBJ != "TRDJ1")

# highlight TRAV1-2
SCpubr::do_DimPlot(sample = seur.mait, cells.highlight = pull(mait.tcr[,"cellid"]), plot.title = "TRAV1-2",
                   legend.position = "none",na.value = "grey90", colors.use = "#a40000")
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_highlightTRAV.jpeg", width=3, height=3)

# barplot TCR usage per cluster
trbv_mait <- unique(mait.tcr$TRBV)
cols.TRBV.extend <- hcl.colors(17, "Temps")
names(cols.TRBV.extend) <- trbv_mait[!trbv_mait %in% names(cols.TRBV)]
cols.TRBV <- c(cols.TRBV, cols.TRBV.extend)
# length(cols.TRBV)
# table(trbv_mait %in% names(cols.TRBV))

plot_grid(
  mait.tcr %>%
    summarise(n=n(), .by=c(cell_annot, TRBV)) %>%
    group_by(cell_annot) %>% filter(sum(n)>10) %>% ungroup() %>%
    ggplot(aes(x=cell_annot, y=n, fill=TRBV))+
    geom_bar(stat="identity", position="fill")+
    scale_fill_manual(values=cols.TRBV, name="")+
    labs(x='',y="%cells", title="βV")+
    theme_cowplot()+
    theme(axis.text.x=element_text(angle=45, hjust=1, size=25),
          axis.title.y=element_text(size=25),
          title=element_text(size=35)),
  mait.tcr %>%
    summarise(n=n(), .by=c(cell_annot, TRBD)) %>%
    group_by(cell_annot) %>% filter(sum(n)>10) %>% ungroup() %>%
    ggplot(aes(x=cell_annot, y=n, fill=TRBD))+
    geom_bar(stat="identity", position="fill")+
    scale_fill_manual(values=cols.TRBD, name="")+
    labs(x='',y="%cells", title="βD")+
    theme_cowplot()+
    theme(axis.text.x=element_text(angle=45, hjust=1, size=25),
          axis.title.y=element_text(size=25),
          title=element_text(size=35)),
  mait.tcr %>%
    summarise(n=n(), .by=c(cell_annot, TRBJ)) %>%
    group_by(cell_annot) %>% filter(sum(n)>10) %>% ungroup() %>%
    ggplot(aes(x=cell_annot, y=n, fill=TRBJ))+
    geom_bar(stat="identity", position="fill")+
    scale_fill_manual(values=cols.TRBJ, name="")+
    labs(x='',y="%cells", title="βJ")+
    theme_cowplot()+
    theme(axis.text.x=element_text(angle=45, hjust=1, size=25),
          axis.title.y=element_text(size=25),
          title=element_text(size=35)),
  ncol=3, rel_widths = c(1.2,1,1))
# ggsave("./scripts-in-progress/human-thymus/HumanThymus_24_PlotGenesofInterest/plots/thymus_mait_tcrusage_betachain.jpeg", width=18, height=6)

## /end ####