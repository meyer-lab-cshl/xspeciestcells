# Purpose: Compare gene lists (GEPs, Cano-Gamez, Rose, Poon signatures)
# Author: Salomé Carcy
# Date: June 2023


# **************
# 1. IMPORT ####
# **************

# Import librairies
library(Seurat)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(pheatmap)

# Import data
genes.gep     <- read.csv("./data/human-thymus/HumanData_17_GEPsOnParkData/genes_per_GEP_df_2023-04-07.csv", row.names=1)
genes.cano    <- readxl::read_excel("./data/human-thymus/HumanData_17_GEPsOnCanogamezData/canogamez_supp.xlsx", skip=3, sheet=1)
genesCD8.rose <- readxl::read_excel("./data/human-thymus/HumanData_22_CompareGeneLists/rose_supp_clusmodulesCD8.xlsx", sheet=1)
genesCD4.rose <- readxl::read_excel("./data/human-thymus/HumanData_22_CompareGeneLists/rose_supp_clusmodulesCD4.xlsx", sheet=1)
genes.poon    <- readxl::read_excel("./data/human-thymus/HumanData_22_CompareGeneLists/poon_supp.xlsx", sheet=6)[,-1]


# Sanity check: reproduce Rose heatmaps
htmp <- function(rosedf, lineage){
  # Get rpkm
  if(lineage=="CD4"){ counts <- rosedf[,9:20] }
  else if(lineage=="CD8"){ counts <- rosedf[,9:24] }
  rownames(counts) <- rosedf$ENTREZID
  # Get metadata
  meta <- rosedf[,c(1,3)] %>%
    column_to_rownames("ENTREZID")
  cols <- RColorBrewer::brewer.pal(5, "Greys")
  names(cols) <- unique(meta$k.5.cluster)
  # Plot
  pheatmap(counts,
           scale="row",
           cluster_rows=F,
           cluster_cols=F,
           color=rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(100)),
           annotation_row=meta,
           annotation_colors=list(k.5.cluster=cols),
           show_rownames=F)
}

htmp(genesCD4.rose, lineage="CD4")
htmp(genesCD8.rose, lineage="CD8")




# *****************************************
# 2. PREPARE GENE LISTS IN LONG FORMAT ####
# *****************************************

genes.long <- list()

# ***************
## 2.1. GEPs ####
head(genes.gep)
gapindf <- gather(genes.gep, key=geneprogram, value=gene, colnames(genes.gep)) %>%
  filter(!is.na(gene)) %>%
  filter(geneprogram %in% c("GEP_1", "GEP_4", "GEP_5", "GEP_6", "GEP_7")) %>%
  mutate(geneprogram=sub("_", "", geneprogram),
         geneprogram=replace(geneprogram, geneprogram=="GEP1", "GEP1_MAIT"),
         geneprogram=replace(geneprogram, geneprogram=="GEP4", "GEP4_Temra"),
         geneprogram=replace(geneprogram, geneprogram=="GEP5", "GEP5_Tnaive"),
         geneprogram=replace(geneprogram, geneprogram=="GEP6", "GEP6_Tcm"))
head(gapindf)
table(gapindf$geneprogram, useNA="ifany")
table(is.na(gapindf$gene))

genes.long[["gapin"]] <- gapindf
# table(is.na(genes.long$gapin$gene)) # sanity check


# **************************
## 2.2. Cano-Gamez data ####
head(genes.cano)
# take a look at the different clusters/gene programs available
table(genes.cano$cluster, useNA="ifany")
# rename the clusters a bit
canodf <- genes.cano %>%
  mutate(cluster=gsub(" ", "", cluster),
         cluster=paste0("CanoGamez | CD4_", cluster)) %>%
  # keep only significant genes & with minimum log2FC
  filter(p_val_adj < 0.05 & LFC > 0) %>%
  # convert to long format
  select(cluster, gene) %>%
  rename(geneprogram=cluster)
head(canodf)
genes.long[["canogamez"]] <- canodf
# table(is.na(genes.long$canogamez$gene)) # sanity check

# Define colors
cols_cano <- c("CanoGamez | CD4_Tnaive"= "#b3e2cd",
               "CanoGamez | CD4_TCM"   = "#f4cae4",
               "CanoGamez | CD4_TEM"   = "#cbd5e8",
               "CanoGamez | CD4_TEMRA" = "#fdcdac",
               "CanoGamez | CD4_nTreg" = "#fbb4ae")


# **************************
## 2.3. Rose data ####
head(genesCD4.rose)
head(genesCD8.rose)

roseCD4 <- genesCD4.rose %>%
  select(SYMBOL, k.5.cluster) %>%
  filter(!is.na(SYMBOL)) %>% # the few NA symbols are miRNA or uncharacterized loci
  rename(gene=SYMBOL, geneprogram=k.5.cluster) %>%
  mutate(geneprogram=ifelse(geneprogram==1, "Rose | CD4_modul1_Tem",
                     ifelse(geneprogram==2, "Rose | CD4_modul2_Tcm/em",
                     ifelse(geneprogram==3, "Rose | CD4_modul3_Tem",
                     ifelse(geneprogram==4, "Rose | CD4_modul4_Tem",
                     ifelse(geneprogram==5, "Rose | CD4_modul5_Tnaive", "?"))))))

roseCD8 <- genesCD8.rose %>%
  select(SYMBOL, k.5.cluster) %>%
  filter(!is.na(SYMBOL)) %>% # the few NA symbols are miRNA or uncharacterized loci
  rename(gene=SYMBOL, geneprogram=k.5.cluster) %>%
  mutate(geneprogram=ifelse(geneprogram==1, "Rose | CD8_modul1_Temra",
                     ifelse(geneprogram==2, "Rose | CD8_modul2_Tcm/em",
                     ifelse(geneprogram==3, "Rose | CD8_modul3_Tem/emra",
                     ifelse(geneprogram==4, "Rose | CD8_modul4_Tnaive",
                     ifelse(geneprogram==5, "Rose | CD8_modul5_Tnaive/cm", "?"))))))

head(rbind(roseCD4, roseCD8))
genes.long[["rose"]] <- rbind(roseCD4, roseCD8)
# table(is.na(genes.long$rose$gene)) # sanity check
# table(genes.long$rose$geneprogram, useNA="ifany")

# Define colors
cols_rose <- c("Rose | CD4_modul1_Tem"      ="#cbd5e8",
               "Rose | CD4_modul2_Tcm/em"   ="#f4cae4",
               "Rose | CD4_modul3_Tem"      ="#cbd5e8",
               "Rose | CD4_modul4_Tem"      ="#cbd5e8",
               "Rose | CD4_modul5_Tnaive"   ="#b3e2cd",
               "Rose | CD8_modul1_Temra"    ="#fdcdac",
               "Rose | CD8_modul2_Tcm/em"   ="#f4cae4",
               "Rose | CD8_modul3_Tem/emra" ="#cbd5e8",
               "Rose | CD8_modul4_Tnaive"   ="#b3e2cd",
               "Rose | CD8_modul5_Tnaive/cm"="#b3e2cd")


# ********************
## 2.4. Poon data ####
head(genes.poon)
colnames(genes.poon) <- gsub(" ", "", colnames(genes.poon))

# Transform to long format gene program by gene program (each gene program has 3 columns: genes symbols, padj, logFC)
poondf <- data.frame()
for (i in seq(1,ncol(genes.poon), 3)){
  subdf <- genes.poon[,i:(i+2)]
  geneprogram <- sub("_.*", "",colnames(subdf)[1])
  colnames(subdf) <- c("gene", "padj", "logFC")
  subdf <- subdf %>%
    filter(padj<0.01 & logFC>0) %>%
    slice_head(n=200) # take 200 top genes
    # slice_max(logFC, n=200)
  subdf$geneprogram <- paste0("Poon | ", geneprogram) # add column with the name of the geneprogram
  poondf <- rbind(poondf, subdf)
}
head(poondf)
table(poondf$geneprogram, useNA="ifany")
table(is.na(poondf$gene))

genes.long[["poon"]] <- poondf[,c("geneprogram", "gene")]

# Define colors
cols_poon <- c("Poon | CyclingTRM"="#f1e2cc",
               "Poon | CD4Naive"="#b3e2cd",
               "Poon | CD4TCM/TFH"="#f4cae4",
               "Poon | CD4TRM"="#fff2ae",
               "Poon | CD4Treg"="#fbb4ae",
               "Poon | CD8MAIT"="#cbd5e8",
               "Poon | CD8Naive"="#b3e2cd",
               "Poon | CD8TEM/TEMRA"="#fdcdac",
               "Poon | CD8TRM"="#fff2ae")


# *******************************
## 2.5. Get it all in one df ####
longdf <- bind_rows(genes.long, .id="dataset")
table(longdf$dataset, useNA="ifany")

cols_alldatasets <- c(cols_cano, cols_rose, cols_poon)



# ******************
# 3. UPSET PLOT ####
# ******************
library(ggupset)

# *********************
## 3.1. Upset plot ####
# Gapin VS 1 other
# jpeg("./data/human-thymus/HumanData_22_CompareGeneLists/gapin_poon_slicehead.jpeg", width=3000, height=1200, res=200)
longdf %>%
  filter(dataset %in% c("gapin", "poon")) %>% # restrict to only 2 datasets
  group_by(gene) %>%
  # keep only genes that are present in >1 dataset (not interested in intersections intra-dataset)
  filter(n_distinct(dataset)>1) %>% 
  summarize(geneprograms = list(geneprogram)) %>%
  # filter(lengths(geneprograms)==2) %>%
  ggplot(aes(x = geneprograms)) +
    geom_bar() +
    scale_x_upset() +
    labs(title="Gapin vs Poon")
# dev.off()


# All datasets together
gapingenes <- unique(genes.long$gapin$gene)
longdf %>%
  filter(gene %in% gapingenes) %>% # make sure it involves a Gapin GEP
  group_by(gene) %>%
  # keep only genes that are present in >1 dataset (not interested in intersections intra-dataset)
  filter(n_distinct(dataset)>1) %>% 
  summarize(geneprograms = list(geneprogram)) %>%
  # filter(lengths(geneprograms)==2) %>%
  ggplot(aes(x = geneprograms)) +
  geom_bar() +
  scale_x_upset(n_intersections=22, sets=unique(longdf$geneprogram)) +
  labs(title="intersections with >20 genes in common")


# *********************
## 3.2. Bar plot ####

# Get all GEP genes
gapindf2 <- gather(genes.gep, key=geneprogram, value=gene, colnames(genes.gep)) %>%
  filter(!is.na(gene)) %>%
  # filter(geneprogram %in% c("GEP_1", "GEP_4", "GEP_5", "GEP_6", "GEP_7")) %>%
  mutate(geneprogram=sub("_", "", geneprogram),
         geneprogram=replace(geneprogram, geneprogram=="GEP1", "GEP1_MAIT"),
         geneprogram=replace(geneprogram, geneprogram=="GEP4", "GEP4_Temra"),
         geneprogram=replace(geneprogram, geneprogram=="GEP5", "GEP5_Tnaive"),
         geneprogram=replace(geneprogram, geneprogram=="GEP6", "GEP6_Tcm"))
table(gapindf2$geneprogram, useNA="ifany")
table(is.na(gapindf2$gene))


# Cycle through the GEPs (plot all overlaps with our GEPs)
plist <- list()
for (gep in unique(gapindf2$geneprogram)){
  print(gep)
  geneslist <- as.vector(gapindf2 %>% filter(geneprogram==gep) %>% pull(gene))
  print(length(geneslist))
  
  p <- longdf %>%
    as_tibble() %>%
    filter(dataset != "gapin") %>%
    filter(gene %in% geneslist) %>%
    group_by(geneprogram) %>%
    count() %>%
    ungroup() %>%
    mutate(total_gene=length(geneslist),
           prop_genes=n*100/total_gene) %>%
    ggplot(aes(x=reorder(geneprogram, -prop_genes), y=prop_genes, fill=geneprogram))+
      geom_bar(stat="identity")+
      labs(x="", y="% GEP genes found in each gene program", title=gep)+
      ylim(c(0,30))+
      scale_fill_manual(values=cols_alldatasets)+
      theme_cowplot()+
      theme(axis.text.x=element_text(angle=90, hjust=1), legend.position="none")
  plist[[gep]] <- p
}
# Plot all of them together
plot_grid(plotlist=plist, nrow=4)
# ggsave("./data/human-thymus/HumanData_22_CompareGeneLists/geneoverlapwithGEP_bars.jpeg", width=20, height=30)


# Cycle through the GEPs (get the top 5 overlaps for each GEP)
df.facet <- data.frame()
for (gep in unique(gapindf2$geneprogram)){
  print(gep)
  geneslist <- as.vector(gapindf2 %>% filter(geneprogram==gep) %>% pull(gene))
  print(length(geneslist))
  
  df <- longdf %>%
    as_tibble() %>%
    filter(dataset != "gapin") %>%
    filter(gene %in% geneslist) %>%
    group_by(geneprogram) %>%
    count() %>%
    ungroup() %>%
    rename(geneprogram2=geneprogram) %>%
    mutate(total_gene=length(geneslist),
           prop_genes=n*100/total_gene,
           geneprogram1=gep) %>%
    top_n(5,prop_genes)
  df.facet <- rbind(df.facet, df)
}

df.facet <- df.facet %>%
  mutate(geneprogram1=gsub("_.*", "", geneprogram1))

# Plot with facet
ggplot(df.facet, aes(x=reorder_within(geneprogram2, -prop_genes, geneprogram1), y=prop_genes, fill=geneprogram2))+
  geom_bar(stat="identity")+
  facet_wrap(~factor(geneprogram1, levels=paste0("GEP", 1:12)), nrow=1, scales="free_x")+
  tidytext::scale_x_reordered() +
  ylim(c(0,30))+
  scale_fill_manual(values=cols_alldatasets)+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1), legend.position="none")+
  labs(x="", y="% GEP genes found in each gene program")
# ggsave("./data/human-thymus/HumanData_22_CompareGeneLists/geneoverlapwithGEP_bars_top5.jpeg", width=15, height=8)





# ********************************
# 4. GENE SCORES COEXPRESSION ####
# ********************************

seur.human <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
print(seur.human) # 78,607 cells (it's the whole seurat object)
source("./scripts-final/colors_universal.R")
# DimPlot(seur.human, reduction="UMAP_50", group.by="new_clusters", cols = cols_integrated)


# clean-up useless columns
colnames(seur.human@meta.data)
seur.human@meta.data[,18:75] <- NULL
seur.human@meta.data[,13:16] <- NULL


# ******************************
## 4.1. Compute gene scores ####
# Get gene lists into a list (lol)
longdf2 <- rbind(gapindf2, canodf, roseCD4, roseCD8, poondf[,c("gene", "geneprogram")])
# table(is.na(longdf2$geneprogram))
# table(is.na(longdf2$gene))
geneprograms.list <- list()
for(gp in unique(longdf2$geneprogram)){
  print(gp)
  genes_in_program <- longdf2 %>%
    filter(geneprogram==gp & gene %in% rownames(seur.human)) %>%
    slice_head(n=200) %>%
    pull(gene)
  print(length(genes_in_program))
  geneprograms.list[[gp]] <- genes_in_program
}
# print(lengths(geneprograms.list))
# Compute cell scores
names(geneprograms.list) <- gsub(" \\| ", "_", names(geneprograms.list))
seur.human <- AddModuleScore(seur.human, name = names(geneprograms.list), features=geneprograms.list)

# Remove the annoying numbers that are being added
colnames(seur.human@meta.data)[14:22] <- stringr::str_sub(colnames(seur.human@meta.data)[14:22], end=-2)
colnames(seur.human@meta.data)[23:49] <- stringr::str_sub(colnames(seur.human@meta.data)[23:49], end=-3)

# Sanity check GEPs
SCpubr::do_FeaturePlot(seur.human, features=colnames(seur.human@meta.data)[14:25], reduction="UMAP_50", viridis_color_map = "B", order=T, ncol=3)
# ggsave("./data/human-thymus/HumanData_22_CompareGeneLists/umaps/umap_allgeps.jpeg", width=20, height=30)

# Checkpoint save
# saveRDS(seur.human, "./data/human-thymus/HumanData_22_CompareGeneLists/seuratobj_gepscores_top200genes.rds")
# seur.human <- readRDS("./data/human-thymus/HumanData_22_CompareGeneLists/seuratobj_gepscores_top200genes.rds")


# ********************************************
## 4.2. Plot co-expression on scatterplot ####
library(GGally)

test <- seur.human@meta.data %>%
  as_tibble() %>%
  pivot_longer(cols=colnames(seur.human@meta.data)[14:25], names_to="gep", values_to="gep_score") %>%
  pivot_longer(cols=colnames(seur.human@meta.data)[26:49], names_to="genexprogram", values_to="genexprogram_score")

ggplot(test,
       aes(x=gep_score, y=genexprogram_score))+
  facet_grid(genexprogram ~ gep, scales="free")+
  # facet_wrap(~genexprogram)+
  geom_point(aes(color=new_clusters), size=0.1)+
  scale_color_manual(values=cols_integrated)+
  theme_classic()
ggsave("./data/human-thymus/HumanData_22_CompareGeneLists/ggpairs_scatters/all_scatters.jpeg", width=20, height=50)


ggplot(seur.human@meta.data, aes(x=GEP1_MAIT, y=Poon_CD8MAIT))+
  geom_point(aes(color=new_clusters))+
  scale_color_manual(values=cols_integrated)+
  theme_classic()




# ******************************************
## 4.3. Functions to plot co-expression ####

# Set color matrix
BlendMatrix <- function(
    n = 10,
    col.threshold = 0.5,
    two.colors = c("#ff0000", "#00ff00"),
    negative.color = "black"
) {
  if (0 > col.threshold || col.threshold > 1) {
    stop("col.threshold must be between 0 and 1")
  }
  C0 <- as.vector(col2rgb(negative.color, alpha = TRUE))
  C1 <- as.vector(col2rgb(two.colors[1], alpha = TRUE))
  C2 <- as.vector(col2rgb(two.colors[2], alpha = TRUE))
  blend_alpha <- (C1[4] + C2[4])/2
  C0 <- C0[-4]
  C1 <- C1[-4]
  C2 <- C2[-4]
  merge.weight <- min(255 / (C1 + C2 +  C0 + 0.01))
  sigmoid <- function(x) {
    return(1 / (1 + exp(-x)))
  }
  blend_color <- function(
    i,
    j,
    col.threshold,
    n,
    C0,
    C1,
    C2,
    alpha,
    merge.weight
  ) {
    c.min <- sigmoid(5 * (1 / n - col.threshold))
    c.max <- sigmoid(5 * (1 - col.threshold))
    c1_weight <- sigmoid(5 * (i / n - col.threshold))
    c2_weight <- sigmoid(5 * (j / n - col.threshold))
    c0_weight <-  sigmoid(5 * ((i + j) / (2 * n) - col.threshold))
    c1_weight <- (c1_weight - c.min) / (c.max - c.min)
    c2_weight <- (c2_weight - c.min) / (c.max - c.min)
    c0_weight <- (c0_weight - c.min) / (c.max - c.min)
    C1_length <- sqrt(sum((C1 - C0) ** 2))
    C2_length <- sqrt(sum((C2 - C0) ** 2))
    C1_unit <- (C1 - C0) / C1_length
    C2_unit <- (C2 - C0) / C2_length
    C1_weight <- C1_unit * c1_weight
    C2_weight <- C2_unit * c2_weight
    C_blend <- C1_weight * (i - 1) * C1_length / (n - 1) + C2_weight * (j - 1) * C2_length / (n - 1) + (i - 1) * (j - 1) * c0_weight * C0 / (n - 1) ** 2 + C0
    C_blend[C_blend > 255] <- 255
    C_blend[C_blend < 0] <- 0
    return(rgb(
      red = C_blend[1],
      green = C_blend[2],
      blue = C_blend[3],
      alpha = alpha,
      maxColorValue = 255
    ))
  }
  blend_matrix <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      blend_matrix[i, j] <- blend_color(
        i = i,
        j = j,
        col.threshold = col.threshold,
        n = n,
        C0 = C0,
        C1 = C1,
        C2 = C2,
        alpha = blend_alpha,
        merge.weight = merge.weight
      )
    }
  }
  return(blend_matrix)
}

# Plot color matrix
Melt <- function(x) {
  if (!is.data.frame(x = x)) {
    x <- as.data.frame(x = x)
  }
  return(data.frame(
    rows = rep.int(x = rownames(x = x), times = ncol(x = x)),
    cols = unlist(x = lapply(X = colnames(x = x), FUN = rep.int, times = nrow(x = x))),
    vals = unlist(x = x, use.names = FALSE)
  ))
}

BlendMap <- function(color.matrix, step=2, xtext='rows', ytext="cols") {
  color.heat <- matrix(
    data = 1:prod(dim(x = color.matrix)) - 1,
    nrow = nrow(x = color.matrix),
    ncol = ncol(x = color.matrix),
    dimnames = list(
      1:nrow(x = color.matrix),
      1:ncol(x = color.matrix)
    )
  )
  
  # xbreaks <- seq.int(from = 0, to = nrow(x = color.matrix), by = step)
  # ybreaks <- seq.int(from = 0, to = ncol(x = color.matrix), by = step)
  color.heat <- Melt(x = color.heat)
  color.heat$rows <- as.numeric(x = as.character(x = color.heat$rows))
  color.heat$cols <- as.numeric(x = as.character(x = color.heat$cols))
  color.heat$vals <- factor(x = color.heat$vals)
  plot <- ggplot(
    data = color.heat,
    mapping = aes_string(x = "rows", y = "cols", fill = 'vals')
  ) +
    geom_raster(show.legend = FALSE) +
    theme(plot.margin = unit(x = rep.int(x = 0, times = 4), units = 'cm')) +
    # scale_x_continuous(breaks = xbreaks, expand = c(0, 0), labels = xbreaks) +
    # scale_y_continuous(breaks = ybreaks, expand = c(0, 0), labels = ybreaks) +
    scale_fill_manual(values = as.vector(x = color.matrix)) +
    labs(x=xtext, y=ytext)+
    theme_cowplot()
  
  if(step!=0){
    xbreaks <- seq.int(from = 0, to = nrow(x = color.matrix), by = step)
    ybreaks <- seq.int(from = 0, to = ncol(x = color.matrix), by = step)
    plot <- plot +
      scale_x_continuous(breaks = xbreaks, expand = c(0, 0), labels = xbreaks) +
      scale_y_continuous(breaks = ybreaks, expand = c(0, 0), labels = ybreaks)
  } else if(step==0){
    plot <- plot+theme(axis.text=element_blank(), axis.ticks=element_blank())
  }
  
  return(plot)
}

# Normalize expression level of 2 features & return also a "blend" expression (corresponding to color matrix)
BlendExpression <- function(data, nlevels=100) {
  if (ncol(x = data) != 2) {
    stop("'BlendExpression' only blends two features")
  }
  features <- colnames(x = data)
  data <- as.data.frame(x = apply(
    X = data,
    MARGIN = 2,
    FUN = function(x) {
      return(round(x = (nlevels-1) * (x - min(x)) / (max(x) - min(x))))
    }
  ))
  data[, 3] <- data[, 1] + data[, 2] * nlevels
  # colnames(x = data) <- c(features, paste(features, collapse = '_'))
  colnames(x = data) <- c(features, "blend")
  for (i in 1:ncol(x = data)) {
    data[, i] <- factor(x = data[, i])
  }
  return(data)
}

PlotCoexpression <- function(seuratobj,
                             features,
                             plotting="blend",
                             pwithmatrix=T,
                             nlevels=10,
                             cols.neg="#737373", cols.pos=c("#74c476", "#fd8d3c"), col.threshold=0.5, colmatrix_stepsize=10,
                             order=T){
  # GET COLOR MATRIX
  cat("\n-> Getting color matrix\n")
  color.matrix <- BlendMatrix(
    two.colors = cols.pos,
    col.threshold = col.threshold,
    negative.color = cols.neg,
    n=nlevels
  )
  
  # DEFINE COLOR LIST FOR PLOTTING
  cat("\n-> Defining colors for plotting\n")
  colors <- list(
    color.matrix[, 1], # red
    color.matrix[1, ], # green
    as.vector(x = color.matrix)
  )
  
  # BLEND EXPRESSION
  cat("\n-> Blending features expression\n")
  df <- seuratobj@meta.data[, features]
  df <- BlendExpression(df, nlevels=nlevels) # 3 columns
  # head(df)
  # GET PLOTTING DATAFRAME
  cat("\n-> Defining plotting DF\n")
  dims <- seuratobj@reductions$UMAP_50@cell.embeddings
  # head(dims)
  df_final <- cbind(df, dims, seuratobj@meta.data[,"new_clusters"])
  colnames(df_final)[6] <- "clusters"
  # head(df_final)
  
  # PLOT
  if(plotting=="feature1"){
    cat("\n-> Plotting feature 1\n")
    if(order==T){df_final <- df_final[order(df_final[,1]),]}
    df_final[,1] <- as.numeric(as.character(df_final[,1])) # transform factors to numbers for plotting
    p <- ggplot(df_final, aes_string(x=colnames(dims)[1], y=colnames(dims)[2], color=colnames(df_final)[1]))+
      geom_point(size=0.1)+
      scale_color_gradient(low=color.matrix[1, 1], high=color.matrix[nrow(color.matrix), 1])
  }
  else if(plotting=="feature2"){
    cat("\n-> Plotting feature 2\n")
    if(order==T){df_final <- df_final[order(df_final[,2]),]}
    df_final[,2] <- as.numeric(as.character(df_final[,2])) # transform factors to numbers for plotting
    p <- ggplot(df_final, aes_string(x=colnames(dims)[1], y=colnames(dims)[2], color=colnames(df_final)[2]))+
      geom_point(size=0.1)+
      scale_color_gradient(low=color.matrix[1, 1], high=color.matrix[1, ncol(color.matrix)])
  }
  else if(plotting=="blend"){
    cat("\n-> Plotting blended features\n")
    if(order==T){df_final <- df_final[order(df_final[,1], df_final[,2]),]} # order points by increasing value of score 1 and score 2
    # df_final[,3] <- as.numeric(as.character(df_final[,3])) # transform factors to numbers for plotting
    # Colors
    cols.use <- as.vector(color.matrix)
    names(cols.use) <- as.character(0:(length(cols.use)-1))
    # Plot
    p <- ggplot(df_final, aes_string(x=colnames(dims)[1], y=colnames(dims)[2], color=colnames(df_final)[3]))+
      geom_point(size=0.3)+
      scale_color_manual(values=cols.use)+
      labs(x="UMAP1", y="UMAP2")+
      theme_cowplot()+
      theme(legend.position="none",
            axis.text=element_blank(),
            axis.title=element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.border=element_rect(color="white", fill=NA, size=1))
    if(pwithmatrix==T){
      cat("\n-> Adding color matrix on plot\n")
      p <- ggdraw(p)+
        draw_plot(BlendMap(color.matrix, step=colmatrix_stepsize, xtext=features[1], ytext=features[2]),
                  0.05,0.06,.25,.25)
    }
  }
  
  return(p)
}




# *****************************************************
## 4.4. Plot co-expression of gene scores on UMAPs ####

# Sanity check everything yellow
# PlotCoexpression(seuratobj=seur.human, features=c("GEP1_MAIT", "GEP1_MAIT"), nlevels=100, colmatrix_stepsize=0)

# GEPs to loop through
gepsall <- colnames(seur.human@meta.data)[14:49]

# GEP1
for(gep in gepsall[gepsall!="GEP1_MAIT"]){
  filename <- paste0("umap_gep1_with_", gsub("\\.", "_", gep), ".jpeg")
  print(filename)
  p <- PlotCoexpression(seuratobj=seur.human, features=c("GEP1_MAIT", gep), nlevels=100, colmatrix_stepsize=0)
  ggsave(filename=file.path("./data/human-thymus/HumanData_22_CompareGeneLists/umaps/gep1_mait", filename), plot=p, width=9, height=8)
}

# GEP4
for(gep in gepsall[gepsall!="GEP4_Temra"]){
  filename <- paste0("umap_gep4_with_", gsub("\\.", "_", gep), ".jpeg")
  print(filename)
  p <- PlotCoexpression(seuratobj=seur.human, features=c("GEP4_Temra", gep), nlevels=100, colmatrix_stepsize=0)
  ggsave(filename=file.path("./data/human-thymus/HumanData_22_CompareGeneLists/umaps/gep4_temra", filename), plot=p, width=9, height=8)
}

# GEP6
for(gep in gepsall[gepsall!="GEP6_Tcm"]){
  filename <- paste0("umap_gep6_with_", gsub("\\.", "_", gep), ".jpeg")
  print(filename)
  p <- PlotCoexpression(seuratobj=seur.human, features=c("GEP6_Tcm", gep), nlevels=100, colmatrix_stepsize=0)
  ggsave(filename=file.path("./data/human-thymus/HumanData_22_CompareGeneLists/umaps/gep6_tcm", filename), plot=p, width=9, height=8)
}


# UMAPs of interest
colnames(seur.human@meta.data)[14] <- "GEP1"
colnames(seur.human@meta.data)[17] <- "GEP4"
colnames(seur.human@meta.data)[19] <- "GEP6"
# GEP1
ggsave(plot=PlotCoexpression(seuratobj=seur.human, features=c("GEP1", "CanoGamez_CD4_TEM"), nlevels=100, colmatrix_stepsize=0),
       filename="./data/human-thymus/HumanData_22_CompareGeneLists/umaps/final_umaps/gep1_canogamez_CD4TEM.pdf", width=9, height=8)
ggsave(plot=PlotCoexpression(seuratobj=seur.human, features=c("GEP1", "Rose_CD8_modul3_Tem.emra"), nlevels=100, colmatrix_stepsize=0),
       filename="./data/human-thymus/HumanData_22_CompareGeneLists/umaps/final_umaps/gep1_rose_CD8modul3.pdf", width=9, height=8)
ggsave(plot=PlotCoexpression(seuratobj=seur.human, features=c("GEP1", "Poon_CD8MAIT"), nlevels=100, colmatrix_stepsize=0),
       filename="./data/human-thymus/HumanData_22_CompareGeneLists/umaps/final_umaps/gep1_poon_CD8mait.pdf", width=9, height=8)

# GEP4
ggsave(plot=PlotCoexpression(seuratobj=seur.human, features=c("GEP4", "CanoGamez_CD4_TEMRA"), nlevels=100, colmatrix_stepsize=0),
       filename="./data/human-thymus/HumanData_22_CompareGeneLists/umaps/final_umaps/gep4_canogamez_CD4TEMRA.pdf", width=9, height=8)
ggsave(plot=PlotCoexpression(seuratobj=seur.human, features=c("GEP4", "Rose_CD8_modul1_Temra"), nlevels=100, colmatrix_stepsize=0),
       filename="./data/human-thymus/HumanData_22_CompareGeneLists/umaps/final_umaps/gep4_rose_CD8modul1.pdf", width=9, height=8)
ggsave(plot=PlotCoexpression(seuratobj=seur.human, features=c("GEP4", "Poon_CD8TEM.TEMRA"), nlevels=100, colmatrix_stepsize=0),
       filename="./data/human-thymus/HumanData_22_CompareGeneLists/umaps/final_umaps/gep4_poon_CD8tem_temra.pdf", width=9, height=8)

# GEP6
ggsave(plot=PlotCoexpression(seuratobj=seur.human, features=c("GEP6", "CanoGamez_CD4_TCM"), nlevels=100, colmatrix_stepsize=0),
       filename="./data/human-thymus/HumanData_22_CompareGeneLists/umaps/final_umaps/gep6_canogamez_CD4Tcm.pdf", width=9, height=8)
ggsave(plot=PlotCoexpression(seuratobj=seur.human, features=c("GEP6", "Rose_CD8_modul2_Tcm.em"), nlevels=100, colmatrix_stepsize=0),
       filename="./data/human-thymus/HumanData_22_CompareGeneLists/umaps/final_umaps/gep6_rose_CD8modul2.pdf", width=9, height=8)
ggsave(plot=PlotCoexpression(seuratobj=seur.human, features=c("GEP6", "Poon_CD4TCM.TFH"), nlevels=100, colmatrix_stepsize=0),
       filename="./data/human-thymus/HumanData_22_CompareGeneLists/umaps/final_umaps/gep6_poon_CD4Tcm.pdf", width=9, height=8)

