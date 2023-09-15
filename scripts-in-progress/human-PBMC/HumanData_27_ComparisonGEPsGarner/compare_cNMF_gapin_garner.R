# Purpose: Compare GEPs from our dataset & from Garner stimulation dataset
# Author: Salomé Carcy
# Date: September 2023


# **************
# 1. IMPORT ####
# **************

# Import librairies
library(Seurat)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(tidytext)

# Import data
setwd("~/Projects/HumanThymusProject/")
cnmf_garner <- read.csv("./data/human-PBMC/HumanData_27_ComparisonGEPsGarner/garner_genes_per_gep_post_rank_threshold_k4.csv", row.names=1)
cnmf_gapin  <- read.csv("./data/human-PBMC/HumanData_27_ComparisonGEPsGarner/limited_nonimputed_genes_per_gep_post_rank_threshold_k12.csv", row.names=1)
head(cnmf_garner)
head(cnmf_gapin)

# seurat objects
seur.human <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_08_28_23.RDS")
seur.pbmc <- subset(seur.human, Tissue=="PBMC")
print(seur.pbmc) # 17,204 genes and 40,986 cells

seur.garner <- readRDS("./data/human-PBMC/HumanData_27_ComparisonGEPsGarner/GSE238138_20h_seurat.rds")
SCpubr::do_DimPlot(seur.garner, reduction="umap", group.by="stimulation", label=T, legend.position="none")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_27_ComparisonGEPsGarner/plots/umap_garner_stim.jpeg", width=8, height=7)




# *****************
# 2. FUNCTIONS ####
# *****************

# **************************************
## 2.1. Compute gene lists overlaps ####

# Define Jaccard function
overlapcoef <- function(a, b, coef="overlapcoef") {
  if(coef=="jaccard"){
    # cat("\n- Computing jaccard -\n")
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
  }
  else if(coef=="jaccardweight"){
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    maxj = min(length(a),length(b))/max(length(a), length(b))
    return (intersection/(union*maxj))
  }
  else if(coef=="overlapcoef"){
    # cat("\n- Computing overlap coefficient -\n")
    intersection = length(intersect(a, b))
    denominator = min(length(a),length(b))
    return (intersection/denominator)
  }
}

# Define function to compute similarity coefficient (overlap, jaccard, ...) between pairs of gene lists.
# H0: draw 1000 random gene lists representative of list A and list B (1000 for each), compute their pairwise coefficient
NullOverlap_double <- function(seuratobj=seur.pbmc,
                               df_geneprograms=dflong,
                               geneprograms_to_test=unique(dflong$geneprogram), # by default GEPs
                               coefficient_to_compute="jaccardweight",
                               nbins=25,
                               nrandom=1000){
  
  # ___________________________________________________________
  # -- 1. Get DF with all genes & binned by expression level --
  cat("\n-- 1. GET ALL GENES & BIN THEM BY EXPRESSION LEVEL --\n")
  allgenes_binned_DF <- data.frame("gene"=rownames(seuratobj),
                                   "totalexpression"=Matrix::rowSums(seuratobj@assays$RNA@data))
  # nrow(allgenes_binned_DF) == nrow(seuratobj) # sanity check
  # table(rownames(allgenes_binned_DF)==allgenes_binned_DF$gene) # sanity check
  # hist(allgenes_binned_DF$totalexpression, breaks=100) # visualize distribution
  
  # Get the bins by percentiles
  allgenes_binned_DF <- allgenes_binned_DF[order(allgenes_binned_DF$totalexpression),]
  allgenes_binned_DF$bin <- ggplot2::cut_number(x = allgenes_binned_DF$totalexpression, n = nbins, labels = FALSE, right = FALSE)
  cat("Number of genes in each bin, from bin 1 to bin", nbins, ":\n")
  print(table(allgenes_binned_DF$bin))
  # ggplot(allgenes_binned_DF)+
  #   geom_density(aes(x=totalexpression, group=bin, color=factor(bin, levels=1:10))) # sanity check
  
  # _____________________________________________________
  # -- 2. Draw random ctrl genes for each gene program --
  cat("\n-- 2. DRAW RANDOM CTRL GENES FOR EACH GENE SET --\n")
  randomgenesets_list <- list()
  
  for(refgeneprog in geneprograms_to_test){
    cat("\n++ Drawing random ctrl genes for", refgeneprog, "\n")
    # Find in which expression bins are our genes from the gene set
    refgeneset <- df_geneprograms[df_geneprograms$geneprogram==refgeneprog, "gene"]
    refgeneset_binned_DF <- allgenes_binned_DF[allgenes_binned_DF$gene %in% refgeneset,]
    # print(nrow(refgeneset_binned_DF)==length(refgeneset))
    # print(paste0("Length gene set:", length(unique(refgeneset))))
    # print(paste0("Nb rows DF:", nrow(refgeneset_binned_DF)))
    
    # ++++++++++++++++++++++++++++++++++
    # Do 1000 random draws of gene lists
    ctrlgenematrix <- matrix(nrow=length(refgeneset), ncol=nrandom)
    for(i in 1:nrandom){
      if(i%%100==0){cat(paste0("\nrandom draw #", i))} # progress bar
      # cat(paste0("\nrandom draw #", i))
      # Loop by expression bin & sample equal number of genes from each bin
      ctrlgeneset <- c()
      for(bin in unique(refgeneset_binned_DF$bin)){
        # print(bin)
        ngenes_to_sample <- sum(refgeneset_binned_DF$bin==bin)
        # print(ngenes_to_sample)
        ctrlgenes_in_same_bin <- sample(x=allgenes_binned_DF[allgenes_binned_DF$bin==bin,"gene"], size=ngenes_to_sample, replace=FALSE)
        # print(length(ctrlgenes_in_same_bin))
        ctrlgeneset <- c(ctrlgeneset, ctrlgenes_in_same_bin)
      }
      # Sanity check
      if(length(ctrlgeneset)!=length(refgeneset)){cat("\nPROBLEM: control gene set of size", length(ctrlgeneset), "while gene program length is", length(refgeneset))}
      
      # Add ctrlgeneset to matrix
      ctrlgenematrix[,i] <- ctrlgeneset
    }
    # Sanity check
    cat("\n\n++ Selection of random ctrl genes over! ++\n")
    cat("Length of reference gene set is", length(unique(refgeneset)), "and length of the random gene sets:", nrow(ctrlgenematrix), "\n")
    # ++++++++++++++++++++++++++++++++++
    
    # Add to the list
    randomgenesets_list[[refgeneprog]] <- ctrlgenematrix
  }
  
  # _____________________________________________________________________________
  # -- 3. Compute random & observed overlap between each pairwise gene program --
  cat("\n-- 3. COMPUTE RANDOM AND OBSERVED OVERLAP --\n")
  
  combinations <- t(combn(geneprograms_to_test,2))
  colnames(combinations) <- c("geneprogram1", "geneprogram2")
  
  # add self combinations
  combinations <- rbind(combinations,
                        t(matrix(rep(geneprograms_to_test, 2), nrow=2, byrow=T)))
  combinations <- as.data.frame(combinations)
  ncombinations <- nrow(combinations)
  
  for(i in 1:nrow(combinations)){
    # cat("\n", combinations[i,1], "vs", combinations[i,2])
    if(i%%25==0){cat(paste0("\ncomparison ", i, "/", ncombinations))} # progress bar
    
    # Get observed overlap
    genesetA <- df_geneprograms %>% filter(geneprogram==combinations[i,"geneprogram1"]) %>% pull(gene)
    genesetB <- df_geneprograms %>% filter(geneprogram==combinations[i,"geneprogram2"]) %>% pull(gene)
    observedpercent <- overlapcoef(genesetA, genesetB, coef=coefficient_to_compute)
    
    # Get vector of random overlaps
    randomgenesetA <- randomgenesets_list[[combinations[i,"geneprogram1"]]]
    randomgenesetB <- randomgenesets_list[[combinations[i,"geneprogram2"]]]
    ctrlpercent <- sapply(1:nrandom, function(x) overlapcoef(randomgenesetA[,x], randomgenesetB[,x], coef=coefficient_to_compute))
    if(length(ctrlpercent) != nrandom){cat("\nPROBLEM: vector of random overlaps of length", length(ctrlpercent), "when there should be", nrandom, "random sets")}
    
    # save info
    combinations[i,"observedoverlap"] <- observedpercent
    combinations[i,"randomoverlap_mean"] <- mean(ctrlpercent)
    combinations[i,"randomoverlap_min"] <- min(ctrlpercent)
    combinations[i,"randomoverlap_max"] <- max(ctrlpercent)
    combinations[i,"pval"] <- sum(ctrlpercent>observedpercent)/nrandom
  }
  
  return(combinations)
}


# *****************************
## 2.2. Plot co-expression ####

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
                             dataset="gapin",
                             pwithmatrix=T,
                             rasterdpi=300,
                             nlevels=100,
                             cols.neg="#969696", cols.pos=c("#74c476", "#fd8d3c"), col.threshold=0.5, colmatrix_stepsize=10,
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
  if(dataset=="gapin"){
    dims <- seuratobj@reductions$UMAP_50@cell.embeddings
    df_final <- cbind(df, dims, seuratobj@meta.data[,"new_clusters"])
  }
  else if(dataset=="garner"){
    dims <- seuratobj@reductions$umap@cell.embeddings
    df_final <- cbind(df, dims, seuratobj@meta.data[,"stimulation"])
  }
  # head(dims)
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
    # rasterise to avoid computer crashing
    p <- ggrastr::rasterise(p, layers="Point", dpi=rasterdpi)
    # add color matrix if specified
    if(pwithmatrix==T){
      cat("\n-> Adding color matrix on plot\n")
      p <- ggdraw(p)+
        draw_plot(BlendMap(color.matrix, step=colmatrix_stepsize, xtext=features[1], ytext=features[2]),
                  0.05,0.06,.25,.25)
    }
  }
  
  return(p)
}




# ****************************
# 3. GENE LIST COMPARISON ####
# ****************************

## 3.1. Transform GEP dataframes into a long df ####

# Gapin
head(cnmf_gapin)
dflong_gapin <- gather(cnmf_gapin, key=geneprogram, value=gene, colnames(cnmf_gapin)) %>%
  filter(!is.na(gene)) %>%
  mutate(geneprogram=paste0(geneprogram, "_gapin"))
dflong_gapin$dataset <- "gapin"
head(dflong_gapin)
table(dflong_gapin$geneprogram, useNA="ifany")
table(is.na(dflong_gapin$gene))

# Garner
head(cnmf_garner)
dflong_garner <- gather(cnmf_garner, key=geneprogram, value=gene, colnames(cnmf_garner)) %>%
  filter(!is.na(gene)) %>%
  mutate(geneprogram=paste0(geneprogram, "_garner"))
dflong_garner$dataset <- "garner"
head(dflong_garner)
table(dflong_garner$geneprogram, useNA="ifany")
table(is.na(dflong_garner$gene))

# Bind rows
dflong <- rbind(dflong_gapin, dflong_garner) %>%
  filter(gene %in% rownames(seur.pbmc))
table(dflong$dataset, useNA="ifany")
# gapin   garner
# 4,068    1,994


## 3.2. Get weighted jaccard index between every pairwise GEP ####
geps_overlap_wjaccard <- NullOverlap_double(seuratobj=seur.pbmc,
                                            df_geneprograms = dflong,
                                            geneprograms_to_test=unique(dflong$geneprogram),
                                            coefficient_to_compute="jaccardweight",
                                            nbins=25, nrandom=1000)

geps_overlap_wjaccard %>%
  as_tibble() %>%
  # keep geps of interest
  filter(geneprogram1 %in% unique(grep("_gapin", geneprogram1, value=T))) %>%
  filter(geneprogram2 %in% unique(grep("_garner", geneprogram2, value=T))) %>%
  # add padj and keep only other dataset's programs that have at least one x% overlap
  group_by(geneprogram1) %>%
  mutate(padj=pval*n_distinct(geneprogram2)) %>%
  ungroup() %>%
  # filter(max_observedoverlap>0.333) %>%
  # distinct() %>%
  # last details (color, padj to plot)
  # left_join(cols_df, by=join_by("geneprogram2"=="geneprogram")) %>%
  mutate(padj_toplot=ifelse(padj==0, "< 0.001",
                            ifelse(padj<0.05, as.character(round(padj, 2)), ""))) %>%
  # PLOT
  ggplot(aes(x=factor(geneprogram1, levels=rev(paste0("GEP", 1:12, "_gapin"))), y=observedoverlap))+
  geom_bar(stat="identity")+
  geom_point(aes(y=randomoverlap_mean), shape="|", size=3)+
  tidytext::scale_x_reordered() +
  geom_text(aes(label=padj_toplot), size=4, angle=0, hjust=-0.2)+
  # facet_grid(geneprogram_cat~geneprogram2, space="free", scales="free")+
  facet_wrap(~factor(geneprogram2, levels=paste0("GEP", 1:4, "_garner")), ncol=6)+
  # facet_manual(~factor(geneprogram2, levels=order_programs), design=rbind(c(1:6), c(7:9,NA,NA,NA), c(10:15), c(16:19,NA,NA)), scales="free_x")+
  ylim(c(0,0.5))+
  coord_flip()+
  theme_cowplot()+
  # scale_fill_manual(values=cols_geneprogcat, name="")+
  labs(y="Weighted Jaccard Index", x="")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_27_ComparisonGEPsGarner/plots/genelists_comparison_wjaccard.jpeg", width=10, height=4)
 

## 3.3. Get list of common genes between GEP5 and mait one
genes_common_gep5_gep1garner <- intersect(dflong[dflong$geneprogram=="GEP5_gapin","gene"], dflong[dflong$geneprogram=="GEP1_garner", "gene"])
length(genes_common_gep5_gep1garner) # 38
genes_common_gep5_gep3garner <- intersect(dflong[dflong$geneprogram=="GEP5_gapin","gene"], dflong[dflong$geneprogram=="GEP3_garner", "gene"])
length(genes_common_gep5_gep3garner) # 35
intersect(genes_common_gep5_gep1garner, genes_common_gep5_gep3garner)

seur.human <- AddModuleScore(seur.human, name=c("maitcommongenes_unstim", "maitcommongenes_cytok"),
                             features=list("maitcommongenes_unstim"=genes_common_gep5_gep1garner,
                                           "maitcommongenes_cytok"=genes_common_gep5_gep3garner),
                             seed=1)
colnames(seur.human@meta.data)[66:67] <- c("maitcommongenes_unstim", "maitcommongenes_cytok")
SCpubr::do_FeaturePlot(seur.human,  features=c("maitcommongenes_unstim", "maitcommongenes_cytok"), reduction="UMAP_50", order=T, ncol=2)
VlnPlot(seur.human, features="maitcommongenes_unstim", group.by="new_clusters")+
  labs(x="", y="gene list score", title="Score of common genes btw GEP5_gapin & unstim GEP1_garner") |
VlnPlot(subset(seur.human, subset= new_clusters %in% 13:14), features="maitcommongenes_unstim", group.by="group.ident")+
  labs(x="", y="gene list score", title="Cells in clusters 13-14 only")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_27_ComparisonGEPsGarner/plots/commongenes_gep5gapin_gep1garner_scoring.pdf", width=15, height=6)



# *************************
# 4. GEPs COEXPRESSION ####
# *************************

# ******************************
## 4.1. Compute gene scores ####

# Get gene lists into a list (lol)
# table(is.na(dflong$geneprogram))
# table(is.na(dflong$gene))
geneprograms_list <- list()
for(gp in unique(dflong$geneprogram)){
  print(gp)
  geneprograms_list[[gp]] <- dflong %>%
    filter(geneprogram==gp & gene %in% rownames(seur.human) & gene %in% rownames(seur.garner)) %>%
    pull(gene)
}
print(lengths(geneprograms_list))

# Keep only genes that are among HVGs
geneprograms_gapinhvg_list <- lapply(geneprograms_list, function(x) x[x %in% VariableFeatures(seur.human)])
names(geneprograms_gapinhvg_list) <- gsub("_", "hvg_", names(geneprograms_gapinhvg_list))
print(lengths(geneprograms_gapinhvg_list))

DefaultAssay(seur.garner) <- "RNA"
seur.garner <- FindVariableFeatures(seur.garner)
hvg_garner_RNA <- VariableFeatures(seur.garner)
geneprograms_garnerhvgRNA_list <- lapply(geneprograms_list, function(x) x[x %in% hvg_garner_RNA])
names(geneprograms_garnerhvgRNA_list) <- gsub("_", "hvg_", names(geneprograms_garnerhvgRNA_list))
print(lengths(geneprograms_garnerhvgRNA_list))

# Compute cell scores
seur.human  <- AddModuleScore(seur.human,  name = names(geneprograms_list), features=geneprograms_list, seed=1)
seur.garner <- AddModuleScore(seur.garner, name = names(geneprograms_list), features=geneprograms_list, seed=1, assay="SCT")
# seur.garner <- AddModuleScore(seur.garner, name = paste0(names(geneprograms_list), "_RNA"), features=geneprograms_list, seed=1, assay="RNA")
seur.human  <- AddModuleScore(seur.human,  name = names(geneprograms_gapinhvg_list), features=geneprograms_gapinhvg_list, seed=1)
seur.garner <- AddModuleScore(seur.garner, name = paste0(names(geneprograms_garnerhvgRNA_list), "_RNA"), features=geneprograms_garnerhvgRNA_list, seed=1, assay="RNA")

# Remove the annoying numbers that are being added
colnames(seur.human@meta.data)[50:58] <- stringr::str_sub(colnames(seur.human@meta.data)[50:58], end=-2)
colnames(seur.human@meta.data)[59:65] <- stringr::str_sub(colnames(seur.human@meta.data)[59:65], end=-3)
colnames(seur.garner@meta.data)[36:44] <- stringr::str_sub(colnames(seur.garner@meta.data)[36:44], end=-2)
colnames(seur.garner@meta.data)[45:51] <- stringr::str_sub(colnames(seur.garner@meta.data)[45:51], end=-3)

colnames(seur.human@meta.data)[68:76] <- stringr::str_sub(colnames(seur.human@meta.data)[68:76], end=-2)
colnames(seur.human@meta.data)[77:83] <- stringr::str_sub(colnames(seur.human@meta.data)[77:83], end=-3)
colnames(seur.garner@meta.data)[68:76] <- stringr::str_sub(colnames(seur.garner@meta.data)[68:76], end=-2)
colnames(seur.garner@meta.data)[77:83] <- stringr::str_sub(colnames(seur.garner@meta.data)[77:83], end=-3)

# Sanity check GEPs
SCpubr::do_FeaturePlot(seur.human,  features=colnames(seur.human@meta.data)[50:60],  reduction="UMAP_50", order=T, ncol=3)
SCpubr::do_FeaturePlot(seur.garner, features=colnames(seur.garner@meta.data)[48:51], reduction="umap", order=T, ncol=2)


# *****************************************************
## 4.2. Plot GEP coexpression on garner seurat obj ####

SCpubr::do_FeaturePlot(seur.garner, features=paste0("GEP", 1:11, "_gapin_RNA"), reduction="umap", order=T, ncol=3)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_27_ComparisonGEPsGarner/plots/umap_garner_GEPsgapin_RNAscore.jpeg", width=20, height=28)
SCpubr::do_FeaturePlot(seur.garner, features=paste0("GEP", 1:11, "hvg_gapin_RNA"), reduction="umap", order=T, ncol=3, plot.title = "HVG only within GEPs")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_27_ComparisonGEPsGarner/plots/umap_garner_GEPsgapin_RNAscore_hvg.jpeg", width=20, height=28)

# GEP3_gapin (naive) & GEP1_garner (unstim MAIT blood)
SCpubr::do_FeaturePlot(seur.garner, features=c("GEP3_gapin", "GEP1_garner"), reduction="umap", order=T, ncol=2)
PlotCoexpression(seuratobj=seur.garner, dataset="garner", colmatrix_stepsize=0, rasterdpi=300, features=c("GEP3_gapin", "GEP1_garner"))

# GEP5_gapin (Th1/Th17) & GEP1_garner (unstim MAIT blood) & GEP3_garner (cytokine stim MAIT)
SCpubr::do_FeaturePlot(seur.garner, features=c("GEP5_gapin", "GEP1_garner", "GEP3_garner"), reduction="umap", order=T, ncol=3)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_27_ComparisonGEPsGarner/plots/umap_garner_GEP5gapin_GEP1garner_GEP3garner.jpeg", width=18, height=7)
PlotCoexpression(seuratobj=seur.garner, dataset="garner", colmatrix_stepsize=0, rasterdpi=300, features=c("GEP5_gapin", "GEP1_garner"))

# GEP6_gapin & GEP4_garner (TCR stim MAIT)
SCpubr::do_FeaturePlot(seur.garner, features=c("GEP6_gapin", "GEP4_garner"), reduction="umap", order=T, ncol=2)
PlotCoexpression(seuratobj=seur.garner, dataset="garner", colmatrix_stepsize=0, rasterdpi=300, features=c("GEP6_gapin", "GEP4_garner"))


# *****************************************************
## 4.3. Plot GEP coexpression on gapin seurat obj ####

SCpubr::do_FeaturePlot(seur.human, features=paste0("GEP", 1:4, "hvg_garner"), reduction="UMAP_50", order=T, ncol=2, plot.title = "HVG only within GEPs")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_27_ComparisonGEPsGarner/plots/umap_gapin_GEPsgarner_hvg.jpeg", width=12, height=14)

# GEP5_gapin (Th1/Th17) & GEP1_garner (unstim MAIT blood) & GEP3_garner (cytokine stim MAIT)
SCpubr::do_FeaturePlot(seur.human, features=c("GEP5_gapin", "GEP1_garner", "GEP3_garner"), reduction="UMAP_50", order=T, ncol=3)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_27_ComparisonGEPsGarner/plots/umap_gapin_GEP5gapin_GEP1garner_GEP3garner.jpeg", width=18, height=7)




# ********************
# 5. METANEIGHBOR ####
# ********************

# ***********************
## 5.1. Prepare data ####

# Counts gapin: keep only PBMC MAITs & the clusters
seur.mait <- readRDS("./data/human-PBMC/HumanData_26_GEPusage_per_lineage_PBMC/blood.MAIT_03_16_23.RDS")
table(seur.mait$group.ident, useNA="ifany")
table(seur.mait$new_clusters, useNA="ifany")
table(seur.mait$new_clusters_MAIT, useNA="ifany")
gapin.counts <- seur.mait@assays$RNA@counts
# gapin.counts[11:15,1:5] # sanity check raw counts
gapin.metadata <- seur.mait@meta.data[,c("new_clusters", "new_clusters_MAIT")]
colnames(gapin.metadata) <- c("mait_clusters_highres", "mait_clusters_lowres")
gapin.metadata$mait_clusters_highres <- paste0("gapin_", gapin.metadata$mait_clusters_highres)
gapin.metadata$study <- "gapin"
head(gapin.metadata)


# Counts garner
garner.counts <- seur.garner@assays$RNA@counts
# garner.counts[16:20,1:5] # sanity check raw counts
garner.metadata <- seur.garner@meta.data[,c("snn_30_res.0.3", "stimulation")]
colnames(garner.metadata) <- c("mait_clusters_highres", "mait_clusters_lowres")
garner.metadata$mait_clusters_highres <- paste0("garner_", garner.metadata$mait_clusters_highres)
garner.metadata$study <- "garner"
head(garner.metadata)


# Get HVGs
gapin.hvg <- VariableFeatures(seur.mait)
DefaultAssay(seur.garner) <- "RNA"
seur.garner[["RNA"]]@meta.features <- data.frame(row.names = rownames(seur.garner[["RNA"]]))
seur.garner <- FindVariableFeatures(seur.garner)
garner.hvg <- VariableFeatures(seur.garner)
hvg.mait <- union(gapin.hvg, garner.hvg)
length(hvg.mait) # 288 (intersect), 3712 (union)

# Subset counts to common genes to both
allgenes <- intersect(rownames(gapin.counts), rownames(garner.counts))
length(allgenes) # 16,058
gapin.counts  <- gapin.counts[allgenes,]
garner.counts <- garner.counts[allgenes,]

# Merge everything into one
seur.total <- merge(CreateSeuratObject(counts=gapin.counts, meta.data=gapin.metadata),
                    CreateSeuratObject(counts=garner.counts, meta.data=garner.metadata)) # 16,058 genes

# Sanity checks
head(seur.total@meta.data)
table(seur.total$mait_clusters_lowres, useNA="ifany")
table(seur.total$mait_clusters_highres, useNA="ifany")
table(seur.total$study, useNA="ifany")


# Convert seurat count matrix to SummarizedExperiment object
library(SummarizedExperiment)
se.total <- SummarizedExperiment(assays=seur.total@assays[["RNA"]]@counts,
                                 colData=seur.total$mait_clusters_lowres)


# _______________________
# METANEIGHBOR
# Run metaneighbor
library(MetaNeighbor)
mtn <- MetaNeighborUS(var_genes=hvg.mait,
                      dat=se.total,
                      study_id=seur.total$study,
                      cell_type=seur.total$mait_clusters_lowres,
                      fast_version=TRUE)


# ****************************
## 5.2. Plot metaneighbor ####

# With everything (check diagonal)
library(gplots)
# jpeg("./scripts-in-progress/human-PBMC/HumanData_27_ComparisonGEPsGarner/plots/mait_clusterslowres_hvg3712union_fastversion_fulltree.jpeg", width=1500, height=1500, res=250)
heatmap.2(mtn,
          # trace
          trace="none",
          # dendrogram
          # Rowv=FALSE,
          # Colv=FALSE,
          # dendrogram="none",
          # superimpose a density histogram on color key
          density.info="none",
          # color scale
          col=rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100)),
          breaks=seq(0,1,length=101),
          key.xlab = "AUROC",
          key.title="",
          keysize = 1.2,
          # text labels
          main="Union of HVGs (3,712 genes)",
          cexRow=0.6,
          cexCol=0.6,
          # margins
          margins=c(6,6))
# dev.off()





# Quick look at gene scores
seur.mait <- AddModuleScore(seur.mait, name=c("maitcommongenes_unstim", "maitcommongenes_cytok", "gep1_garner", "gep5_gapin"),
                            features=list("maitcommongenes_unstim"=genes_common_gep5_gep1garner,
                                          "maitcommongenes_cytok"=genes_common_gep5_gep3garner,
                                          "gep1_garner"=dflong[dflong$geneprogram=="GEP1_garner","gene"],
                                          "gep5_gapin"=dflong[dflong$geneprogram=="GEP5_gapin","gene"]),
                            seed=1)
colnames(seur.mait@meta.data)[42:45] <- c("maitcommongenes_unstim", "maitcommongenes_cytok", "gep1_garner", "gep5_gapin")
SCpubr::do_FeaturePlot(seur.mait,  features=c("maitcommongenes_unstim", "maitcommongenes_cytok"), order=F, ncol=2)
SCpubr::do_FeaturePlot(seur.mait,  features=c("gep1_garner", "gep5_gapin"), order=F, ncol=2)





