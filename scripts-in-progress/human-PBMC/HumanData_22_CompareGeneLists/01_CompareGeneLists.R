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
library(tidytext)
library(pheatmap)

# Import data
genes.gep     <- read.csv("./data/human-thymus/HumanData_17_GEPsOnParkData/genes_per_GEP_df_2023-04-07.csv", row.names=1)
genes.cano    <- readxl::read_excel("./data/human-thymus/HumanData_17_GEPsOnCanogamezData/canogamez_supp.xlsx", skip=3, sheet=1)
genesCD8.rose <- readxl::read_excel("./data/human-thymus/HumanData_22_CompareGeneLists/rose_supp_clusmodulesCD8.xlsx", sheet=1)
genesCD4.rose <- readxl::read_excel("./data/human-thymus/HumanData_22_CompareGeneLists/rose_supp_clusmodulesCD4.xlsx", sheet=1)
genes.poon    <- readxl::read_excel("./data/human-thymus/HumanData_22_CompareGeneLists/poon_supp.xlsx", sheet=6)[,-1]

# seurat object
seur.human <- readRDS("./data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS")
seur.pbmc <- subset(seur.human, Tissue=="PBMC")
print(seur.pbmc) # 17,204 genes and 41,238 cells

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
  mutate(geneprogram=sub("_", "", geneprogram)) # %>%
  # filter(geneprogram %in% c("GEP_1", "GEP_4", "GEP_5", "GEP_6", "GEP_7")) %>%
  # mutate(geneprogram=replace(geneprogram, geneprogram=="GEP1", "GEP1_MAIT"),
  #        geneprogram=replace(geneprogram, geneprogram=="GEP4", "GEP4_Temra"),
  #        geneprogram=replace(geneprogram, geneprogram=="GEP5", "GEP5_Tnaive"),
  #        geneprogram=replace(geneprogram, geneprogram=="GEP6", "GEP6_Tcm"))
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
  filter(p_val_adj < 0.05 & LFC > 0) %>% # lowest log2FC is 0.25
  # convert to long format
  select(cluster, gene) %>%
  rename(geneprogram=cluster)
head(canodf)
genes.long[["canogamez"]] <- canodf
# table(is.na(genes.long$canogamez$gene)) # sanity check
# table(genes.long$canogamez$gene %in% rownames(seur.human))

# Define colors
cols_cano <- c("CanoGamez | CD4_Tnaive"= "#b3e2cd",
               "CanoGamez | CD4_TCM"   = "#f4cae4",
               "CanoGamez | CD4_TEM"   = "#cbd5e8",
               "CanoGamez | CD4_TEMRA" = "#fdcdac",
               "CanoGamez | CD4_nTreg" = "#fbb4ae")
## end cano-gamez ####


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
# table(genes.long$rose$gene %in% rownames(seur.human))
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
## end rose ####


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
    # OPTION 1: full list of Poon genes
    filter(padj<0.01 & logFC>0.25)
    # OPTION 2: short list of Poon genes
    # slice_head(n=504) # take 500 top genes
  cat(geneprogram, "\n")
  cat("# genes:", nrow(subdf), "\n")
  cat("min log2FC:", min(subdf$logFC), "\n\n")
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
## end poon ####


# *******************************
## 2.5. Get it all in one df ####
longdf <- bind_rows(genes.long, .id="dataset")
table(longdf$dataset, useNA="ifany")
# canogamez     gapin      poon      rose 
#   387         19,407    31,852     1,672

# Remove any genes that we don't have in seur.human anyway
longdf <- longdf %>%
  filter(gene %in% rownames(seur.pbmc)) %>%
  filter(geneprogram != "GEP7")
table(longdf$dataset, useNA="ifany")
# canogamez     gapin      poon      rose 
#   341         15,378    22,246     1,542
# table(longdf$geneprogram)

# Create a colors vector & dataframe
cols_alldatasets <- c(cols_cano, cols_rose, cols_poon)
cols_df <- data.frame("geneprogram"=c(names(cols_cano), names(cols_rose), names(cols_poon)),
                      "geneprogram_cat"=c("Tnaive", "Tcm", "Tem", "Temra", "Treg", # canogamez
                                          "Tem", "Tcm", "Tem", "Tem", "Tnaive", # rose CD4
                                          "Temra", "Tcm", "Tem", "Tnaive", "Tnaive", # rose CD8
                                          "Trm", "Tnaive", "Tcm", "Trm", "Treg", "CD8 MAIT", "Tnaive", "Temra", "Trm"))
cols_df <- cols_df %>%
  mutate("color"=ifelse(geneprogram_cat=="Tnaive", "#b3e2cd",
                 ifelse(geneprogram_cat=="Tcm", "#f4cae4",
                 ifelse(geneprogram_cat=="Tem", "#cbd5e8",
                 ifelse(geneprogram_cat=="Temra", "#fdcdac",
                 ifelse(geneprogram_cat=="Treg", "#fbb4ae",
                 ifelse(geneprogram_cat=="Trm", "#4292c6", "#8c6bb1")))))))
cols_geneprogcat <- unique(cols_df$color)
names(cols_geneprogcat) <- unique(cols_df$geneprogram_cat)

## end bind ####




# ******************************************************
# 3. PROPORTION GEP GENES FOUND IN OTHER GENE LISTS ####
# ******************************************************
# library(ggupset)

# *********************
## 3.1. Upset plot ####
# # Gapin VS 1 other
# # jpeg("./data/human-thymus/HumanData_22_CompareGeneLists/gapin_poon_slicehead.jpeg", width=3000, height=1200, res=200)
# longdf %>%
#   filter(dataset %in% c("gapin", "poon")) %>% # restrict to only 2 datasets
#   group_by(gene) %>%
#   # keep only genes that are present in >1 dataset (not interested in intersections intra-dataset)
#   filter(n_distinct(dataset)>1) %>% 
#   summarize(geneprograms = list(geneprogram)) %>%
#   # filter(lengths(geneprograms)==2) %>%
#   ggplot(aes(x = geneprograms)) +
#     geom_bar() +
#     scale_x_upset() +
#     labs(title="Gapin vs Poon")
# # dev.off()
# 
# 
# # All datasets together
# gapingenes <- unique(genes.long$gapin$gene)
# longdf %>%
#   filter(gene %in% gapingenes) %>% # make sure it involves a Gapin GEP
#   group_by(gene) %>%
#   # keep only genes that are present in >1 dataset (not interested in intersections intra-dataset)
#   filter(n_distinct(dataset)>1) %>% 
#   summarize(geneprograms = list(geneprogram)) %>%
#   # filter(lengths(geneprograms)==2) %>%
#   ggplot(aes(x = geneprograms)) +
#   geom_bar() +
#   scale_x_upset(n_intersections=22, sets=unique(longdf$geneprogram)) +
#   labs(title="intersections with >20 genes in common")

## -end 3.1 ####


# *********************
## 3.2. Bar plot ####

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

# Get all GEP genes
# gapindf_all <- gather(genes.gep, key=geneprogram, value=gene, colnames(genes.gep)) %>%
#   filter(!is.na(gene)) %>%
#   mutate(geneprogram=sub("_", "", geneprogram)) # %>%
#   # mutate(geneprogram=replace(geneprogram, geneprogram=="GEP1", "GEP1_MAIT"),
#   #        geneprogram=replace(geneprogram, geneprogram=="GEP4", "GEP4_Temra"),
#   #        geneprogram=replace(geneprogram, geneprogram=="GEP5", "GEP5_Tnaive"),
#   #        geneprogram=replace(geneprogram, geneprogram=="GEP6", "GEP6_Tcm"))
# table(gapindf_all$geneprogram, useNA="ifany")
# table(is.na(gapindf_all$gene))

gapindf_all <- gapindf


# Cycle through the GEPs (plot all overlaps with our GEPs)
plist <- list()
for (gep in unique(gapindf_all$geneprogram)){
  # get list of genes for specific GEP
  print(gep)
  geneslist <- as.vector(gapindf_all %>% filter(geneprogram==gep) %>% pull(gene))
  print(length(geneslist))
  
  # Plot %of GEP genes found in other datasets
  # p <- longdf %>%
  #   as_tibble() %>%
  #   filter(dataset != "gapin") %>%
  #   filter(gene %in% geneslist) %>%
  #   group_by(geneprogram) %>%
  #   count() %>%
  #   ungroup() %>%
  #   mutate(total_gene=length(geneslist),
  #          prop_genes=n*100/total_gene) %>%
  #   ggplot(aes(x=reorder(geneprogram, -prop_genes), y=prop_genes, fill=geneprogram))+
  #     geom_bar(stat="identity")+
  #     labs(x="", y="% GEP genes found in each gene program", title=gep)+
  #     ylim(c(0,100))+
  #     scale_fill_manual(values=cols_alldatasets)+
  #     theme_cowplot()+
  #     theme(axis.text.x=element_text(angle=90, hjust=1), legend.position="none")
  
  # Plot overlap coefficient between GEP and other dataset gene set
  p <- longdf %>%
    as_tibble() %>%
    filter(dataset != "gapin") %>%
    group_by(geneprogram) %>%
    summarise(overlap=overlapcoef(a=geneslist, b=gene)) %>%
    ggplot(aes(x=reorder(geneprogram, -overlap), y=overlap, fill=geneprogram))+
      geom_bar(stat="identity")+
      labs(x="", y="overlap coef btw GEP/gene set", title=gep)+
      ylim(c(0,1))+
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
for (gep in unique(gapindf_all$geneprogram)){
  # get list of genes for specific GEP
  print(gep)
  geneslist <- as.vector(gapindf_all %>% filter(geneprogram==gep) %>% pull(gene))
  print(length(geneslist))
  
  df <- longdf %>%
    as_tibble() %>%
    filter(dataset != "gapin") %>%
    group_by(geneprogram) %>%
    summarise(overlap=overlapcoef(a=geneslist, b=gene)) %>%
    mutate(geneprogram1=gep) %>%
    rename(geneprogram2=geneprogram) %>%
    top_n(5, overlap)
  df.facet <- rbind(df.facet, df)
}

df.facet <- df.facet %>%
  mutate(geneprogram1=gsub("_.*", "", geneprogram1))

# Plot with facet
ggplot(df.facet, aes(x=reorder_within(geneprogram2, -overlap, geneprogram1), y=overlap, fill=geneprogram2))+
  geom_bar(stat="identity")+
  facet_wrap(~factor(geneprogram1, levels=paste0("GEP", 1:12)), nrow=1, scales="free_x")+
  tidytext::scale_x_reordered() +
  ylim(c(0,1))+
  scale_fill_manual(values=cols_alldatasets)+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1), legend.position="none")+
  labs(x="", y="% GEP genes found in each gene program")
# ggsave("./data/human-thymus/HumanData_22_CompareGeneLists/geneoverlapwithGEP_bars_top5.jpeg", width=15, height=8)

## end bar plot ####


# ******************
## 3.3. Heatmap ####

df.heatmp <- data.frame()

# Cycle through the gene sets (plot all jaccard index)
for (geneset in unique(longdf$geneprogram)){
  # get list of genes for specific geneset
  print(geneset)
  geneslist <- as.vector(longdf %>% filter(geneprogram==geneset) %>% pull(gene))
  print(length(geneslist))
  
  # Plot overlap coefficient between GEP and other dataset gene set
  df.heatmp <- rbind(df.heatmp,
                     longdf %>%
                       as_tibble() %>%
                       group_by(geneprogram) %>%
                       summarise(overlap=overlapcoef(a=geneslist, b=gene, coef="jaccardweight")) %>%
                       rename(geneprogram2=geneprogram) %>%
                       mutate(geneprogram1=geneset))
}

dim(df.heatmp)
nrow(df.heatmp)==35*35
table(df.heatmp$overlap==1) # should be 35 true

# Get a matrix
df.heatmp <- as.data.frame(pivot_wider(df.heatmp, names_from=geneprogram2, values_from=overlap))
rownames(df.heatmp) <- df.heatmp$geneprogram1
df.heatmp$geneprogram1 <- NULL
df.heatmp <- df.heatmp[,rownames(df.heatmp)]
df.heatmp <- as.matrix(df.heatmp)

# replace 1 with NAs
df.heatmp[df.heatmp==1] <- NA

# Plot heatmap
library(ComplexHeatmap)
Heatmap(df.heatmp,
        cluster_rows=F,
        cluster_columns = F,
        col=rev(colorRampPalette(RColorBrewer::brewer.pal(10, "Spectral"))(100)))

## end heatmp ####



# ****************************************
## 3.4. Bar plot with null hypothesis ####

# H0: any random list of 504 genes (length of GEP1) would have 60% overlap with Poon_CD8MAIT
# H1: GEP1 genes are specifically overlapping with Poon_CD8MAIT, with a 60% overlap (more than random)

# Define function
NullOverlap <- function(seuratobj=seur.pbmc,
                        df_geneprograms=longdf,
                        geneprograms_ref=unique(longdf[longdf$dataset=="gapin", "geneprogram"]), # by default GEPs
                        geneprograms_others=unique(longdf[longdf$dataset!="gapin", "geneprogram"]), # by default all others
                        coefficient_to_compute="overlapcoef",
                        nbins=10,
                        nrandom=1000){
  
  # ___________________________________________________________
  # -- 1. Get DF with all genes & binned by expression level --
  cat("\n-- 1. GET ALL GENES & BIN THEM BY EXPRESSION LEVEL --\n")
  allgenes_binned_DF <- data.frame("gene"=rownames(seuratobj),
                                   "totalexpression"=rowSums(seuratobj@assays$RNA@data))
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
  
  # ___________________________________________________________________________________________________________
  # -- 2. Draw random ctrl genes for each GEP & get observed overlap & see if significant compared to random --
  cat("\n-- 2. DRAW RANDOM CTRL GENES FOR EACH GEP AND GET PVAL --\n")
  # df.ref <- df_geneprograms[df_geneprograms$dataset=="gapin",]
  overlaps_per_gep_list <- list()
  
  # Loop through the "reference" gene programs (the ones we want to compare everything else against)
  for(refgeneprog in geneprograms_ref){
    cat("\n++ Drawing random ctrl genes for", refgeneprog, "\n")
    # Find in which expression bins are our genes from the gene set
    refgeneset <- df_geneprograms[df_geneprograms$geneprogram==refgeneprog, "gene"]
    refgeneset_binned_DF <- allgenes_binned_DF[allgenes_binned_DF$gene %in% refgeneset,]
    # table(nrow(refgeneset_binned_DF)==length(refgeneset))
    
    # ++++++++++++++++++++++++++++++++++
    # Do 1000 random draws of gene lists
    ctrlgenelist <- list()
    for(i in 1:nrandom){
      # Progress bar
      if(i%%100==0){cat(paste0("\nrandom draw #", i))}
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
      
      # Add ctrlgeneset to list
      ctrlgenelist[[i]] <- ctrlgeneset
    }
    # Sanity check
    cat("\n\n++ Selection of random ctrl genes over! ++\n")
    cat("Length of reference gene set is", length(unique(refgeneset)), "and length of the random gene sets:\n")
    print(table(sapply(ctrlgenelist, length)))
    # ++++++++++++++++++++++++++++++++++
    
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Get the % of gene overlap with every other gene program
    cat("\n\n++ Compute %overlap of random gene sets with other gene programs\n")
    overlap_null <- data.frame()
    # othergeneprograms <- unique(df_geneprograms[df_geneprograms$dataset!="gapin", "geneprogram"])
    othergeneprograms <- geneprograms_others
    # Sanity check
    # cat("There are [", length(othergeneprograms), "] gene programs from other datasets (should be 24)\n")
    for(geneprog in othergeneprograms){
      # cat("\n-> computing %random and %observed overlap with:", geneprog)
      # get genes from other gene program (CanoGamez, Rose or Poon)
      othergeneprogram_genes <- df_geneprograms %>% filter(geneprogram==geneprog) %>% pull(gene)
      # cat(length(othergeneprogram_genes)) # sanity check
      # get % gene overlap for each random draw with the gene program (% of genes from random gene set that can be found in gene program X)
      # ctrlpercent <- sapply(ctrlgenelist, function(x) sum(x %in% othergeneprogram_genes)*100/length(x))
      ctrlpercent <- sapply(ctrlgenelist, function(x) overlapcoef(x,othergeneprogram_genes, coef=coefficient_to_compute) ) # get overlap coefficient
      # get % gene overlap of GEP with the gene program (% of genes from GEP1 that can be found in gene program X)
      # observedpercent <- sum(refgeneset %in% othergeneprogram_genes)*100/length(refgeneset)
      observedpercent <- overlapcoef(refgeneset, othergeneprogram_genes, coef=coefficient_to_compute) # get overlap coefficient
      # print(observedpercent)
      df.temp <- data.frame("geneprogram"=geneprog, "randomoverlap"=ctrlpercent, "observedoverlap"=observedpercent)
      # add rows
      overlap_null <- rbind(overlap_null, df.temp)
    }
    # table(overlap_null$geneprogram) # should all be 1000 (for the 1000 draws)
    # table(overlap_null$observedoverlap)
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # Summarize overlap
    cat("\n\n++ Summarize %overlap of random/GEP gene sets with other gene programs\n")
    final_df <- overlap_null %>%
      as_tibble() %>%
      group_by(geneprogram) %>%
      summarize(observedoverlap=unique(observedoverlap),
                mean_randomoverlap=mean(randomoverlap),
                min_randomoverlap=min(randomoverlap),
                max_randomoverlap=max(randomoverlap),
                nb_random_above_observed=sum(randomoverlap>observedoverlap)) %>%
      mutate(pval=nb_random_above_observed/nrandom,
             padj=pval*length(unique(othergeneprograms)),
             refgenesetsize=length(unique(refgeneset))) %>%
      arrange(-observedoverlap)
    
    # sanity check
    cat("Gene program with highest overlap with", refgeneprog, "is [", pull(final_df[1,"geneprogram"]), "] with [", pull(final_df[1,"observedoverlap"]), "] overlap\n")
    
    # Add to list
    overlaps_per_gep_list[[refgeneprog]] <- final_df
    cat("\n ----------------------------------------------- \n")
    } # end of refgeneprog loop
  
  # Collapse list into one big dataframe <3
  bigdf <- bind_rows(overlaps_per_gep_list, .id="refgeneprog")
  
  return(bigdf)
}

NullOverlap_double <- function(seuratobj=seur.pbmc,
                               df_geneprograms=longdf,
                               geneprograms_to_test=unique(longdf$geneprogram), # by default GEPs
                               coefficient_to_compute="overlapcoef",
                               nbins=25,
                               nrandom=1000){
  
  # ___________________________________________________________
  # -- 1. Get DF with all genes & binned by expression level --
  cat("\n-- 1. GET ALL GENES & BIN THEM BY EXPRESSION LEVEL --\n")
  allgenes_binned_DF <- data.frame("gene"=rownames(seuratobj),
                                   "totalexpression"=rowSums(seuratobj@assays$RNA@data))
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
    # table(nrow(refgeneset_binned_DF)==length(refgeneset))
    
    # ++++++++++++++++++++++++++++++++++
    # Do 1000 random draws of gene lists
    ctrlgenematrix <- matrix(nrow=length(refgeneset), ncol=nrandom)
    for(i in 1:nrandom){
      if(i%%100==0){cat(paste0("\nrandom draw #", i))} # progress bar
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



# GET OVERLAP (RANDOM & OBSERVED)
# df_geneprograms should contain 3 columns: dataset (gapin, rose, etc.); geneprogram (GEP1, GEP2, ... GEP12); gene
genesets_overlap <- NullOverlap(seuratobj=seur.pbmc,
                                df_geneprograms = longdf,
                                geneprograms_ref=unique(longdf[longdf$dataset=="gapin","geneprogram"]),
                                geneprograms_others=unique(longdf[longdf$dataset=="gapin","geneprogram"]),
                                coefficient_to_compute="jaccard",
                                nbins=25, nrandom=1000)
# Sanity check
table(genesets_overlap$refgeneprog) # should all be 24 (for the 24 gene programs)

# Plot with facet
ggplot(genesets_overlap %>%
         left_join(cols_df, by="geneprogram") %>%
         mutate(padj_toplot=ifelse(padj==0, "< 0.001", ifelse(padj<0.05, as.character(round(padj, 2)), ""))),
       aes(x=reorder_within(geneprogram, -observedoverlap, refgeneprog), y=observedoverlap))+
  # geom_hline(yintercept = 10, color="lightgrey", linetype=2)+
  geom_bar(stat="identity", aes(fill=geneprogram_cat))+
  facet_wrap(~factor(refgeneprog, levels=paste0("GEP", 1:12)), ncol=3, scales="free_x")+
  geom_point(aes(y=mean_randomoverlap), shape="_", size=3)+
  tidytext::scale_x_reordered() +
  # ggrepel::geom_text_repel(aes(label=padj_toplot), size=4, direction="y", nudge_y=1, force=4) +
  geom_text(aes(label=padj_toplot), size=4, angle=90, hjust=-0.2)+
  ylim(c(0,1))+
  scale_fill_manual(values=cols_geneprogcat, name="")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10), #legend.position="none",
        panel.grid.major.y=element_line(colour="lightgrey", linetype=2),
        plot.margin=margin(10,10,10,70))+
  labs(x="", y="% GEP genes found in each gene program")
# ggsave("./data/human-thymus/HumanData_22_CompareGeneLists/geneoverlapwithGEP_bars2.jpeg", width=15, height=20)


# Plot with only GEPs of interest
genesets_overlap %>%
  left_join(cols_df, by="geneprogram") %>%
  group_by(refgeneprog) %>%
  top_n(10, observedoverlap) %>%
  ungroup() %>%
  mutate(padj_toplot=ifelse(padj==0, "< 0.001",
                            ifelse(padj<0.05, as.character(round(padj, 2)), ""))) %>%
  filter(refgeneprog %in% c("GEP1", "GEP4", "GEP5", "GEP6", "GEP11")) %>%
  # filter(padj<0.05) %>%
  distinct() %>%
ggplot(aes(x=reorder_within(geneprogram, -observedoverlap, refgeneprog), y=observedoverlap))+
  # geom_hline(yintercept = 10, color="lightgrey", linetype=2)+
  geom_bar(stat="identity", aes(fill=geneprogram_cat))+
  facet_wrap(~factor(refgeneprog, levels=paste0("GEP", 1:12)), nrow=1, scales="free_x")+
  geom_point(aes(y=mean_randomoverlap), shape="_", size=3)+
  tidytext::scale_x_reordered() +
  # ggrepel::geom_text_repel(aes(label=padj_toplot), size=4, direction="y", nudge_y=1, force=4) +
  geom_text(aes(label=padj_toplot), size=4, angle=90, hjust=-0.2)+
  ylim(c(0,1))+
  scale_fill_manual(values=cols_geneprogcat, name="")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=12), #legend.position="none",
        panel.grid.major.y=element_line(colour="lightgrey", linetype=2),
        plot.margin=margin(10,10,10,70))+
  # labs(x="", y="% GEP genes found in each gene program")
  labs(x="", y="Jaccard index")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_22_CompareGeneLists/plots/geneoverlapcoeff_bars_top10_refGapin_self_jaccard.jpeg", width=14, height=8)


# Plot with only GEPs of interest and Poon only
geps_of_interest <- c("GEP1", "GEP4", "GEP5", "GEP6", "GEP11", "GEP12")
df.facet2bis <- genesets_overlap %>%
  left_join(cols_df, by="geneprogram") %>%
  # keep only geps of interest and Poon
  filter(refgeneprog %in% geps_of_interest) %>%
  filter(geneprogram %in% grep("Poon", geneprogram, value=T)) %>%
  mutate(padj_toplot=ifelse(padj==0, "< 0.001",
                            ifelse(padj<0.05, as.character(round(padj, 2)), "")),
         refgeneprog=factor(refgeneprog, levels=geps_of_interest)) %>%
  distinct()
ggplot(df.facet2bis,
       aes(x=refgeneprog, y=observedoverlap))+
  # geom_hline(yintercept = 10, color="lightgrey", linetype=2)+
  geom_bar(stat="identity", aes(fill=geneprogram_cat))+
  facet_wrap(~geneprogram, nrow=2, scales="free_x")+
  geom_point(aes(y=mean_randomoverlap), shape="|", size=3)+
  tidytext::scale_x_reordered() +
  # geom_text(aes(label=padj_toplot), size=4, angle=90, hjust=-0.2)+
  geom_text(aes(label=padj_toplot), size=4, angle=0, hjust=-0.2)+
  ylim(c(0,100))+
  scale_fill_manual(values=cols_geneprogcat, name="")+
  coord_flip()+
  theme_cowplot()+
  # theme(axis.text.x=element_text(angle=45, hjust=1, size=12), #legend.position="none",
  #       panel.grid.major.y=element_line(colour="lightgrey", linetype=2),
  #       plot.margin=margin(10,10,10,70))+
  labs(x="", y="% GEP genes found in each Poon et al. gene program")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_22_CompareGeneLists/plots/geneoverlaptoGEP_bars_all_fullPoon_Poononly.jpeg", width=15, height=6)




# ****************************************************************
# 4. PROPORTION GENES FROM OTHER DATASETS FOUND IN OTHER GEPs ####
# ****************************************************************

# *********************************************
## 4.1. First look (without H0 generation) ####

# Cycle through the GEPs (get the top 5 overlaps for each gene program)
othergeneprograms <- unique(longdf[longdf$dataset!="gapin", "geneprogram"])
df.facet3 <- data.frame()
for (geneprog in othergeneprograms){
  # get list of genes for specific gene program
  print(geneprog)
  geneslist <- as.vector(longdf %>% filter(geneprogram==geneprog) %>% pull(gene))
  print(length(geneslist))
  
  df <- longdf %>%
    as_tibble() %>%
    filter(dataset == "gapin") %>%
    filter(gene %in% geneslist) %>%
    group_by(geneprogram) %>%
    count() %>%
    ungroup() %>%
    rename(geneprogram2=geneprogram) %>%
    mutate(total_gene=length(geneslist),
           prop_genes=n*100/total_gene,
           geneprogram1=geneprog) %>%
    top_n(5,prop_genes)
  df.facet3 <- rbind(df.facet3, df)
}

# Plot with facet (quick first look)
ggplot(df.facet3 %>% filter(geneprogram1 %in% unique(grep("CanoGamez|Rose", df.facet3$geneprogram1, value=T))) %>% data.frame(),
       aes(x=reorder_within(geneprogram2, -prop_genes, geneprogram1), y=prop_genes, fill=geneprogram2))+
  geom_bar(stat="identity")+
  facet_wrap(~geneprogram1, nrow=1, scales="free_x")+
  tidytext::scale_x_reordered() +
  ylim(c(0,100))+
  scale_fill_manual(values=cols_alldatasets)+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1), legend.position="none")+
  labs(x="", y="% genes found in each GEP")

## end first look ####


# ****************************************
## 4.2. Bar plot with null hypothesis ####

genesets_overlap2 <- NullOverlap(seuratobj=seur.pbmc,
                                 df_geneprograms = longdf,
                                 geneprograms_ref=unique(longdf[longdf$dataset!="gapin","geneprogram"]),
                                 geneprograms_others=unique(longdf[longdf$dataset=="gapin" & longdf$geneprogram!="GEP7", "geneprogram"]),
                                 nbins=25, nrandom=1000)

# Plot all with facet
geps_of_interest <- c("GEP1", "GEP4", "GEP5", "GEP6", "GEP11", "GEP12")
ggplot(genesets_overlap2 %>%
         left_join(cols_df, by=join_by("refgeneprog" == "geneprogram")) %>%
         mutate(padj_toplot=ifelse(padj==0, "< 0.001", ifelse(padj<0.05, as.character(round(padj, 2)), ""))) %>%
         filter(geneprogram %in% geps_of_interest) %>%
         # filter(refgeneprog %in% unique(grep("CanoGamez|Rose", genesets_overlap2$refgeneprog, value=T))) %>%
         mutate(geneprogram=factor(geneprogram, levels=geps_of_interest)),
       aes(x=geneprogram, # reorder_within(geneprogram, -observedoverlap, refgeneprog),
           y=observedoverlap))+
  # geom_hline(yintercept = 10, color="lightgrey", linetype=2)+
  geom_bar(stat="identity", aes(fill=geneprogram_cat))+
  facet_wrap(~refgeneprog, ncol=6, scales="free_x")+
  geom_point(aes(y=mean_randomoverlap), shape="|", size=3)+
  tidytext::scale_x_reordered() +
  # geom_text(aes(label=padj_toplot), size=4, angle=90, hjust=-0.2)+
  geom_text(aes(label=padj_toplot), size=4, angle=0, hjust=-0.2)+
  ylim(c(0,100))+
  scale_fill_manual(values=cols_geneprogcat, name="")+
  coord_flip()+
  theme_cowplot()+
  # theme(# axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10),
  #       # panel.grid.major.y=element_line(colour="lightgrey", linetype=2),
  #       # panel.grid.major.x=element_line(colour="lightgrey", linetype=2),
  #       plot.margin=margin(10,10,10,70))+
  labs(x="", y="% genes found in each GEP")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_22_CompareGeneLists/plots/geneoverlaptoGEP_bars_all_Poontop500.jpeg", width=17, height=9)


# Plot only gene programs that have at least 30% overlap with one of the GEPs
ggplot(genesets_overlap2 %>%
         left_join(cols_df, by=join_by("refgeneprog" == "geneprogram")) %>%
         mutate(padj_toplot=ifelse(padj==0, "< 0.001", ifelse(padj<0.05, as.character(round(padj, 2)), ""))) %>%
         filter(geneprogram %in% geps_of_interest) %>%
         # filter(refgeneprog %in% unique(grep("CanoGamez|Rose", genesets_overlap2$refgeneprog, value=T))) %>%
         mutate(geneprogram=factor(geneprogram, levels=geps_of_interest)) %>%
         group_by(refgeneprog) %>% mutate(max_observedoverlap=max(observedoverlap)) %>% ungroup() %>% filter(max_observedoverlap>33.3),
       aes(x=geneprogram, # reorder_within(geneprogram, -observedoverlap, refgeneprog),
           y=observedoverlap))+
  # geom_hline(yintercept = 10, color="lightgrey", linetype=2)+
  geom_bar(stat="identity", aes(fill=geneprogram_cat))+
  facet_wrap(~refgeneprog, ncol=6, scales="free_x")+
  geom_point(aes(y=mean_randomoverlap), shape="|", size=3)+
  tidytext::scale_x_reordered() +
  # geom_text(aes(label=padj_toplot), size=4, angle=90, hjust=-0.2)+
  geom_text(aes(label=padj_toplot), size=4, angle=0, hjust=-0.2)+
  ylim(c(0,100))+
  scale_fill_manual(values=cols_geneprogcat, name="")+
  coord_flip()+
  theme_cowplot()+
  # theme(# axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10),
  #       # panel.grid.major.y=element_line(colour="lightgrey", linetype=2),
  #       # panel.grid.major.x=element_line(colour="lightgrey", linetype=2),
  #       plot.margin=margin(10,10,10,70))+
  labs(x="", y="% genes found in each GEP")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_22_CompareGeneLists/plots/geneoverlaptoGEP_bars_all_Poontop500_short.jpeg", width=17, height=7)


## end 4/ ####


# ***************************************
# 5. OVERLAP COEFF OF OTHER DATASETS ####
# ***************************************

# *************************
## 5.1. Poon reference ####

# GET OVERLAP (RANDOM & OBSERVED)
# df_geneprograms should contain 3 columns: dataset (gapin, rose, etc.); geneprogram (GEP1, GEP2, ... GEP12); gene
overlap_poon <- NullOverlap(seuratobj=seur.pbmc,
                            df_geneprograms = longdf,
                            geneprograms_ref=unique(longdf[longdf$dataset=="poon","geneprogram"]),
                            geneprograms_others=unique(longdf[longdf$dataset!="poon", "geneprogram"]),
                            coefficient_to_compute="jaccard",
                            nbins=25, nrandom=1000)

# Sanity check
table(overlap_poon$refgeneprog) # should all be 27 (for the 27 gene programs)

# Plot
overlap_poon %>%
  left_join(cols_df, by="geneprogram") %>%
  group_by(refgeneprog) %>%
  top_n(10, observedoverlap) %>%
  ungroup() %>%
  mutate(padj_toplot=ifelse(padj==0, "< 0.001",
                            ifelse(padj<0.05, as.character(round(padj, 2)), ""))) %>%
  # filter(padj<0.05) %>%
  distinct() %>%
ggplot(aes(x=reorder_within(geneprogram, -observedoverlap, refgeneprog), y=observedoverlap))+
  # geom_hline(yintercept = 10, color="lightgrey", linetype=2)+
  geom_bar(stat="identity", aes(fill=geneprogram_cat))+
  facet_wrap(~refgeneprog, nrow=1, scales="free_x")+
  geom_point(aes(y=mean_randomoverlap), shape="_", size=3)+
  tidytext::scale_x_reordered() +
  # ggrepel::geom_text_repel(aes(label=padj_toplot), size=4, direction="y", nudge_y=1, force=4) +
  geom_text(aes(label=padj_toplot), size=4, angle=90, hjust=-0.2)+
  ylim(c(0,1))+
  scale_fill_manual(values=cols_geneprogcat, name="")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=12), #legend.position="none",
        panel.grid.major.y=element_line(colour="lightgrey", linetype=2),
        plot.margin=margin(10,10,10,70))+
  # labs(x="", y="% GEP genes found in each gene program")
  labs(x="", y="Jaccard index")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_22_CompareGeneLists/plots/geneoverlapcoeff_bars_top10_sig_refPoon_jaccard2.jpeg", width=20, height=8)


# *************************
## 5.2. Cano-Gamez reference ####

# GET OVERLAP (RANDOM & OBSERVED)
# df_geneprograms should contain 3 columns: dataset (gapin, rose, etc.); geneprogram (GEP1, GEP2, ... GEP12); gene
overlap_canogamez <- NullOverlap(seuratobj=seur.pbmc,
                                 df_geneprograms = longdf,
                                 geneprograms_ref=unique(longdf[longdf$dataset=="canogamez","geneprogram"]),
                                 geneprograms_others=unique(longdf[longdf$dataset!="canogamez", "geneprogram"]),
                                 coefficient_to_compute="jaccard",
                                 nbins=25, nrandom=1000)

# Sanity check
table(overlap_canogamez$refgeneprog) # should all be 31 (for the 31 other gene programs)

# Plot
overlap_canogamez %>%
  left_join(cols_df, by="geneprogram") %>%
  group_by(refgeneprog) %>%
  top_n(10, observedoverlap) %>%
  ungroup() %>%
  mutate(padj_toplot=ifelse(padj==0, "< 0.001",
                            ifelse(padj<0.05, as.character(round(padj, 2)), ""))) %>%
  filter(padj<0.05) %>%
  distinct() %>%
  ggplot(aes(x=reorder_within(geneprogram, -observedoverlap, refgeneprog), y=observedoverlap))+
  # geom_hline(yintercept = 10, color="lightgrey", linetype=2)+
  geom_bar(stat="identity", aes(fill=geneprogram_cat))+
  facet_wrap(~refgeneprog, nrow=1, scales="free_x")+
  geom_point(aes(y=mean_randomoverlap), shape="_", size=3)+
  tidytext::scale_x_reordered() +
  # ggrepel::geom_text_repel(aes(label=padj_toplot), size=4, direction="y", nudge_y=1, force=4) +
  geom_text(aes(label=padj_toplot), size=4, angle=90, hjust=-0.2)+
  ylim(c(0,1))+
  scale_fill_manual(values=cols_geneprogcat, name="")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=12), #legend.position="none",
        panel.grid.major.y=element_line(colour="lightgrey", linetype=2),
        plot.margin=margin(10,10,10,70))+
  # labs(x="", y="% GEP genes found in each gene program")
  labs(x="", y="Jaccard index")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_22_CompareGeneLists/plots/geneoverlapcoeff_bars_top10_sig_refCanoGamez_jaccard2.jpeg", width=14, height=8)



# *************************
## 5.3. Rose reference ####

# GET OVERLAP (RANDOM & OBSERVED)
# df_geneprograms should contain 3 columns: dataset (gapin, rose, etc.); geneprogram (GEP1, GEP2, ... GEP12); gene
overlap_rose <- NullOverlap(seuratobj=seur.pbmc,
                            df_geneprograms = longdf,
                            geneprograms_ref=unique(longdf[longdf$dataset=="rose","geneprogram"]),
                            geneprograms_others=unique(longdf[longdf$dataset!="rose", "geneprogram"]),
                            coefficient_to_compute="jaccard",
                            nbins=25, nrandom=1000)

# Sanity check
table(overlap_rose$refgeneprog) # should all be 26 (for the 26 other gene programs)

# Plot
overlap_rose %>%
  left_join(cols_df, by="geneprogram") %>%
  group_by(refgeneprog) %>%
  top_n(10, observedoverlap) %>%
  ungroup() %>%
  mutate(padj_toplot=ifelse(padj==0, "< 0.001",
                            ifelse(padj<0.05, as.character(round(padj, 2)), ""))) %>%
  # filter(padj<0.05) %>%
  distinct() %>%
  ggplot(aes(x=reorder_within(geneprogram, -observedoverlap, refgeneprog), y=observedoverlap))+
  # geom_hline(yintercept = 10, color="lightgrey", linetype=2)+
  geom_bar(stat="identity", aes(fill=geneprogram_cat))+
  facet_wrap(~refgeneprog, nrow=2, scales="free")+
  geom_point(aes(y=mean_randomoverlap), shape="_", size=3)+
  tidytext::scale_x_reordered() +
  # ggrepel::geom_text_repel(aes(label=padj_toplot), size=4, direction="y", nudge_y=1, force=4) +
  geom_text(aes(label=padj_toplot), size=4, angle=90, hjust=-0.2)+
  ylim(c(0,1))+
  scale_fill_manual(values=cols_geneprogcat, name="")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=12), #legend.position="none",
        panel.grid.major.y=element_line(colour="lightgrey", linetype=2),
        plot.margin=margin(10,10,10,70))+
  # labs(x="", y="% GEP genes found in each gene program")
  labs(x="", y="Jaccard")
# ggsave("./scripts-in-progress/human-PBMC/HumanData_22_CompareGeneLists/plots/geneoverlapcoeff_bars_top10_sig_refRose_jaccard2.jpeg", width=16, height=16)

