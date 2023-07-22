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
    slice_head(n=500) # take 500 top genes
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
               "Poon | CyclingTRM"="#fff2ae")


# *******************************
## 2.5. Get it all in one df ####
longdf <- bind_rows(genes.long, .id="dataset")
table(longdf$dataset, useNA="ifany")

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
gapindf_all <- gather(genes.gep, key=geneprogram, value=gene, colnames(genes.gep)) %>%
  filter(!is.na(gene)) %>%
  mutate(geneprogram=sub("_", "", geneprogram)) # %>%
  # mutate(geneprogram=replace(geneprogram, geneprogram=="GEP1", "GEP1_MAIT"),
  #        geneprogram=replace(geneprogram, geneprogram=="GEP4", "GEP4_Temra"),
  #        geneprogram=replace(geneprogram, geneprogram=="GEP5", "GEP5_Tnaive"),
  #        geneprogram=replace(geneprogram, geneprogram=="GEP6", "GEP6_Tcm"))
table(gapindf_all$geneprogram, useNA="ifany")
table(is.na(gapindf_all$gene))


# Cycle through the GEPs (plot all overlaps with our GEPs)
plist <- list()
for (gep in unique(gapindf_all$geneprogram)){
  print(gep)
  geneslist <- as.vector(gapindf_all %>% filter(geneprogram==gep) %>% pull(gene))
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
for (gep in unique(gapindf_all$geneprogram)){
  print(gep)
  geneslist <- as.vector(gapindf_all %>% filter(geneprogram==gep) %>% pull(gene))
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
  # ylim(c(0,30))+
  scale_fill_manual(values=cols_alldatasets)+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1), legend.position="none")+
  labs(x="", y="% GEP genes found in each gene program")
# ggsave("./data/human-thymus/HumanData_22_CompareGeneLists/geneoverlapwithGEP_bars_top5.jpeg", width=15, height=8)


# ****************************************
## 3.3. Bar plot with null hypothesis ####

# H0: any random list of 504 genes (length of GEP1) would have 20% overlap with Poon_CD8MAIT
# H1: GEP1 genes are specifically overlapping with Poon_CD8MAIT, with a 20% overlap (more than random)

# Define function
NullOverlap <- function(seuratobj=seur.pbmc, df_geneprograms=longdf, geps, nbins=10, nrandom=1000){
  
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
  gapindf_all <- df_geneprograms[df_geneprograms$dataset=="gapin",]
  overlaps_per_gep_list <- list()
  
  for(gep in geps){
    cat("\n++ Drawing random ctrl genes for", gep, "\n")
    # Find in which expression bins are our genes from the gene set
    gepgeneset <- gapindf_all[gapindf_all$geneprogram==gep, "gene"]
    gepgeneset_binned_DF <- allgenes_binned_DF[allgenes_binned_DF$gene %in% gepgeneset,]
    # table(nrow(gepgeneset_binned_DF)==length(gepgeneset))
    
    # ++++++++++++++++++++++++++++++++++
    # Do 1000 random draws of gene lists
    ctrlgenelist <- list()
    for(i in 1:nrandom){
      # Progress bar
      if(i%%100==0){cat(paste0("\nrandom draw #", i))}
      # Loop by expression bin & sample equal number of genes from each bin
      ctrlgeneset <- c()
      for(bin in unique(gepgeneset_binned_DF$bin)){
        # print(bin)
        ngenes_to_sample <- sum(gepgeneset_binned_DF$bin==bin)
        # print(ngenes_to_sample)
        ctrlgenes_in_same_bin <- sample(x=allgenes_binned_DF[allgenes_binned_DF$bin==bin,"gene"], size=ngenes_to_sample, replace=FALSE)
        # print(length(ctrlgenes_in_same_bin))
        ctrlgeneset <- c(ctrlgeneset, ctrlgenes_in_same_bin)
      }
      # Sanity check
      if(length(ctrlgeneset)!=length(gepgeneset)){cat("\nPROBLEM: control gene set of size", length(ctrlgeneset), "while GEP length is", length(gepgeneset))}
      
      # Add ctrlgeneset to list
      ctrlgenelist[[i]] <- ctrlgeneset
    }
    # Sanity check
    cat("\n\n++ Selection of random ctrl genes over! ++\n")
    cat("Length of GEP is", length(unique(gepgeneset)), "and length of the random gene sets:\n")
    print(table(sapply(ctrlgenelist, length)))
    # ++++++++++++++++++++++++++++++++++
    
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Get the % of gene overlap with every other gene program
    cat("\n\n++ Compute %overlap of random gene sets with other gene programs\n")
    overlap_null <- data.frame()
    othergeneprograms <- unique(df_geneprograms[df_geneprograms$dataset!="gapin", "geneprogram"])
    # Sanity check
    # cat("There are [", length(othergeneprograms), "] gene programs from other datasets (should be 24)\n")
    for(geneprog in othergeneprograms){
      # cat("\n-> computing %random and %observed overlap with:", geneprog)
      # get genes from gene program (CanoGamez, Rose or Poon)
      geneprogram_genes <- df_geneprograms %>% filter(geneprogram==geneprog) %>% pull(gene)
      # cat(length(geneprogram_genes)) # sanity check
      # get % gene overlap for each random draw with the gene program (% of genes from random gene set that can be found in gene program X)
      ctrlpercent <- sapply(ctrlgenelist, function(x) sum(x %in% geneprogram_genes)*100/length(x))
      # get % gene overlap of GEP with the gene program (% of genes from GEP1 that can be found in gene program X)
      observedpercent <- sum(gepgeneset %in% geneprogram_genes)*100/length(gepgeneset)
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
             gepsize=length(unique(gepgeneset))) %>%
      arrange(-observedoverlap)
    
    # sanity check
    cat("Gene program with highest overlap with", gep, "is [", pull(final_df[1,"geneprogram"]), "] with [", pull(final_df[1,"observedoverlap"]), "% ] overlap\n")
    
    # Add to list
    overlaps_per_gep_list[[gep]] <- final_df
    cat("\n ----------------------------------------------- \n")
    } # end of GEP loop
  
  # Collapse list into one big dataframe <3
  bigdf <- bind_rows(overlaps_per_gep_list, .id="gep_program")
  
  return(bigdf)
}

# GET OVERLAP (RANDOM & OBSERVED)
# df_geneprograms should contain 3 columns: dataset (gapin, rose, etc.); geneprogram (GEP1, GEP2, ... GEP12); gene
seur.pbmc <- subset(seur.human, Tissue=="PBMC")
print(seur.pbmc) # 17,204 genes and 41,238 cells

genesets_overlap <- NullOverlap(seuratobj=seur.pbmc,
                                df_geneprograms = longdf,
                                geps=unique(longdf[longdf$dataset=="gapin","geneprogram"]),
                                nbins=25, nrandom=1000)
# Sanity check
table(genesets_overlap$gep_program) # should all be 24 (for the 24 gene programs)

# df.facet2 <- genesets_overlap %>%
#   left_join(cols_df, by="geneprogram") %>%
#   group_by(gep_program) %>%
#   top_n(5, observedoverlap) %>%
#   ungroup() %>%
#   mutate(padj_toplot=ifelse(padj==0, "< 0.001",
#                             ifelse(padj<0.05, as.character(round(padj, 2)), "")))

# Plot with facet
ggplot(genesets_overlap %>%
         left_join(cols_df, by="geneprogram") %>%
         mutate(padj_toplot=ifelse(padj==0, "< 0.001", ifelse(padj<0.05, as.character(round(padj, 2)), ""))),
       aes(x=reorder_within(geneprogram, -observedoverlap, gep_program), y=observedoverlap))+
  # geom_hline(yintercept = 10, color="lightgrey", linetype=2)+
  geom_bar(stat="identity", aes(fill=geneprogram_cat))+
  facet_wrap(~factor(gep_program, levels=paste0("GEP", 1:12)), ncol=3, scales="free_x")+
  geom_point(aes(y=mean_randomoverlap), shape="_", size=3)+
  tidytext::scale_x_reordered() +
  # ggrepel::geom_text_repel(aes(label=padj_toplot), size=4, direction="y", nudge_y=1, force=4) +
  geom_text(aes(label=padj_toplot), size=4, angle=90, hjust=-0.2)+
  ylim(c(0,50))+
  scale_fill_manual(values=cols_geneprogcat, name="")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10), #legend.position="none",
        panel.grid.major.y=element_line(colour="lightgrey", linetype=2),
        plot.margin=margin(10,10,10,70))+
  labs(x="", y="% GEP genes found in each gene program")
# ggsave("./data/human-thymus/HumanData_22_CompareGeneLists/geneoverlapwithGEP_bars2.jpeg", width=15, height=20)


# Plot with only GEPs of interest
df.facet2 <- genesets_overlap %>%
  left_join(cols_df, by="geneprogram") %>%
  group_by(gep_program) %>%
  top_n(10, observedoverlap) %>%
  ungroup() %>%
  mutate(padj_toplot=ifelse(padj==0, "< 0.001",
                            ifelse(padj<0.05, as.character(round(padj, 2)), "")))
ggplot(df.facet2 %>% filter(gep_program %in% c("GEP1", "GEP4", "GEP5", "GEP6", "GEP11")),
       aes(x=reorder_within(geneprogram, -observedoverlap, gep_program), y=observedoverlap))+
  # geom_hline(yintercept = 10, color="lightgrey", linetype=2)+
  geom_bar(stat="identity", aes(fill=geneprogram_cat))+
  facet_wrap(~factor(gep_program, levels=paste0("GEP", 1:12)), nrow=1, scales="free_x")+
  geom_point(aes(y=mean_randomoverlap), shape="_", size=3)+
  tidytext::scale_x_reordered() +
  # ggrepel::geom_text_repel(aes(label=padj_toplot), size=4, direction="y", nudge_y=1, force=4) +
  geom_text(aes(label=padj_toplot), size=4, angle=90, hjust=-0.2)+
  ylim(c(0,50))+
  scale_fill_manual(values=cols_geneprogcat, name="")+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=12), #legend.position="none",
        panel.grid.major.y=element_line(colour="lightgrey", linetype=2),
        plot.margin=margin(10,10,10,70))+
  labs(x="", y="% GEP genes found in each gene program")
ggsave("./data/human-thymus/HumanData_22_CompareGeneLists/pdf_plots/geneoverlapwithGEP_bars_top10_pval.jpeg", width=15, height=8)





