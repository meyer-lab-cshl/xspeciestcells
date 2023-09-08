## functions ####
getGEPgenes <- function(i, geps, threshold, plot=FALSE) {   
  thresh <- threshold[i]
  gepDfs <- as.data.frame(sort(geps[,i], decreasing = T))
  gepDfs$rank <- seq(1, nrow(gepDfs))
  names(gepDfs)[1] <- "signal"
  
  columns = c("rank", "slope")
  ranks = data.frame(matrix(nrow = 0, ncol = length(columns))) 
  names(ranks) <- columns
  
  for (j in seq(2:nrow(gepDfs))){
    curr = predict(fit, j, deriv = 0)
    prev = predict(fit, j-1, deriv = 0)
    ranks[j-1, "rank"] = paste(j, j-1, sep = "_")
    ranks[j-1,"slope"] = prev$y - curr$y
  }
  rankSelect <- min(which((ranks$slope > thresh) == FALSE))
  newThresh <- predict(fit, rankSelect, deriv = 0)
  if (plot) {
    plot(gepDfs$rank,
         gepDfs$signal, 
         xlab = "rank",
         ylab = paste("gep", i, "signal", sep = " "))
    fit <- smooth.spline(gepDfs$rank, gepDfs$signal)
    lines(fit, col = "red", lwd = 2)
    abline(v = rankSelect, col = "gray60", lty = 2)
    abline(h = newThresh, col = "gray60", lty = 2, )
    mtext(paste("thresh = ",
                formatC(round(newThresh$y, digits = 5), format = "e"),
                "\n rank = ", rankSelect, sep =""), side=3)
  }
  return(rankSelect)
}




## data ####
geps_nbc <- read.delim("data/cnmf/NoBatchCorrection/NoBatchCorrection.gene_spectra_score.k_12.dt_0_2.txt",
                   row.names = 1)
geps_nbc <- t(geps_nbc)

geps_bc <- read.delim("data/cnmf/BatchCorrection/BatchCorrection.gene_spectra_score.k_12.dt_0_2.txt",
                   row.names = 1)
geps_bc <- t(geps_bc)

seur <- readRDS("data/raw_data/human_data/seurat_filtered_harmony_08_28_23.RDS")

## analysis ####
threshold <- rep(1.101706e-06, 12)

#geps_bc_sig <- lapply(1:ncol(geps_bc), getGEPgenes, geps_bc, threshold)
#geps_nbc_sig <- lapply(1:ncol(geps_nbc), getGEPgenes, geps_nbc, threshold)


## non-batch corrected ####
plot=FALSE
gepDfs <- c()
rankSelect <- c()
for (i in seq(1:ncol(geps_nbc))){
  gepname <- paste("gep", i, sep = "_")
  gepDfs[[gepname]] <- as.data.frame(sort(geps_nbc[,i], decreasing = T))
  gepDfs[[gepname]]$rank <- seq(1, nrow(gepDfs[[gepname]]))
  names(gepDfs[[gepname]])[1] <- "signal"
  
  
  columns = c("rank","slope")
  ranks = data.frame(matrix(nrow = 0, ncol = length(columns))) 
  names(ranks) <- columns
  for (j in seq(2:nrow(gepDfs[[gepname]]))){
    curr = predict(fit, j, deriv = 0)
    prev = predict(fit, j-1, deriv = 0)
    ranks[j-1, "rank"] = paste(j, j-1, sep = "_")
    ranks[j-1,"slope"] = prev$y - curr$y
  }
  rankSelect[[i]] <- min(which((ranks$slope > thresh) == FALSE))
  newThresh <- predict(fit, rankSelect[[i]], deriv = 0)
  
  if (plot) {
    plot(gepDfs[[gepname]]$rank,
         gepDfs[[gepname]]$signal,
         xlab = "rank",
         ylab = paste("gep", i, "signal", sep = " "))
    fit <- smooth.spline(gepDfs[[gepname]]$rank,gepDfs[[gepname]]$signal)
    lines(fit, col = "red", lwd = 2)
    abline(v = rankSelect[[i]], col = "gray60", lty = 2)
    abline(h = newThresh, col = "gray60", lty = 2)
    mtext(paste("thresh = ", formatC(round(newThresh$y, digits = 5),
                                     format = "e"), "\n rank = ", 
                rankSelect[[i]], sep =""), side=3)
  }
}

gene_cutoffs_nbc = data.frame(rank = seq(1:max(unlist(rankSelect))))
for (g in seq(1:12)){
  geps_sort <- as.data.frame(sort(geps_nbc[,g], decreasing = T))
  top_genes <- as.data.frame(rownames(geps_sort)[1:rankSelect[[g]]])
  colnames(top_genes) <- paste("GEP", g, sep = "")
  top_genes$rank <- seq(1:rankSelect[[g]])
  gene_cutoffs_nbc = merge(gene_cutoffs_nbc, top_genes, by = "rank", all = T)
}
gene_cutoffs_list_nbc <- lapply(as.list(gene_cutoffs_nbc[,-1]), 
                                function(x) x[!is.na(x)])


## batch corrected ####
plot=FALSE
gepDfs <- c()
rankSelect <- c()
for (i in seq(1:ncol(geps_bc))){
  gepname <- paste("gep", i, sep = "_")
  gepDfs[[gepname]] <- as.data.frame(sort(geps_bc[,i], decreasing = T))
  gepDfs[[gepname]]$rank <- seq(1, nrow(gepDfs[[gepname]]))
  names(gepDfs[[gepname]])[1] <- "signal"
  
  
  columns = c("rank","slope")
  ranks = data.frame(matrix(nrow = 0, ncol = length(columns))) 
  names(ranks) <- columns
  for (j in seq(2:nrow(gepDfs[[gepname]]))){
    curr = predict(fit, j, deriv = 0)
    prev = predict(fit, j-1, deriv = 0)
    ranks[j-1, "rank"] = paste(j, j-1, sep = "_")
    ranks[j-1,"slope"] = prev$y - curr$y
  }
  rankSelect[[i]] <- min(which((ranks$slope > thresh) == FALSE))
  newThresh <- predict(fit, rankSelect[[i]], deriv = 0)
  
  if (plot) {
    plot(gepDfs[[gepname]]$rank,
         gepDfs[[gepname]]$signal,
         xlab = "rank",
         ylab = paste("gep", i, "signal", sep = " "))
    fit <- smooth.spline(gepDfs[[gepname]]$rank,gepDfs[[gepname]]$signal)
    lines(fit, col = "red", lwd = 2)
    abline(v = rankSelect[[i]], col = "gray60", lty = 2)
    abline(h = newThresh, col = "gray60", lty = 2)
    mtext(paste("thresh = ", formatC(round(newThresh$y, digits = 5),
                                     format = "e"), "\n rank = ", 
                rankSelect[[i]], sep =""), side=3)
  }
}


gene_cutoffs_bc = data.frame(rank = seq(1:max(unlist(geps_bc_sig))))
for (g in seq(1:12)){
  geps_sort <-   as.data.frame(sort(geps_bc[,g], decreasing = T))
  top_genes <- as.data.frame(rownames(geps_sort)[1:rankSelect[[g]]])
  colnames(top_genes) <- paste("GEP", g, sep = "")
  top_genes$rank <- seq(1:rankSelect[[g]])
  gene_cutoffs_bc = merge(gene_cutoffs_bc, top_genes, by = "rank", all = T)
}
gene_cutoffs_list_bc <- lapply(as.list(gene_cutoffs_bc[,-1]),
                               function(x) x[!is.na(x)])

## Visualise ####

seur <- AddModuleScore(seur, name = "GEP_bc",
                            features=gene_cutoffs_list_bc)
seur <- AddModuleScore(seur, name = "GEP_nbc",
                            features=gene_cutoffs_list_nbc)


SCpubr::do_FeaturePlot(seur, reduction="UMAP_50",
                       features=paste0("GEP_bc", 1:12),
                       ncol=4, order=TRUE,
                       viridis_color_map = "B", raster=TRUE)
ggsave("data/cnmf/BatchCorrection/cNMF_geps_batchcorrected.pdf", width=30,
       height=20)

SCpubr::do_FeaturePlot(seur, reduction="UMAP_50",
                       features=paste0("GEP_nbc", 1:12),
                       ncol=4, order=TRUE,
                       viridis_color_map = "B", raster=TRUE)
ggsave("data/cnmf/NoBatchCorrection/cNMF_geps_nonbatchcorrected.pdf", width=30,
       height=20)

