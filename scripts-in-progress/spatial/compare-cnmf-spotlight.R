#################
## libraries ####
#################
library(Seurat)
library(ggplot2)
library(dplyr)
library(tblhelpr)
library(NMF)
library(ComplexHeatmap)
library(corrplot)


############
## data ####
############
colvalues <- colorRampPalette(brewer.pal(12,"Paired"))(18)
#  'paletteMartin' from the 'colorBlindness' package)
paletteMartin <- c("#000000", "#004949", "#009292", "#ff6db6", "#ffb6db",
                            "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
                            "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")
                            
colhisto <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d')


## learned topics from SPOTlight (integrate-clusers-spatially.R)
spotlight <- readRDS("results/NMF.modls.thymus.batchC.allslices.rds")


## learned profiles from cnmf (Laurent)
cnmf <- readr::read_delim("results/All_cells_cNMF.spectra.k_19.dt_0_1.consensus.txt")
colnames(cnmf)[1] <- "topic"

################
## analysis ####
################

## generate consensus topics for spotlight (simple median?) ####
topics <- lapply(seq_along(spotlight), function(x) {
  tmp <- basis(spotlight[[x]]$NMF)
  colnames(tmp) <- paste0("Topic", seq_len(ncol(tmp)))
  tmp %>%
    as.data.frame %>%
    rownames_to_column("gene") %>%
    mutate(slice = names(spotlight)[x])
}) %>% 
  bind_rows %>%
  as_tibble

topics_consensus <- topics %>%
  group_by(gene) %>%
  summarise(across(starts_with("Topic"), median))
colnames(topics_consensus) <- c("gene",
                                paste0("topic_spotlight", 1:(ncol(topics_consensus)-1)))

## format cnmf
cnmf <- cnmf %>%
  mutate(topic=paste0("topic_cnmf", 1:n())) %>%
  transpose_tibble(col_names = "topic", id_col = "gene")

## filter and scale
spotlight_common <- topics_consensus %>%
  filter(gene %in% cnmf$gene) %>%
  arrange(gene)

cnmf_common <- cnmf %>%
  filter(gene %in% spotlight_common$gene) %>%
  arrange(gene)

common <- left_join(spotlight_common, cnmf_common) %>%
  as.data.frame %>%
  column_to_rownames("gene") %>%
  data.matrix

spotlight_scale <- spotlight_common %>%
  pivot_longer(starts_with("topic"), values_to="value", names_to="topic") %>%
  group_by(gene) %>%
  mutate(value=(value-mean(value)/sd(value))) %>%
  pivot_wider(names_from="topic", values_from = "value")

cnmf_scale <- cnmf_common %>%
  pivot_longer(starts_with("topic"), values_to="value", names_to="topic") %>%
  group_by(gene) %>%
  mutate(value=(value-mean(value)/sd(value))) %>%
  pivot_wider(names_from="topic", values_from = "value")


common_scale <- left_join(spotlight_scale, cnmf_scale) %>%
  as.data.frame %>%
  column_to_rownames("gene") %>%
  data.matrix

Heatmap(common_scale,
        cluster_columns = FALSE,
        show_row_names = FALSE)

cor_common <- cor(common)

corrplot(cor_common,
         order="AOE")

corc <- cor_common[grepl("cnmf", rownames(cor_common)),
                   grepl("spotlight", colnames(cor_common))]
corrplot(corc)
