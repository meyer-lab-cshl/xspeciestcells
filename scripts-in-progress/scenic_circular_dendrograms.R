#################
## libraries ####
#################

library(readr)
library(tidyverse)
library(dplyr)
library(readxl)
library(igraph)
library(ggraph)
library(rgexf)
library(RColorBrewer)
library(scales)
library(ggnewscale)

#x11(type="cairo")
#################
## functions ####
#################


scenic2network <- function(weights, attributes, needsRoot=FALSE) {
  weights <- replace(weights, is.na(weights), values = 0) # replace NAs with 0
  weights[weights < 3] <- 0
  weights$dummy <- 10
  
  names(attributes)[1] <- "regulon"
  names(attributes)[2] <- "target"
  names(attributes)[3] <- "attribute"
  all_regulons <- unique(attributes$regulon)
  all_targets <- unique(attributes$target)
  
  # dataframe to store all regulon-target pairwise comparisons
  pairwise_source_target <- data.frame(matrix(ncol = 3, nrow = 0)) 
  colnames(pairwise_source_target) <- c("regulon", "target", "weight")
  
  # generate all regulon-target pairwise comparisons
  for (i in all_regulons){
    for (j in all_targets){
      tmp <- data.frame(regulon = i, target = j, weight = weights[j,i])
      pairwise_source_target <- rbind(pairwise_source_target, tmp)
    }
  }
  
  # subset df to regulon-target pairs != 0 as their reproducible weight
  non_zero_weights <-
    pairwise_source_target[which(pairwise_source_target$weight != 0),]
  
  if (needsRoot) {
    df_weights <- tibble::tibble(regulon='dummy', 
                                 target=unique(non_zero_weights$regulon),
                                 weight=10)
    non_zero_weights <- rbind(non_zero_weights, df_weights)
  }
  
  # add metadata
  subset_attr <- attributes[,c("target", "attribute")]
  df_attr <- tibble::tibble(target='dummy', attribute="Other")
  subset_att <- rbind(subset_attr, df_attr)
  non_zero_weights_with_attrs <- merge(non_zero_weights, subset_attr,
                                       by = "target", all.x = TRUE)
  
  
  # a dataframe where each row is a node in the network and the value is the 
  # vertex number (helps with reordering and aesthetics)
  target_to_att_index_map <- unique(non_zero_weights_with_attrs[,c('target', 
                                                                   'attribute')])
  missing_nodes <- unique(non_zero_weights_with_attrs$regulon)[
    which(unique(non_zero_weights_with_attrs$regulon) %in% 
            target_to_att_index_map$target == FALSE)]
  missing_nodes_df <- data.frame(target = missing_nodes,
                                 attribute = rep("TF", length(missing_nodes)))
  target_to_att_index_map <- rbind(target_to_att_index_map, missing_nodes_df)
  
  # order by attribute
  target_to_att_index_map <-
    target_to_att_index_map[order(target_to_att_index_map$attribute,
                                  decreasing = TRUE),]
  
  # build the network graph from the data frame
  network <- graph_from_data_frame(d = non_zero_weights, directed = TRUE)
  rownames(target_to_att_index_map) <- target_to_att_index_map$target
  target_to_att_index_map <- target_to_att_index_map[names(as.list(V(network))),]
  network <- set_vertex_attr(network, name = "group", index = V(network),
                             target_to_att_index_map$attribute)
}


plotCircDend <- function(network, labelsize=2) {
  # plots circular dendrogram
  lay <- create_layout(network, layout = 'dendrogram', circular = TRUE)
  ## internal nodes (inside circle radius of 1, rounding errors, so use smaller r)
  internal <- sqrt(lay$x^2 + lay$y^2) < 0.8
  leaves <- sqrt(lay$x^2 + lay$y^2) > 0.8
  
  node_angles <- atan((lay$y[leaves]/lay$x[leaves]))*180/pi
  node_angles[is.na(node_angles)] <- 0
  
  lay$angle[internal] <- 0
  lay$angle[leaves] <- node_angles
  
  lay$hjust <- ifelse(lay$x < 0 , 1.2, - 0.2)
  
  ggraph(lay) +
    geom_edge_diagonal() +
    geom_node_point(aes(filter = leaf, color = group)) +
    scale_color_manual(values = colors_list, guide="none") +
    labs(color="Type") +
    new_scale_colour() +
    geom_node_label(aes(label = ifelse(leaf == FALSE, name, NA), color=group),
                    repel = TRUE, size=labelsize) +
    geom_node_text(aes(label = ifelse(leaf == FALSE, NA, name),
                       color=group), 
                   angle = lay$angle, hjust=lay$hjust, size=labelsize) +
    scale_color_manual(values = colors_list, guide="none") +
    coord_fixed(clip = "off") +
    theme(panel.background = element_blank(),
          legend.key =  element_rect(fill="white"),
          legend.box.background = element_blank(),
          plot.margin = unit(rep(40, 4), "points"))
}


#############
## setup ####
#############

colors <- tibble::tibble(attribute=c("TF", "Cytokine", "Chemokine", "Other",
                                     "Cytokine Receptor", "Cytotoxicity",
                                     "Chemokine Receptor", "Integrin",
                                     "NK receptor"),
                         color=c("#FFD382", "#964A35", "#964A35", "gray",
                                 "#8EB63B", "#C7B354", "#001959", "#45204C",
                                 "#FFACAC"))

colors_list <- as.list(colors$colors)
names(colors_list) <- colors$attribute


##########################
## INPUTS Description ####
##########################
# * weights *
# this the reproducible counts summarized across all seed runs for every regulon
# and every reproducible target gene from each SCENIC run; to
# obtain this csv, you will need to run: 
# STEP1_reproducible_regulons_and_gene_targets_across_runs.py (Tonya can provide 
# if needed)
# column names = regulon; row names = target genes

# * attributes *
# 3 column csv file:
#   - column 1 = name of regulon (column name is "regulon")
#   - column 2 = name of target gene (column name is "target")
#   - column 3 = name of attribute/characteristic/class (column name is 
#   "attribute")
# ex: column 1 tells you the regulon name, column 2 tells you which target gene 
# the regulon is acting upon (arrow coming out of regulon node pointing to 
# target gene node),
# column 3 tells you the attribute or class; if I listed CEBPD in column 1, 
# that means it is a regulon, if in column 2, I list
# GZMK, this means CEBPD regulates the target gene GZMK (so arrow coming from 
# CEBPD node and pointing to GZMK node), if I list "Cytotoxicity"in column 3, 
# that means I am give GZMK the class name "Cytotoxicity"
# (colors the node and network arrow line all the same color for this class)
 
################
## Analysis ####
################

## GEP 1
weights_gep1 <- read.csv("data/scenic_gep_data_with_threshold_50target_70regulon/GEP1_seurat_filtered_harmony_02_15_23_raw_counts_multirun_04192023_regulons_pruned_target_gene_counts.csv",
                         row.names = 1)
attributes_gep1 <- read.csv("data/Scenic_gep_data_with_threhold_30target_70regulon/GEP_1_pruned.csv", row.names = 1)

network_gep1 <- scenic2network(weights_gep1, attributes_gep1)
p_gep1 <- plotCircDend(network_gep1)
ggsave(plot=p_gep1, "results/scenic_gep1.pdf", width=5.5, height=5.5)

## GEP 4
weights_gep4 <- read.csv("data/scenic_gep_data_with_threshold_50target_70regulon/GEP4_seurat_filtered_harmony_02_15_23_raw_counts_multirun_04192023_regulons_pruned_target_gene_counts.csv",
                    row.names = 1)
attributes_gep4 <- read.csv("data/Scenic_gep_data_with_threhold_30target_70regulon/GEP_4_pruned.csv", row.names = 1)

network_gep4 <- scenic2network(weights_gep4, attributes_gep4, needsRoot = TRUE)
p_gep4 <- plotCircDend(network_gep4)
ggsave(plot=p_gep4, "results/scenic_gep4.pdf", width=5.5, height=5.5)

## GEP 6
weights_gep6 <- read.csv("data/scenic_gep_data_with_threshold_50target_70regulon/GEP6_seurat_filtered_harmony_02_15_23_raw_counts_multirun_04192023_regulons_pruned_target_gene_counts.csv",
                         row.names = 1)
attributes_gep6 <- read.csv("data/Scenic_gep_data_with_threhold_30target_70regulon/GEP_6_pruned.csv", row.names = 1)

network_gep6 <- scenic2network(weights_gep6, attributes_gep6)
p_gep6 <- plotCircDend(network_gep6)
ggsave(plot=p_gep6, "results/scenic_gep6.pdf", width=5.5, height=5.5)

