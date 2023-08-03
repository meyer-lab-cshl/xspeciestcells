library(readr)
library(tidyverse)
library(dplyr)
library(readxl)
library(igraph)
library(rgexf)
library(RColorBrewer)
library(scales)

#x11(type="cairo")


##################################################################################
########################### DEFINE VARIABLES AND INPUTS ##########################
##################################################################################
# this the reproducible counts summarized across all seed runs for every regulon and every reproducible target gene from each SCENIC run; to
# obtain this csv, you will need to run: STEP1_reproducible_regulons_and_gene_targets_across_runs.py (Tonya can provide if needed)
# column names = regulon; row names = target genes
weights <- read.csv("data/scenic_gep_data_with_threshold_50target_70regulon/GEP4_seurat_filtered_harmony_02_15_23_raw_counts_multirun_04192023_regulons_pruned_target_gene_counts.csv", row.names = 1)

# 3 column csv file:
#   - column 1 = name of regulon (column name is "regulon")
#   - column 2 = name of target gene (column name is "target")
#   - column 3 = name of attribute/characteristic/class (column name is "attribute")
# ex: column 1 tells you the regulon name, column 2 tells you which target gene the regulon is acting upon (arrow coming out of regulon node pointing to target gene node),
# column 3 tells you the attribute or class; if I listed CEBPD in column 1, that means it is a regulon, if in column 2, I list
# GZMK, this means CEBPD regulates the target gene GZMK (so arrow coming from CEBPD node and pointing to GZMK node), if I list "Cytotoxicity"in column 3, that means I am give GZMK the class name "Cytotoxicity"
# (colors the node and network arrow line all the same color for this class)
attributes <- read.csv("data/Scenic_gep_data_with_threhold_30target_70regulon/GEP_4_pruned.csv", row.names = 1)



##################################################################################
##################################### RUN CODE ###################################
######################## SEE LINE 85 TO SET MANUAL COLORS ########################
##################################################################################

weights <- replace(weights, is.na(weights), values = 0) # replace NAs with 0
weights[weights < 3] <- 0
weights$dummy <- 10

names(attributes)[1] <- "regulon"
names(attributes)[2] <- "target"
names(attributes)[3] <- "attribute"
all_regulons <- unique(attributes$regulon)
all_targets <- unique(attributes$target)
pairwise_source_target <- data.frame(matrix(ncol = 3, nrow = 0)) # dataframe to store all regulon-target pairwise comparisons
colnames(pairwise_source_target) <- c("regulon", "target", "weight")

# generate all regulon-target pairwise comparisons
for (i in all_regulons){
  for (j in all_targets){
    tmp <- data.frame(regulon = i, target = j, weight = weights[j,i])
    pairwise_source_target <- rbind(pairwise_source_target, tmp)
  }
}

# subset dataframe to regulon-target pairs that do not have 0 as their reproducible weight
non_zero_weights <- pairwise_source_target[which(pairwise_source_target$weight != 0),]

df_weights <- tibble::tibble(regulon='dummy', target=unique(non_zero_weights$regulon), weight=10)
non_zero_weights <- rbind(non_zero_weights, df_weights)

# add metadata
subset_attr <- attributes[,c("target", "attribute")]
df_attr <- tibble::tibble(target='dummy', attribute="Other")
subset_att <- rbind(subset_attr, df_attr)
non_zero_weights_with_attrs <- merge(non_zero_weights, subset_attr, by = "target", all.x = T)


# a dataframe where each row is a node in the network and the value is the vertex number (helps with reordering and aesthetics)
target_to_att_index_map <- unique(non_zero_weights_with_attrs[,c('target', 'attribute')])
missing_nodes <- unique(non_zero_weights_with_attrs$regulon)[which(unique(non_zero_weights_with_attrs$regulon) %in% target_to_att_index_map$target == FALSE)]
missing_nodes_df <- data.frame(target = missing_nodes, attribute = rep("TF", length(missing_nodes)))
target_to_att_index_map <- rbind(target_to_att_index_map, missing_nodes_df)


# Auto select colors you want; if you want to manually set colors, go to line 85
coul  <- brewer.pal(length(unique(target_to_att_index_map$attribute)), "Set3")
names(coul) <- levels(as.factor(target_to_att_index_map$attribute))
target_to_att_index_map$color <- NA
target_to_att_index_map$red <- NA
target_to_att_index_map$green <- NA
target_to_att_index_map$blue <- NA
for (colorAdd in names(coul)){
  print(colorAdd)
  target_to_att_index_map[which(target_to_att_index_map$attribute == colorAdd), "color"] <-  coul[colorAdd]
  target_to_att_index_map[which(target_to_att_index_map$attribute == colorAdd), "red"] <- col2rgb(coul[colorAdd])[1]
  target_to_att_index_map[which(target_to_att_index_map$attribute == colorAdd), "green"] <- col2rgb(coul[colorAdd])[2]
  target_to_att_index_map[which(target_to_att_index_map$attribute == colorAdd), "blue"] <- col2rgb(coul[colorAdd])[3]
}

c("")
# manually update colors (repeat this for each color you want to change)
target_to_att_index_map[which(target_to_att_index_map$attribute == "TF"), "color"] <- "#FFD382"
target_to_att_index_map[which(target_to_att_index_map$attribute == "Cytokine"), "color"] <- "#964A35"
target_to_att_index_map[which(target_to_att_index_map$attribute == "Chemokine"), "color"] <- "#964A35"
target_to_att_index_map[which(target_to_att_index_map$attribute == "Other"), "color"] <- "gray"
target_to_att_index_map[which(target_to_att_index_map$attribute == "Cytokine Receptor"), "color"] <- "#8EB63B"
target_to_att_index_map[which(target_to_att_index_map$attribute == "Cytotoxicity"), "color"] <- "#C7B354"
target_to_att_index_map[which(target_to_att_index_map$attribute == "Chemokine Receptor"), "color"] <- "#001959"
target_to_att_index_map[which(target_to_att_index_map$attribute == "Integrin"), "color"] <- "#45204C"
target_to_att_index_map[which(target_to_att_index_map$attribute == "NK receptor"), "color"] <- "#FFACAC"

for (rownumber in 1:nrow(target_to_att_index_map)){
  print(rownumber)
  hex_values <- target_to_att_index_map[rownumber,"color"]
  rgb_list <- col2rgb(hex_values)
  target_to_att_index_map[rownumber, "red"] <- rgb_list[1]
  target_to_att_index_map[rownumber, "green"] <- rgb_list[2]
  target_to_att_index_map[rownumber, "blue"] <- rgb_list[3]
}

# manually update the edges by RGB colors
#target_to_att_index_map[which(target_to_att_index_map$attribute == "Cytokine"), "red"] <- 80
#target_to_att_index_map[which(target_to_att_index_map$attribute == "Cytokine"), "green"] <- 200
#target_to_att_index_map[which(target_to_att_index_map$attribute == "Cytokine"), "blue"] <- 120
#target_to_att_index_map[which(target_to_att_index_map$attribute == "TF"), "red"] <- 0
#target_to_att_index_map[which(target_to_att_index_map$attribute == "TF"), "green"] <- 0
#target_to_att_index_map[which(target_to_att_index_map$attribute == "TF"), "blue"] <- 0

coul <- target_to_att_index_map$color
names(coul) <- target_to_att_index_map$attribute
# check palette
pal <- rgb(target_to_att_index_map$red, target_to_att_index_map$green, target_to_att_index_map$blue, max= 255)
show_col(pal)

# order by attribute
target_to_att_index_map <- target_to_att_index_map[order(target_to_att_index_map$attribute, decreasing = T),]

# buid the network graph from the data frame
network <- graph_from_data_frame(d = non_zero_weights, directed = T)
rownames(target_to_att_index_map) <- target_to_att_index_map$target
target_to_att_index_map <- target_to_att_index_map[names(as.list(V(network))),]
network <- set_vertex_attr(network, name = "group", index = V(network), target_to_att_index_map$attribute)

# get the number of degrees of connection for each node coming **out** of the node (ignoring pointing in -- can be changes with the "mode =" argument)
deg <- degree(network, mode = "out", normalized = F) # doesn't matter of us because I changed the plot to not use degree as size of node
#par(mar=c(0,0,0,0)+.1)

for (relationship in seq(length(E(network)))){
  regulon = get.edgelist(network)[relationship,1]
  target = get.edgelist(network)[relationship,2]
  get_info = target_to_att_index_map[which(target_to_att_index_map$target == target),]
  alpha = E(network)$weight[relationship]*0.05
  #E(network)[relationship]$color = rgb(get_info$red, get_info$green, get_info$blue, alpha = alpha, max = 255)
  E(network)[relationship]$color = rgb(get_info$red, get_info$green, get_info$blue, max = 255)

}
# update colors
coul <- target_to_att_index_map$color
names(coul) <- target_to_att_index_map$attribute
# order by attribute
order <- target_to_att_index_map[order(target_to_att_index_map$attribute, decreasing = T),"target"]

# plots the circular network graph
plot(network, vertex.size = 0, edge.width = E(network)$weight*0.25,
     vertex.color = coul,
     layout = layout_in_circle(network, order = t(as.data.frame(as.list(V(network))))[order,]),
     edge.arrow.size= 0.5, magin = 0)

legend("topleft",bty = "n",
       legend=levels(as.factor(target_to_att_index_map$attribute)),
       fill=coul[unique(names(coul))][levels(as.factor(target_to_att_index_map$attribute))], border=NA, cex = 0.75, y.intersp = 0.75)


# plots circular dendrogram


lay <- create_layout(network, layout = 'dendrogram', circular = T)
head(lay)
ggraph(lay) +
  geom_edge_diagonal() +
  geom_node_point(aes(filter = leaf, color = group)) +
  scale_color_manual(limits = as.factor(unique(lay$group)), values = coul[unique(names(coul))]) +
  geom_node_label(aes(label = ifelse(leaf == FALSE, name, NA)), repel = TRUE) +
  #geom_node_label(aes(label = name), repel = TRUE) +
  #geom_node_text(aes(label = name), repel = TRUE, check_overlap = TRUE, angle = 45, hjust = 1.5, vjust = 1.5) +  # Add this line for labels
  coord_fixed() +
  theme(panel.background = element_blank())

lay <- create_layout(network, layout = 'dendrogram', circular = T)
head(lay)
ggraph(lay) +
  geom_edge_diagonal() +
  geom_node_point(aes(filter = leaf, color = group)) +
  #scale_color_manual(limits = as.factor(unique(lay$group)), values = coul[unique(names(coul))]) +
  #geom_node_label(aes(label = ifelse(leaf == FALSE, name, NA)), repel = TRUE) +
  geom_node_label(aes(label = name), repel = TRUE) +
  #geom_node_text(aes(label = name), repel = TRUE, check_overlap = TRUE, angle = 45, hjust = 1.5, vjust = 1.5) +  # Add this line for labels
  coord_fixed() +
  theme(panel.background = element_blank())

network[[2]]
