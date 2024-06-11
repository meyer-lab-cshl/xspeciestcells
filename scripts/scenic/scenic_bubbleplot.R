#################
## libraries ####
#################

library(readr)
library(tidyverse)
library(dplyr)
library(readxl)
library(RColorBrewer)
library(patchwork)
#############
## setup ####
#############

colors <- tibble::tibble(attribute=c("TF", "cytokine", "chemokine", "other",
                                     "cytokine receptor", "cytotoxicity",
                                     "chemokine receptor", "integrin",
                                     "NK receptor"),
                         color=c("#ec7014", "#964A35", "#964A35", "gray",
                                 "#8EB63B", "#C7B354", "#001959", "#45204C",
                                 "#FFACAC"))

#"#FFD382"
colors_list <- as.list(colors$color)
names(colors_list) <- colors$attribute

#################
## functions ####
#################

scenicBubble <- function(attributes, ordered) {
  attributes_target <- attributes %>%
    select(Target, Attributes) %>%
    distinct
  
  if (ordered == "hierarchical") {
    attributes_matrix <- as.matrix(table(attributes[,1:2]))
    rowclust <- hclust(dist(attributes_matrix, "manhattan"))$order
    colclust <- hclust(dist(t(attributes_matrix), "manhattan"))$order
    
    attributes_matrix <- attributes_matrix[rowclust,colclust]
    
    # Hierarchically clustered
    attributes_df <- attributes_matrix %>%
      as.data.frame %>%
      mutate(Target= fct_inorder(Target),
             TF=fct_inorder(TF)) %>%
      filter(Freq != 0) %>%
      left_join(attributes_target) %>%
      as_tibble
  } else {
    # ordered by categories
    attributes_df <- table(attributes[,1:2]) %>%
      as.data.frame %>%
      filter(Freq != 0) %>%
      left_join(attributes_target) %>%
      as_tibble
  }
  
  tf_tf <- attributes_df  %>%
    rowwise %>%
    mutate(xy=paste(sort(c(as.character(TF), as.character(Target))), collapse="_")) %>%
    group_by(xy) %>%
    summarize(count = n()) %>%
    filter(count > 1) %>%
    ungroup() %>%
    separate(xy, into = c("TF", "Target"))
  
  tf_tf <- tf_tf %>%
    bind_rows(tibble(TF=tf_tf$Target, Target=tf_tf$TF, count=tf_tf$count))
  
  attributes_df <- attributes_df %>%
    left_join(tf_tf) %>%
    mutate(type= case_when(Attributes == "TF" & Target %in% unique(TF) ~ "TF - unidirectional",
                           TRUE ~ "Other"),
           type = case_when(!is.na(count) | TF == Target ~ "TF - bidirectional", 
                            TRUE ~ type))
  
  # Proportion of genes in regulon
  targets_per_tf <- attributes_df %>%
    select(TF) %>%
    group_by(TF) %>%
    summarise(count=n())
  
  tf_per_target <- attributes_df %>%
    select(Target) %>%
    group_by(Target) %>%
    summarise(count=n()) %>%
    left_join(distinct(attributes_df[, c(2,4)])) %>%
    arrange(Attributes, -count) %>%
    mutate(Target= fct_inorder(Target))
  
  p_tf_per_target <- ggplot(tf_per_target,
                            aes(y=count, x=Target, fill=Attributes)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=colors_list) +
    scale_y_continuous(breaks=seq(0,12,2)) +
    cowplot::theme_cowplot() +
    labs(y="# TFs")+
    theme(axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")
  
  if (ordered != "hierarchical") {
    attributes_df <- attributes_df %>%
      mutate(Target= factor(Target, levels=levels(tf_per_target$Target)),
             TF=fct_inorder(TF)) 
  }
  
  
  p_reg <- ggplot(attributes_df)
  p_reg <- p_reg +
    geom_point(aes(x=Target, y=TF, color=type), size=4) +
    labs(y="Transcription Factors",
         x="Targets") +
    scale_color_manual(values=c("darkgrey", "#7570b3",  "#1b9e77"))+
    coord_fixed() +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
          legend.position = "bottom",
          panel.grid.major = element_line(color="lightgrey", size=0.2),
          panel.border = element_rect(color="black", linewidth=1),
          axis.line = element_blank())
  
  
  
  p_tf_per_target /p_reg  + plot_layout(heights = unit(c(3, -1), c("cm", "null")))
}

##########################
## INPUTS Description ####
##########################

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

## Regulons conventional T cells
t_conv <- read.csv("results/Pruned_Genes_Conventional_regulons_w_attributes.csv")
t_inn <- read.csv("results/Pruned_Genes_Innate_regulons_w_attributes.csv")

p_conv <- scenicBubble(t_conv, ordered="count")
p_inn <- scenicBubble(t_inn, ordered="count")

ggsave(plot=p_conv, "results/scenic_Tconv.pdf",
       width=9, height=7)
ggsave(plot=p_inn, "results/scenic_Tinn.pdf",
       width=9, height=7)
####

attributes_target <- attributes %>%
  select(Target, Attributes) %>%
  distinct

attributes_tf <- attributes %>%
  #select(TF, Attributes) %>%
  distinct %>%
  rename(Attributes_TF=Attributes)

attributes_TF_TF <- attributes %>%
  filter(Target %in% unique(TF))

attributes_uniqueTF <- attributes %>%
  filter(!Target %in% unique(TF))

attributes_TF_TF <- attributes %>%
  mutate(TF_TF = case_when(Target %in% unique(TF) ~ "TF",
                           TRUE ~ "other"))

attributes_df <- attributes_df %>%
  mutate(type = case_when((tf_tf$Target == Target &  tf_tf$TF == TF) | 
                            (tf_tf$TF == Target &  tf_tf$Target == TF) ~ "TF-TF",
                          TRUE ~ "Other"),
         TF_fill = 0.5,
         Target_fill = 0.5,
         Target_fill=case_when(type == "TF" ~ 0,
                               TRUE ~ Target_fill),
         TF_fill= case_when(Target_fill == 0 ~ 1,
                            TRUE ~ TF_fill))

att <- attributes_df %>%
  mutate(TF_num = as.numeric(TF),
         Target_num = as.numeric(fct_inorder(Target)))


p <- ggplot() +
  geom_scatterpie(data=att, aes(x=Target_num, y=TF_num),
                  cols=c("TF_fill", "Target_fill"), pie_scale=0.4, color=NA) +
  labs(y="Transcription Factors",
       x="Targets") +
  scale_y_continuous(breaks=sort(unique(att$TF_num)),
                     labels=levels(att$TF)) + 
  scale_x_continuous(breaks=sort(unique(att$Target_num)),
                     labels=levels(fct_inorder(att$Target))) + 
  scale_fill_manual(values=c("#ec7014", "darkgrey"), guide="none") +
  coord_fixed() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        panel.grid.minor = element_blank())


