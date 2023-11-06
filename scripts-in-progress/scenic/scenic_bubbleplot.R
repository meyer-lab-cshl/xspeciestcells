#################
## libraries ####
#################

library(readr)
library(tidyverse)
library(dplyr)
library(readxl)
library(RColorBrewer)

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
attributes <- read.csv("results/Pruned_Genes_Adaptive_regulons_w_attributes.csv")

attributes_TF_TF <- attributes %>%
  filter(Target %in% unique(TF))

attributes_uniqueTF <- attributes %>%
  filter(!Target %in% unique(TF))

attributes_matrix <- as.matrix(table(attributes[,1:2]))
rowclust <- hclust(dist(attributes_matrix, "manhattan"))$order
colclust <- hclust(dist(t(attributes_matrix), "manhattan"))$order

attributes_matrix <- attributes_matrix[rowclust,colclust]

dummy_target <- attributes %>%
  select(Target, Attributes) %>%
  distinct %>%
  rename(Attributes_target=Attributes)

dummy_tf <- attributes %>%
  select(TF, Attributes) %>%
  distinct %>%
  rename(Attributes_TF=Attributes)

attributes_df <- attributes_matrix %>%
  as.data.frame %>%
  filter(Freq != 0) %>%
  left_join(dummy_target) %>%
  left_join(dummy_tf) %>%
  mutate(type = case_when(Attributes_TF == Attributes_target ~ "TF",
                          TRUE ~ "other"))

p_reg <- ggplot(attributes_df)
p_reg +
  geom_point(aes(x=Target, y=TF, color=type), size=4) +
  labs(y="Transcription Factors",
       x="Targets") +
  scale_color_manual(values=c("darkgrey", "#ec7014"), guide="none")+
  coord_fixed() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))



