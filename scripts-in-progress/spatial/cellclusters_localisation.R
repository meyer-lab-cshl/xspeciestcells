library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(pals)
library(ggh4x)

# Import
summary.annotated <- read.csv("~/Projects/HumanThymusProject/data/human-thymus/SpatialData/proportion_cluster_across_regions.csv")

# cluster levels
clusters_levels1=c("thyCD4_ISP", "thyCD4_DPp", "thyCD4_DPq", "thyCD4_ccr9", "thyCD4_ccr7", "thyCD4_Tagonist", "thyCD4_Treg",
                  "thyCD8_DP", "thyCD8_ccr9", "thyCD8_idk", "thyCD8_ccr7", "thyCD8_cd8aa1", "thyCD8_cd8aa2",
                  "thyMAIT_DP", "thyMAIT_cd8aa", "thyMAIT_ccr9", "thyMAIT_tbc", "thyMAIT_ccr7", "thyMAIT_IFNsig", "thyMAIT_effector",
                  "thyNKT_cd8aa", "thyNKT_ccr9", "thyNKT_ccr7", "thyNKT_effFOSJUN", "thyNKT_effector",
                  "thyGDT_DP", "thyGDT_immat_cycl", "thyGDT_immat", "thyGDT_ccr9", "thyGDT_IFNsig", "thyGDT_effector")

# colors
colhisto <- brewer.pal(7, "Dark2")
location_levels <- c("Cortex", "Junction", "Medulla", "Interlobular", "Inflammatory", "Hemoglobin rich", "Hassall associated")

# design <- rbind(c(1:7), c(8:13,NA), c(14:19,NA), c(20:26), c(27:31, NA, NA))
ggplot(summary.annotated,
       aes(x=factor(cluster, level=clusters_levels1), y=proportion_cluster_across_regions, fill=factor(location, level=location_levels))) +
  geom_bar(stat='identity', position="fill") +
  scale_fill_manual(values=colhisto, name="Location") +
  labs(x="", y="% cells in each location")+
  theme_cowplot()+
  # facet_manual(~cluster, design) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
ggsave("~/Projects/HumanThymusProject/data/human-thymus/SpatialData/proportion_clusters_across_locations.jpeg", width=12, height=6)

# ggsave(plot=p_bar,
#        filename = 'results/plots/proportion_cluster_across_regions.pdf',
#        width = 15, height=12)