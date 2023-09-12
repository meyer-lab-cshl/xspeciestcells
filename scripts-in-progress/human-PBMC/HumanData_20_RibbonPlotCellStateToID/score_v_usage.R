library(ggalluvial)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)

seur.geps <- readRDS("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/seuratobj_gepscores_allgenes_2023-08-04.rds")
#usages <- read_delim("data/human-thymus/HumanData_20_RibbonPlotCellStateToID/imputed_cNMF.usages.k_12.dt_0_02.consensus.txt", delim="\t",
#                     col_names = c("cellid", paste0("GEP", 1:12)), skip=1)

usages <- read_delim("data/human-thymus/HumanData_20_RibbonPlotCellStateToID/imputed_cNMF_allcells.usages.k_13.dt_0_02.consensus.txt", delim="\t",
                     col_names = c("cellid", paste0("GEP", 1:13)), skip=1)

usages_wide <- usages %>%
  mutate(usage_assigned = paste0("GEP", apply(usages[, 2:13], 1, which.max)))

usages <- usages_wide %>%
  pivot_longer(starts_with("GEP"), values_to="usage", names_to="name")

gep_pbmc <- c("GEP1", "GEP4", "GEP5", "GEP6", "GEP8", "GEP12")
#gep_pbmc <- c("GEP1", "GEP4", "GEP5", "GEP6", "GEP8", "GEP12")
df_wide <- as.data.frame(seur.geps@meta.data[seur.geps@meta.data$Tissue=="PBMC",gep_pbmc])
df <- df_wide %>%
  rownames_to_column("cellid") %>%
  pivot_longer(cols=all_of(gep_pbmc), 
               names_to="name", values_to="score")
# Check the distribution
ggplot(df)+
  geom_histogram(aes(x=score), bins = 100) +
  facet_wrap(~name, ncol=4)+
  # _______
  # geom_violin(aes(x=gep, y=score), width=1)+
  # geom_boxplot(aes(x=gep, y=score), outlier.shape = NA, width=0.05)+
  # geom_jitter(aes(x=gep, y=score), size=0.1, width = 0.05)+
  # _______
  labs(y="GEP score", title="Raw GEP score")



combined <- df %>%
  left_join(usages, by=c("cellid", "name")) %>%
  as_tibble

combined <- combined %>%
  filter(name %in% paste0("GEP", c(1,4,5,6))) %>% 
  filter(usage_assigned %in% paste0("GEP", c(1,4,5,6))) %>%
  group_by(cellid) %>%
  mutate(assign_TRUE = usage_assigned == name) %>%
  ungroup

ggplot(combined) +
  geom_density(aes(usage,  color=usage_assigned==name)) +
    facet_wrap(~name, scales="free")

tt <- combined %>%
  arrange(name, usage)

ggplot(combined) +
  geom_point(aes(x=cellid, y=usage, color=usage_assigned==name)) +
  facet_wrap(~name, scales="free")+
  theme(axis.text.x = element_blank())

ggplot(combined) +
  geom_density(aes(score_pmax)) +
  facet_wrap(~name, scales="free") 

ggplot(combined) +
  geom_density(aes(score)) +
  facet_wrap(~name, scales="free") 

ggplot(combined) +
  geom_point(aes(x=usage, y=score_pmax, color=usage_assigned==name)) +
  scale_color_brewer(type="qual") +
  facet_wrap(usage_assigned==name~name, nrow=2)


combined %>%
  filter(name== "GEP1") %>%
  select(usage, assign_TRUE) %>%
  pivot_wider(values_from="usage", names_from="assign_TRUE")
  
## ribbon plot ####
usages_wide <- usages_wide %>% 
  filter(grepl("PBMC", cellid)) %>%
  column_to_rownames("cellid") %>%
  #select(paste0("GEP", c(1,4,5,6))) 
  select(paste0("GEP", c(2,4,5,6,7))) 


usages_wide$score_asign <- apply(usages_wide, 1, max)
usages_wide$gep_assign <- colnames(usages_wide)[apply(usages_wide, 1, which.max)]
usages_wide$score_2ndmax <- sapply(1:nrow(usages_wide), function(x) {
  names(sort(usages_wide[x,1:4], decreasing = T)[2])
  #colnames(usages_wide[,1:4])[which(tmp == usages_wide[x,1:4])]
})
usages_wide$score_2ndmax_asign <- colnames(usages_wide)[usages_wide[,1:4] == usages_wide$score_2ndmax]


table(rownames(seur.geps@meta.data[seur.geps$Tissue=="PBMC",]) == rownames(usages_wide))
seur.geps@meta.data$gep_assign[seur.geps$Tissue=="PBMC"] <- usages_wide$gep_assign

cellids <- read_csv("scripts-in-progress/human-PBMC/HumanData_20_RibbonPlotCellStateToID/cell.ID.after.removal.csv")[,2]

cd8_pbmc <- seur.geps@meta.data %>%
  filter(group.ident == "CD8_PBMC")
cd8_pbmc_remove <- rownames(cd8_pbmc)[!rownames(cd8_pbmc) %in% cellids$x]

counts <- seur.geps@meta.data %>%
  rownames_to_column("cellid") %>%
  as_tibble() %>%
  filter(Tissue=="PBMC") %>%
  filter(!cellid %in% cd8_pbmc_remove) %>%
  # get nb of cells per gep assignment
  group_by(cell.ident, gep_assign) %>%
  summarise(ncells=n()) %>%
  # get %cells in each gep assignment
  ungroup() %>%
  group_by(cell.ident) %>%
  mutate(totalcells=sum(ncells),
         freq = ncells*100/totalcells) %>%
  ungroup() %>%
  # rename a few variables
  mutate(gep_assign=toupper(gep_assign),
         cell.ident=replace(cell.ident, cell.ident=="NKT", "iNKT"))

counts %>%
  mutate(gep_assign=replace(gep_assign, !gep_assign%in%c("GEP2", "GEP4", "GEP5", "GEP6", "GEP7"), "other"),
         gep_assign = factor(gep_assign, levels=c("GEP2", "GEP4", "GEP5", "GEP6", "GEP7", "other"))) %>%
  filter(gep_assign != "other") %>%
  ggplot(aes(axis1=cell.ident, axis2=gep_assign, y=freq)) +
  geom_alluvium(aes(fill=cell.ident))+
  geom_stratum()+
  geom_text(stat="stratum", aes(label=after_stat(stratum)), size=8)+
  scale_fill_manual(values=c("#2ca25f", "#dd1c77", "#045a8d", "#8856a7", "#bdc9e1"), name="cell type")+
  # scale_x_discrete(limits=c("Cell type", "GEP")) + theme_classic()
  theme_void()+
  theme(legend.position="none")

library(GGally)
ggpairs(usages_wide[,1:4])
       

ggplot(usages_wide) +
  geom_bar(aes(score_2ndmax)) +
  facet_wrap(~gep_assign, scales="free")
