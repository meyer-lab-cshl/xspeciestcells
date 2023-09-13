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
  as_tibble %>%
  mutate(assign_TRUE = usage_assigned == name) %>%
  ungroup

p_score <- ggplot(combined) +
  geom_density(aes(score)) +
  facet_wrap(~name, scales="free") 
ggsave(plot=p_score , "data/human-thymus/HumanData_20_RibbonPlotCellStateToID/score_dist.pdf")


p_usagevscore <- ggplot(combined) +
  geom_point(aes(x=usage, y=score, color=usage_assigned==name)) +
  scale_color_brewer(type="qual") +
  facet_wrap(usage_assigned==name~name, nrow=2)
ggsave(plot=p_usagevscore , "data/human-thymus/HumanData_20_RibbonPlotCellStateToID/usage_v_score.pdf")

## Remove GEP12 as only 1 cell assigned
combined <- combined %>%
  filter(name %in% paste0("GEP", c(1,4,5,6))) %>% 
  filter(usage_assigned %in% paste0("GEP", c(1,4,5,6))) %>%
  group_by(cellid) %>%
  mutate(assign_TRUE = usage_assigned == name) %>%
  ungroup


## find usage thresholds ####
assigned_true <- combined %>%
  filter(assign_TRUE)

assigned_false <- combined %>%
  filter(!assign_TRUE)

intersection.points <- sapply(unique(combined$name), function(gep) {
  limits <- range(filter(combined, name==gep)$usage)
  
  at <- density(filter(assigned_true, name == gep)$usage, from=limits[1],
                to=limits[2], n=2^10)
  af <- density(filter(assigned_false, name == gep)$usage, from=limits[1],
                to=limits[2], n=2^10)
  dens.diff <- at$y-af$y
  median(at$x[which(diff(dens.diff > 0) != 0) + 1])
})

intsec <- tibble(name=names(intersection.points), intersection=intersection.points)

## check that intersects found correctly
p_usage_thr <- ggplot() +
  geom_density(data=combined, aes(usage,  color=usage_assigned==name)) + 
  geom_vline(data=intsec, aes(xintercept = intersection)) +
  facet_wrap(~name, scales="free")
ggsave(plot=p_usage_thr , "data/human-thymus/HumanData_20_RibbonPlotCellStateToID/usage_with_percentilethr.pdf")


## find percentile of intersect threshold
intersection_percentile <- sapply(unique(combined$name), function(gep) { 
  sum(filter(combined, name == gep)$usage < filter(intsec, name==gep)$intersection)/nrow(filter(combined, name == gep))
})

intsec_percentile <- tibble(name=names(intersection_percentile), intersection_percentile=intersection_percentile)
write_csv(intsec_percentile, "data/human-thymus/HumanData_20_RibbonPlotCellStateToID/score_percentile_thr.csv")

intsec_score <- combined %>%
  inner_join(intsec_percentile, by="name") %>%
  group_by(name) %>%
  summarise(intersection_score=quantile(score, unique(intersection_percentile)))

## check where intersect scores based on intersect usage percentile are located
p_scores <- ggplot() +
    geom_histogram(data=combined, aes(x=score), bins = 100)+ 
    geom_vline(data=intsec_score, aes(xintercept = intersection_score)) +
    facet_wrap(~name, scales="free",  ncol=4) 
ggsave(plot=p_scores, "data/human-thymus/HumanData_20_RibbonPlotCellStateToID/score_with_percentilethr.pdf")


combined <- combined %>%
  # get back to wide
  select(cellid, name, score) %>%
  pivot_wider(names_from=name, values_from = score, values_fill=0) %>%
  # thresholds
  mutate(GEP1_threshold=ifelse(GEP1>intsec_score$intersection_score[intsec_score$name=="GEP1"], T, F),
         GEP4_threshold=ifelse(GEP4>intsec_score$intersection_score[intsec_score$name=="GEP4"], T, F),
         GEP5_threshold=ifelse(GEP5>intsec_score$intersection_score[intsec_score$name=="GEP5"], T, F),
         GEP6_threshold=ifelse(GEP6>intsec_score$intersection_score[intsec_score$name=="GEP6"], T, F)) %>%
  # pivot longer
  pivot_longer(cols=ends_with("threshold"), names_to="name", values_to="pass") %>%
  mutate(pass=as.numeric(pass),
         name=gsub("_threshold", "", name)) %>%
  # keep only lines with GEPs that passed threshold (1 gep assigned)
  group_by(cellid) %>%
  filter(pass==1 & sum(pass)==1) %>%
  # keep only columns of interest
  select(cellid, name) %>%
  dplyr::rename(score_assigned=name) %>%
  distinct() %>%
  inner_join(combined, by="cellid") %>%
  select(cellid, name, score, usage, score_assigned, usage_assigned)


p_corr <- ggplot(filter(combined,name!="GEP8")) +
  geom_hex(aes(x=usage, y=score), bins=100) +
  scale_fill_continuous(type = "viridis") +
  facet_wrap(usage_assigned==score_assigned~name, nrow=2)


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

ggsave(plot=p_corr, "data/human-thymus/HumanData_20_RibbonPlotCellStateToID/usage_score_thr_correspondence.pdf")
