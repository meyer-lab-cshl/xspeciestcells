seur.geps <- readRDS("./data/human-thymus/HumanData_20_RibbonPlotCellStateToID/seuratobj_gepscores_allgenes_2023-08-04.rds")
usages <- read_delim("data/human-thymus/HumanData_20_RibbonPlotCellStateToID/imputed_cNMF.usages.k_12.dt_0_02.consensus.txt", delim="\t",
                     col_names = c("cellid", paste0("GEP", 1:12)), skip=1)

usages <- usages %>%
  mutate(usage_assigned = paste0("GEP", apply(usages[, 2:13], 1, which.max))) %>%
  pivot_longer(starts_with("GEP"), values_to="usage", names_to="name")

df <- as.data.frame(seur.geps@meta.data[seur.geps@meta.data$Tissue=="PBMC",gep_pbmc])
df <- df %>%
  rownames_to_column("cellid") %>%
  pivot_longer(cols=all_of(gep_pbmc), 
               names_to="name", values_to="score")
# Check the distribution
ggplot()+
  geom_histogram(aes(x=score), bins = 100)+
  facet_wrap(~gep, ncol=4)+
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
  group_by(name) %>%
  mutate(score_nn = ifelse(score < 0, 0, score),
         score_pmax = score_nn/max(score_nn)) %>%
  group_by(cellid) %>%
  mutate(score_assigned = max(score_pmax)) %>%
  filter(name %in% paste0("GEP", c(1,4,5,6)))

ggplot(combined) +
  geom_density(aes(usage)) +
    facet_wrap(~name, scales="free")

tt <- combined %>%
  arrange(name, usage) %>%
  mutate(cellid=fct_inorder(cellid))

ggplot(tt) +
  geom_point(aes(x=cellid, y=usage, color=usage_assigned==name)) +
  facet_wrap(~name, scales="free")+
  theme(axis.text.x = element_blank())

ggplot(combined) +
  geom_density(aes(score_pmax)) +
  facet_wrap(~name, scales="free") 

ggplot(combined) +
  geom_point(aes(x=usage, y=score_pmax, color=usage_assigned==name)) +
  scale_color_brewer(type="qual") +
  facet_wrap(usage_assigned==name~name, nrow=2)

             