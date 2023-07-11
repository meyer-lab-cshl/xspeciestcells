###
# Purpose: Combine metaneighbor bubble plots (NKT-MAIT-GD)
# Date: May 10th 2023
# Author: Salomé Carcy
###




## 1. IMPORT ####
library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)


# Import data
mtn.inkt <- readRDS("./data/cross-species/04_Metaneighbor_nkt/nkt_ms-hu_mtnslowversion_DF.rds")
mtn.mait <- readRDS("./data/cross-species/04_Metaneighbor_mait/mait_ms-hu_mtnslowversion_DF.rds")
mtn.gdt  <- readRDS("./data/cross-species/04_Metaneighbor_gdt/gdt_mssagar-hu_mtnslowversion_DF_2023-07-10.rds")



## 2. iNKT PLOT ####

# PROPORTION OF HUMAN NKT CELLS IN EACH CLUSTER
inkt.bpX <- ggplot(data=mtn.inkt%>% select(human,propcells_human) %>% distinct(),
                   aes(x=factor(human, levels=paste0("iNKT_c", 0:6)), y=propcells_human))+
  geom_bar(stat="identity", fill="#bdbdbd") + theme_cowplot()+
  scale_x_discrete(position="top")+
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100))+
  labs(y="%cells")+
  theme(axis.text.x  = element_text(size = 15, angle=45, hjust=0),
        axis.title.x = element_blank(),
        axis.line.x  = element_blank(),
        axis.text.y  = element_text(size=15),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size=20),
        legend.position = "none")

# PROPORTION OF MOUSE NKT CELLS IN EACH CLUSTER
order_inkt <- c("Stage0", "iNKTp", "iNKT2", "iNKT17", "iNKT1")
inkt.bpY <- ggplot(data=mtn.inkt%>% select(mouse,propcells_mouse) %>% distinct(),
                   aes(x=factor(mouse, levels=rev(order_inkt)), y=propcells_mouse))+
  geom_bar(stat="identity", fill="#bdbdbd") +
  scale_x_discrete(position="top") +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100))+
  labs(y="%cells")+ coord_flip() + theme_cowplot()+
  theme(axis.text    = element_text(size=15),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size=20),
        axis.title.y = element_blank(),
        axis.line.y  = element_blank(),
        legend.position = "none")

# BUBBLE PLOT
inkt.hm <- ggplot(mtn.inkt, aes(x=factor(human, levels=paste0("iNKT_c", 0:6)),
                                y=factor(mouse, levels=rev(order_inkt)))) +
  geom_point(aes(size = auroc, color= auroc))+
  geom_text(data=mtn.inkt %>% filter(auroc>0.65) %>% mutate(across("auroc", \(x) round(x,2))), aes(label=auroc), color="white")+
  scale_size_continuous(limits=c(0,1), breaks=seq(0.2,0.8, by=0.2), range = c(1, 15))+
  scale_color_gradient2(low="#d9d9d9", mid="white", high="#a50f15", midpoint=0.5, limits=c(0,1), name="AUROC", breaks=seq(0,1, by=0.2))+
  labs(x="Human clusters", y="Mouse clusters", size="AUROC")+
  theme_cowplot()+
  theme(legend.position  = "bottom",
        legend.key.width = unit(0.8, 'cm'),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        axis.text.y = element_text(size=15))


# COMBINE
inkt.tot <- (inkt.bpX+plot_spacer() + plot_layout(widths = c(5, 1))) / (inkt.hm + inkt.bpY + plot_layout(widths = c(5, 1))) + plot_layout(heights = c(1, 5))
# ggsave("./data/cross-species/04_Metaneighbor_nkt/nkt_ms-hu_metaneighbor_bubbleplot4.svg", width=9, height=8)




## 3. MAIT PLOT ####

# PROPORTION OF HUMAN MAIT CELLS IN EACH CLUSTER
mait.bpX <- ggplot(data=mtn.mait %>% select(human,propcells_human) %>% distinct(),
                   aes(x=factor(human, levels=paste0("MAIT_c", 0:6)), y=propcells_human))+
  geom_bar(stat="identity", fill="#bdbdbd") + theme_cowplot()+
  scale_x_discrete(position="top")+
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100))+
  labs(y="%cells")+
  theme(axis.text.x  = element_text(size = 15, angle=45, hjust=0),
        axis.title.x = element_blank(), axis.line.x=element_blank(),
        axis.text.y  = element_text(size=15),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size=20),
        legend.position = "none")


# PROPORTION OF MOUSE MAIT CELLS IN EACH CLUSTER
order_mait <- c("MAIT0", "Cluster 7", "CyclingG2M", "CyclingS", "MAIT1", "MAIT17a", "MAIT17b")
mait.bpY <- ggplot(data=mtn.mait %>% select(mouse,propcells_mouse) %>% distinct(),
                   aes(x=factor(mouse, levels=rev(order_mait)), y=propcells_mouse))+
  geom_bar(stat="identity", fill="#bdbdbd") +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100))+
  scale_x_discrete(position="top") +
  labs(y="%cells")+ coord_flip() + theme_cowplot()+
  theme(axis.text    = element_text(size=15),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size=20),
        axis.title.y = element_blank(),
        axis.line.y  = element_blank(),
        legend.position = "none")



# BUBBLE PLOT
# library(scales)
mait.hm <- ggplot(mtn.mait, aes(x=factor(human, levels=paste0("MAIT_c", 0:6)),
                                y=factor(mouse, levels=rev(order_mait)))) +
  geom_point(aes(size=auroc, color= auroc))+
  geom_text(data=mtn.mait %>% filter(auroc>0.65) %>% mutate(across("auroc", \(x) round(x,2))), aes(label=auroc), color="white")+
  scale_size_continuous(limits=c(0,1), breaks=seq(0.2,0.8, by=0.2), range = c(1, 15))+
  scale_colour_gradient2(low="#d9d9d9", mid="white", high="#a50f15", midpoint=0.5, limits=c(0,1), name="AUROC", breaks=seq(0,1, by=0.2))+
  labs(x="Human clusters",y="Mouse clusters", size="AUROC")+
  theme_cowplot()+
  theme(legend.position  = "bottom",
        legend.key.width = unit(0.8, 'cm'),
        axis.text   = element_text(size=15),
        axis.title  = element_text(size=20),
        axis.text.x = element_text(angle=45, hjust=1))


# COMBINE
mait.tot <- (mait.bpX+plot_spacer() + plot_layout(widths = c(5, 1))) / (mait.hm + mait.bpY + plot_layout(widths = c(5, 1))) + plot_layout(heights = c(1, 5))
# ggsave("./data/cross-species/04_Metaneighbor_mait/mait_ms-hu_metaneighbor_bubbleplot3.svg", width=9, height=8)




## 4. GDT PLOT ####

# PROPORTION OF HUMAN GDT CELLS IN EACH CLUSTER
gdt.bpX <- ggplot(data=mtn.gdt %>% select(human,propcells_human) %>% distinct(),
                   aes(x=factor(human, levels=paste0("GD_c", 0:7)), y=propcells_human))+
  geom_bar(stat="identity", fill="#bdbdbd") + theme_cowplot()+
  scale_x_discrete(position="top")+
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100))+
  labs(y="%cells")+
  theme(axis.text.x  = element_text(size = 15, angle=45, hjust=0),
        axis.title.x = element_blank(),
        axis.line.x  = element_blank(),
        axis.text.y  = element_text(size=15),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size=20),
        legend.position = "none")


# PROPORTION OF MOUSE GDT CELLS IN EACH CLUSTER
order_gdt <- c("cKIT+ DN1", "DN2", "DN3", "Pre-selected GD", "Post-selected GD",
               "Pan GD (mainly CD24+)", "CD122+ GD", "CD24- GD")
gdt.bpY <- ggplot(data=mtn.gdt%>% select(mouse,propcells_mouse) %>% distinct(),
                   aes(x=factor(mouse, levels=rev(order_gdt)), y=propcells_mouse))+
  geom_bar(stat="identity", fill="#bdbdbd") +
  scale_x_discrete(position="top") +
  scale_y_continuous(limits=c(0,100), breaks=c(0,50,100))+
  labs(y="%cells")+ coord_flip() + theme_cowplot()+
  theme(axis.text    = element_text(size=15),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size=20),
        axis.title.y = element_blank(),
        axis.line.y  = element_blank(),
        legend.position = "none")


# BUBBLE PLOT
gdt.hm <- ggplot(mtn.gdt, aes(x=factor(human, levels=paste0("GD_c", 0:7)),
                                y=factor(mouse, levels=rev(order_gdt)))) +
  geom_point(aes(size=auroc, color=auroc))+
  geom_text(data=mtn.gdt %>% filter(auroc>0.65) %>% mutate(across("auroc", \(x) round(x,2))), aes(label=auroc), color="white")+
  scale_size_continuous(limits=c(0,1), breaks=seq(0.2,0.8, by=0.2), range = c(1, 15))+
  scale_color_gradient2(low="#d9d9d9", mid="white", high="#a50f15", midpoint=0.5, limits=c(0,1), name="AUROC", breaks=seq(0,1, by=0.2))+
  labs(x="Human clusters", y="Mouse clusters", size="AUROC")+
  theme_cowplot()+
  theme(legend.position  = "bottom",
        legend.key.width = unit(0.8, 'cm'),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        axis.text.y = element_text(size=15))


# COMBINE
gdt.tot <- (gdt.bpX+plot_spacer() + plot_layout(widths = c(5, 1))) / (gdt.hm + gdt.bpY + plot_layout(widths = c(5, 1))) + plot_layout(heights = c(1, 5))
# ggsave("./data/cross-species/04_Metaneighbor_gdt/gdt_ms-hu_metaneighbor_bubbleplot.jpeg", width=12, height=9)




## 3. COMBINE EVERYTHING ####
inkt.tot | mait.tot | gdt.tot
ggsave("./data/cross-species/nktmaitgdt_ms-hu_metaneighbor_bubbleplot3.pdf", width=27, height=8)
