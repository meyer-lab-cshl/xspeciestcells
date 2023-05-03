# Purpose: run UMAP with different parameters
# Author: Salomé Carcy
# Date: March 2023


# **************
# 1. IMPORT ####
# **************

# Import librairies
library(Seurat)
library(cowplot)
library(ggplot2)

path <- "~/Projects/HumanThymusProject/data/human-thymus/HumanData_18_UMAPparam/"

# Import data
seur <- readRDS("~/Projects/HumanThymusProject/data/raw_data/human_data/seurat_filtered_harmony_02_15_23.RDS") # removed low-quality cells

cols_integrated <- c("0" = "#f4c40f", "1" = "#b75347", "2" = "#d8443c", "3" = "#e09351", "4" = "#2b9b81", 
                     "5" = "#421401", "6" = "#92c051", "7" = "#9f5691", "8" = "#17154f", "9" = "#74c8c3", 
                     "10" = "#5a97c1", "11" = "gold", "12" = "#a40000", "13" = "#72bcd5", "14" = "grey50",
                     "15" = "orange", "16" = "blueviolet", "17" = "#0a2e57", "18" = "#bdbdbd")
DimPlot(seur, reduction="UMAP_50", group.by = "new_clusters", label=T, repel=T, cols=cols_integrated)+
  labs(title=paste0("n.neighbors=30  |  min.dist=0.3"))+
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        legend.position="none",
        panel.border=element_rect(colour="black", fill=NA, linewidth=1),
        axis.title=element_text(hjust=0))
# ggsave("~/Projects/HumanThymusProject/umapp/umap_n30_mindist0.3.jpeg", width=6, height=6)



# *****************
# 2. RUN UMAPs ####
# *****************

ndim <- 20 # nb of PCs chosen by Laurent
nneighbors <- c(10, 20, 30, 50, 100)
mindist <- c(0.1, 0.3, 0.5, 0.7)

# plist <- list()

for(neighb in nneighbors){
  # plist[[paste0("neighb_", neighb)]] <- list()
  
  for(dist in mindist){
    print(paste0("n.neighbors=", neighb, " | min.dist=", dist))
    seurobj <- seur
    seurobj@reductions$initial_umap <- NULL
    seurobj@reductions$UMAP_50 <- NULL
    seurobj <- RunUMAP(seur, reduction = "harmony", dims = 1:ndim, seed.use=42, reduction.name = "newumap", reduction.key="UMAP",
                       n.neighbors=neighb,
                       min.dist=dist)
    p <- DimPlot(seurobj, reduction="newumap", group.by = "new_clusters", label=T, repel=T, cols=cols_integrated)+
      labs(title=paste0("n.neighbors=", neighb, "  |  min.dist=", dist))+
      theme(axis.ticks=element_blank(),
            axis.text=element_blank(),
            legend.position="none",
            panel.border=element_rect(colour="black", fill=NA, linewidth=1),
            axis.title=element_text(hjust=0))
    filename <- paste0("~/Projects/HumanThymusProject/umapp/umap_n", neighb, "_mindist", dist, ".jpeg")
    ggsave(plot=p, filename=filename, width=6, height=6)
    # plist[[paste0("neighb_", neighb)]][[paste0("mindist_", dist)]] <- p
    # plist[[paste0("neighb_", neighb)]][[paste0("mindist_", dist)]] <- paste0("n.neighbors=", neighb, "  |  min.dist=", dist)
    
  }
}




# ******************
# 3. CREATE GIF ####
# ******************

## list file names and read in
imgs <- list.files("~/Projects/HumanThymusProject/umapp/mindist", full.names = TRUE)
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 1 frame per second
img_animated <- image_animate(img_joined, fps = 1)

## view animated image
img_animated

## save to disk
image_write(image = img_animated,
            path = "~/Projects/HumanThymusProject/umapp/mindist.gif")


