###
# Purpose: Define colors for every figure/plot
# Date: June 2023
# Author: Salomé Carcy
###


# Integrated UMAP clusters
cols_integrated <- c("0" = "#f4c40f", "1" = "#b75347", "2" = "#d8443c", "3" = "#e09351", "4" = "#2b9b81", 
                     "5" = "#421401", "6" = "#92c051", "7" = "#9f5691", "8" = "#17154f", "9" = "#74c8c3", 
                     "10" = "#5a97c1", "11" = "gold", "12" = "#a40000", "13" = "#72bcd5", "14" = "grey50",
                     "15" = "orange", "16" = "blueviolet", "17" = "#0a2e57", "18" = "#bdbdbd")
# Park clusters
cols_park <- c("DN(early)" = "#78c679",
               "DN(P)" = "#41ab5d",
               "DN(Q)" = "#238443",
               "γδT" = "#92c051",
               "DP(P)" = "#b75347",
               "DP(Q)" = "#d8443c",
               "αβT(entry)" = "#e09351",
               "CD8+T"= "#5a97c1",
               "CD8αα(I)" = "#421401",
               "CD8αα(II)" = "#0a2e57",
               "CD4+T"= "gold",
               "T(agonist)" = "#9f5691",
               "Treg(diff)" = "#9f5691",
               "Treg" = "blueviolet",
               "Th17" = "#a40000",
               "NKT" = "#72bcd5")

# Cell lineages
cols_lineages <- c("CD4"  = "#74c476",
                   "CD8"  = "#df65b0",
                   "GD"   = "#08519c",
                   "MAIT" = "#9ecae1",
                   "NKT"  = "#9e9ac8")

# Cell states
cols_cellstate <- c("Tnaive"= "#b3e2cd",
                    "Tcm"   = "#f4cae4",
                    "Th17"  = "#cbd5e8",
                    "Temra" = "#fdcdac",
                    "Treg" = "#fbb4ae")

# PBMC cluster colors by lineage
cols_pbmc_cd4  <- c("0" = "#DF6D27FF", "1" = "#E9BE99FF", "2" = "grey40", "3" = "#7EF547", "4" = "grey70", "5" = "#a40000")
cols_pbmc_cd8  <- c("0" = "#DF6D27FF", "1" = "#E9BE99FF", "2" = "grey40", "3" = "grey70", "4" = "#a40000", "5" = "gold")
cols_pbmc_gdt  <- c("0" = "#DF6D27FF", "1" = "grey40", "2" = "#AB6969", "3" = "#a40000", "4" = "gold")
cols_pbmc_nkt  <- c("0" = "#DF6D27FF", "1" = "#5B8DB9FF", "2" = "grey40", "3" = "#a40000")
cols_pbmc_mait <- c("0" = "grey90", "1" = "grey70", "2" = "grey40", "3" = "#a40000")
