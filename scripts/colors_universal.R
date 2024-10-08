###
# Purpose: Define colors for every figure/plot
# Date: June 2023
# Author: Salomé Carcy
###


# Integrated UMAP clusters
cols_integrated <- c("0" = "#f4c40f", "1" = "#b75347", "2" = "#d8443c", "3" = "#e09351", "4" = "#2b9b81", 
                     "5" = "#421401", "6" = "#92c051", "7" = "#9f5691", "8" = "#17154f", "9" = "#74c8c3", 
                     "10" = "#5a97c1", "11" = "gold", "12" = "#a40000", "13" = "#72bcd5", "14" = "grey50",
                     "15" = "orange", "16" = "blueviolet", "17" = "#0a2e57")
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
                   "iNKT"  = "#9e9ac8")

# Cell states
cols_cellstate <- c("Tnaive"= "#b3e2cd",
                    "Tcm"   = "#f4cae4",
                    "Th17"  = "#cbd5e8",
                    "Temra" = "#fdcdac",
                    "Treg" = "#fbb4ae")

# GEP assignment
cols_GEPs <-  c("gep1"  =  "#f5bb50",
                "gep2"  =  "#ada43b",
                "gep3"  =  "#b0799a",
                "gep4"  =  "#f6b3b0",
                "gep5"  =  "#bf3729",
                "gep6"  =  "#17154f",
                "gep7"  =  "#355828",
                "gep8"  =  "#e48171",
                "gep9"  =  "#2f357c",
                "gep10" =  "#6c5d9e",  
                "gep11" =  "#e69b00")

# PBMC cluster colors by lineage
cols_thym_cd4  <- c(
  "CD4_thymus_c0" = "#DF6D27FF",
  "CD4_thymus_c1" = "#E9BE99FF",
  "CD4_thymus_c2" = "grey40",
  "CD4_thymus_c3" = "grey70",
  "CD4_thymus_c4" = "#a40000",
  "CD4_thymus_c5" = "gold",
  "CD4_thymus_c6" = "#72bcd5"
)
cols_thym_cd8  <- c(
  "CD8_thymus_c0" = "#DF6D27FF",
  "CD8_thymus_c1" = "grey40",
  "CD8_thymus_c2" = "#AB6969",
  "CD8_thymus_c3" = "#a40000",
  "CD8_thymus_c4" = "gold",
  "CD8_thymus_c5" = "#9f5691"
)
cols_thym_gdt  <- c(
  "GDT_thymus_c0" = "#d8443c",
  "GDT_thymus_c1" = "#e09351",
  "GDT_thymus_c2" = "gold",
  "GDT_thymus_c3" = "#9f5691",
  "GDT_thymus_c4" = "#72bcd5",
  "GDT_thymus_c5" = "blueviolet",
  "GDT_thymus_c6" = "olivedrab2",
  "GDT_thymus_c7" = "grey50"
)
cols_thym_nkt  <- c(
  "iNKT_thymus_c0" = "#d8443c",
  "iNKT_thymus_c1" = "#e09351",
  "iNKT_thymus_c2" = "gold",
  "iNKT_thymus_c3" = "#74c8c3",
  "iNKT_thymus_c4" = "#5a97c1",
  "iNKT_thymus_c5" = "#a40000",
  "iNKT_thymus_c6" = "#72bcd5"
)
cols_thym_mait <- c(
  "MAIT_thymus_c0" = "#d8443c",
  "MAIT_thymus_c1" = "#e09351",
  "MAIT_thymus_c2" = "gold",
  "MAIT_thymus_c3" = "#74c8c3",
  "MAIT_thymus_c4" = "#a40000",
  "MAIT_thymus_c5" = "#5a97c1",
  "MAIT_thymus_c6" = "orange"
)

# PBMC cluster colors by lineage
cols_pbmc_cd4  <- c("0" = "#DF6D27FF", "1" = "#E9BE99FF", "2" = "grey40", "3" = "#7EF547", "4" = "grey70", "5" = "#a40000")
cols_pbmc_cd8  <- c("0" = "#DF6D27FF", "1" = "#E9BE99FF", "2" = "grey70", "3" = "#a40000", "4" = "gold")
cols_pbmc_gdt  <- c("0" = "#DF6D27FF", "1" = "grey40", "2" = "#AB6969", "3" = "#a40000", "4" = "gold")
cols_pbmc_nkt  <- c("0" = "#DF6D27FF", "1" = "#5B8DB9FF", "2" = "grey40", "3" = "#a40000")
cols_pbmc_mait <- c("0" = "grey90", "1" = "grey70", "2" = "grey40", "3" = "#a40000")
