

# Download data
library(GEOquery)
gse = getGEO("GSE126030")[[1]]

# Transform it into a SummarizedExperiment object
library(SummarizedExperiment)
se <- as(gse, "SummarizedExperiment") # doesn't contain any expression
colData(se)

# # Get the RNA files (thymus, spleen, lymph)
# getGEOSuppFiles("GSM5417469", baseDir="~/Projects/HumanThymusProject/data/cross-species/04_Metaneighbor_gdt/li2023", makeDirectory=FALSE)
# getGEOSuppFiles("GSM5417470", baseDir="~/Projects/HumanThymusProject/data/cross-species/04_Metaneighbor_gdt/li2023", makeDirectory=FALSE)
# getGEOSuppFiles("GSM5417471", baseDir="~/Projects/HumanThymusProject/data/cross-species/04_Metaneighbor_gdt/li2023", makeDirectory=FALSE)
# 
# 
# # Import data into Seurat
# library(Seurat)
# seu.thymus <- Read10X(data.dir="./data/cross-species/04_Metaneighbor_gdt/li2023/Thymus_RNA/")

