Title: Gene lists overlap
Date: 2023-07-22
Category: Human project
Summary: Compare gene lists btw GEPs, Cano-Gamez, Rose, Poon
Tags: scRNAseq, effector programs

[TOC]

# 1. Context
From our analyses, we believe that PBMCs are into one of 3 effector programs:
- GEP1 (Th1/Th17)
- GEP4 (Temra)
- GEP6 (Tcm)
The question is whether these are new programs or do they match/overlap with other publications?


# 2. Gene lists prep
**GEP genes**
- provided by Laurent
- btw 500-4000 genes per GEP

**Cano-Gamez gene lists**
- from scRNAseq on CD4 PBMCs
- https://www.nature.com/articles/s41467-020-15543-y
- selected genes with padj<0.05 and log2FC>0 (min log2FC is ~0.25)
- btw 17-183 genes per program

**Rose gene lists**
- from bulk RNAseq on CD4/CD8 PBMCs
- https://www.nature.com/articles/s42003-023-04747-9
- genes provided in supptable already have padj<0.05
- btw 70-363 genes per program

**Poon gene lists**
- from scRNAseq on T cells across multiple tissues (including PBMCs)
- https://www.nature.com/articles/s41590-022-01395-9
- selected top 500 marker genes of each cluster (top marker genes with lowest padj)

<br/>

# 3. % GEP genes found in other gene programs

## 3.1. With top 500 marker genes from Poon dataset
Initial analysis: selected top 500 marker genes for each cluster (with lowest padj).
Generated a null hypothesis by selecting randomly an equal number of genes as the size of the GEP (e.g. for GEP1, selecting randomly 504 genes in same expression bins). Creating 1000 random gene sets of equal size of the GEP, and checking their overlap with each gene program. P-value is calculated by #random sets with higher overlap than the observed overlap divided by 1000. Adjusted p-value is pval multiplied by nb of geneprograms (nb of comparisons), so 24.

Barplot with top 10 overlaps:
<img src="geneoverlapwithGEP_bars_top10_pval.jpeg" height="400"/>

As correctly mentioned by Laurent, the max overlap we get is ~35%. So that means in best case scenario we have only a third of GEP1 genes that are found in Poon_CD8MAIT gene list... Quite low (what are the other two thirds of genes doing?...). **However** it's a freaking bias because the GEP gene lists are much longer than the other gene programs: so even if all of the Poon_CD8MAIT genes were found in GEP1, we would only find a portion of GEP1 genes in Poon_CD8MAIT.

To note, GEP1 and Poon_CD8MAIT are same gene set size. However I selected only the top 500 marker genes with lowest padj in Poon, so maybe the 2/3 of GEP1 genes that we can't make overlap will be found in the rest of Poon_CD8MAIT marker genes...


## 3.2. With all marker genes from Poon dataset
Instead of cutting # marker genes per cluster/gene program in Poon dataset, I took all of marker genes with padj<0.01 and log2FC>0.25 (to have similar threshold to Cano-Gamez). Now the Poon gene lists are between 2400-6500 genes per program (so mostly longer than the GEPs).

Barplot with top 10 overlaps:
<img src="geneoverlapwithGEP_bars_top10_fullPoon.jpeg" height="400"/>

We've got much better overlap now (60-75% with GEP1/GEP4). But now we're overrepresented by Poon gene lists (because they're much longer so they're gonna get better overlap with the GEPs). Let's try subseting to only "significant" (padj<0.05) overlaps:
<img src="geneoverlapwithGEP_bars_top10_fullPoon_sig.jpeg" height="200"/>

Still this is not the best analysis as most of the gene programs have too few genes compared to the size of the GEPs. We need to do the analysis the other way around.


# 4. % genes in other gene programs found in GEPs

Had to adapt the `NullHypothesis` function so that it wouldn't consider the GEPs as reference, but could take in any gene programs as reference as input.

Anyways, a bit annoying because even though Cano-Gamez and Rose have shorter gene lists than GEPs, that's not the case for Poon (>2400 genes per program). The "PBMC GEPs" have between 504-1700 genes per gene program of course. So I split the analysis in 2:

1. **Proportion of Cano-Gamez or Rose genes found in GEPs**
<img src="geneoverlaptoGEP_bars_all_noPoon.jpeg" height="400"/>

2. **Proportion of GEP genes found in Poon gene sets**
<img src="geneoverlaptoGEP_bars_all_fullPoon_Poononly.jpeg" height="400"/>

That's so annoying to think about. Let's try again selecting in Poon the top 500  (more precisely 504) genes with lowest padj, so that we can compute the %Poon genes found in GEPs:
<img src="geneoverlaptoGEP_bars_all_Poontop500.jpeg" height="400"/>

