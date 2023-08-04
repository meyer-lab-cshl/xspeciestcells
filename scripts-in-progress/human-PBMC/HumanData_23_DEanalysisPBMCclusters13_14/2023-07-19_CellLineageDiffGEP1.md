Title: DE genes btw lineages clusters 13-14
Date: 2023-07-21
Category: Human project
Summary: DE analysis on PBMCs to detect lineage-specific signature in clusters 13-14
Tags: scRNAseq, effector programs

[TOC]

# 1. Context
From clustering, cNMF, DE analysis we observe that cell state transcriptomic signal is stronger than lineage.
Laurent's lab exposed GDT cells or MAIT cells to different combinations of cytokines:
- IL12+IL18
- IL15
- IL15+IL17
We would expect a small proportion of GDT, and most of MAITs to react (cells in GEP1 express receptors to these cytokines). However, only GDT proliferated in presence of these cytokines, so clearly there must be still some differences btw MAIT and GDTs in the "same" cell state (or what we think is the same cell state).


<br/>

# 2. DE analysis

## 2.1. LRT test on all PBMCs in clusters 13-14
First I performed a LRT test on all PBMCs present in clusters 13-14.

- grouping cells by batch/donor, lineage, and keeping only groups with 100 cells or more
- full model: `batchdonor_id + lineage_id`
- reduced model: `~ batchdonor_id`
- threshold: padj < 0.01

Heatmap:
<img src="All_heatmap_clust13-14_padj0_01.pdf" height="400"/>
<img src="il12rb2_gdmait.jpeg" height="200"/>

We can see that IL12RB2 is significantly ⬆️ in GDTs compared to MAITs. No other cytokine receptor of IL12,15,17,18 seems DE.

Heatmap on cytokine receptors of interest specifically (careful, this doesn't prove any significant DE):
<img src="All_heatmap_clust13-14_cytokinereceptors.pdf" height="400"/>


The limitation with the LRT is that the genes deemed "significant" and genes that are better modelled when `lineage_id` is included in the model. So we don't know if it's information about MAIT/GD or any other lineage that's useful to predicting the genes' expression level. We can look at patterns with `DEGreport::degPatterns()` though:

<img src="degpatterns_clust13-14_padj0_01.jpeg" height="400"/>

I saved this list of genes & their "group" assignment as a `.csv` file. A lot of cytotoxic genes are found in groups 5-6:

<img src="degpatterns_clust13-14_padj0_01_goenrich_groups5-6.jpeg" height="200"/>


**Big limitation**: here I included all batches, but in batches G and H we know that all CD8+ cells were sorted, thus part of them must be MAITs...

<br/>

## 2.2. LRT test on all PBMCs in clusters 13-14, including batches E5 and I11
- grouping cells by batch/donor, lineage, and keeping only groups with 10 cells or more (to get some CD4 and CD8 groups)
- full model: `batchdonor_id + lineage_id`
- reduced model: `~ batchdonor_id`
- threshold: padj < 0.05

<img src="All_heatmap_clust13-14_batchEI_padj0_05.pdf" height="400"/>

About half of the genes (28/68) are the same as found in [2.1.]

<br/>

## 2.3. LRT test on MAITs/GDTs in clusters 13-14
- grouping cells by batch/donor, lineage, keeping only MAITs/GDTs, and only groups with 100 cells or more
- full model: `batchdonor_id + lineage_id`
- reduced model: `~ batchdonor_id`
- threshold: padj < 0.05

<img src="MAITvsGDT_heatmap_clust13-14_padj0_05.pdf" height="400"/>