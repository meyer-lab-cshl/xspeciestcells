Title: Cell lineage signature
Date: 2023-07-21
Category: Human project
Summary: DE analysis on PBMCs to detect lineage-specific signature
Tags: scRNAseq, effector programs

[TOC]

# 1. Context
From clustering, cNMF, DE analysis we observe that cell state transcriptomic signal is stronger than lineage. But are there some lineage-specific signatures?


# 2. DE analysis

## 2.1. LRT test on all PBMCs
First I performed a LRT test on all PBMCs

- grouping cells by batch/donor, cluster, lineage and keeping only groups with 50 cells or more
- full model: `batchdonor_id + cluster_id + lineage_id`
- reduced model: `~ batchdonor_id + cluster_id`
- threshold: padj < 0.01

Heatmap:
<img src="allPBMCs_allbatches_heatmap_padj0_01.pdf" height="400"/>

There are multiple problems with this approach though:
1. some CD8s in batches G and H are MAITs in disguise, because of the sorting strategy;
2. every batch has a different combination of lineages, so that means any DE genes found between lineages, we don't know if it's a _technical_ or _biological_ signal;


## 2.2. LRT test on PBMCs in batches E and I
I repeated the analysis by restricting to cells in batches E and I (donors 5 and 11), the only ones where all lineages are present at once, and the same sorting strategy was applied.

Heatmap:
<img src="allPBMCs_batchE-I_heatmap_padj0_01.pdf" height="400"/>

There is still another limitation to this approach: the LRT test decides to take one lineage as reference (in my case, CD4), whereas I'm interested in genes that are specifically upregulated in CD4, or CD8, etc.
```R
> resultsNames(dds)
[1] "Intercept"               "batchdonor_id_I11_vs_E5" "lineage_id_CD8_vs_CD4"
[4] "lineage_id_GD_vs_CD4"    "lineage_id_MAIT_vs_CD4"  "lineage_id_NKT_vs_CD4"
[7] "cluster_id_11_vs_10"     "cluster_id_12_vs_10"     "cluster_id_13_vs_10"
[10] "cluster_id_14_vs_10"     "cluster_id_15_vs_10"     "cluster_id_16_vs_10" 
[13] "cluster_id_17_vs_10"     "cluster_id_3_vs_10"      "cluster_id_6_vs_10"
[16] "cluster_id_7_vs_10"      "cluster_id_9_vs_10"     
```

We could overcome that by looking at the clusters of genes that exhibit a specific pattern (of being upregulated in some lineages versus others):
<img src="allPBMCs_batchE-I_DEgenescluster_padj0_05.pdf" height="400"/>


## 2.3. LRT test on PBMCs in batches E and I, pairwise contrasts
I repeated the analysis by testing all pairwise contrasts. I built a function where you can specify one lineage (e.g. CD4) and it will go through a loop to return a list of upregulated genes for each pairwise contrast (e.g. CD4vsCD8, CD4vsMAIT, CD4vsNKT, CD4vsGD).

Observations:
- had to be very lax on thresholds (padj<0.05 and log2FC>0.1) otherwise it was hard to get DE genes between NKT and others
- very few common upregulated genes, so had to keep genes that were upregulated in contrast to _at least_ 2 other lineages (e.g. genes upregulated in CD4vsCD8 and CD4vsMAIT, but not in CD4vsNKT)
- only "true" lineage-specific genes were TRDC for GDs and CD8B for CD8s

Heatmap:
<img src="allPBMCs_batchE-I_lineagespecific_upregulatedgenes_heatmap_padj0_05.pdf" height="400"/>

Hypotheses:
- not enough resolution with single cell (genes not well detected) to find lineage-specific genes (that are expressed in 1 lineage and none others)
- or maybe there are no "lineage-specific" genes, and there are only lineage-specific signatures (unique combination of genes specific to one lineage, but there may be partial overlap with the gene signature of another lineage)

Next steps:
- incorporate batches C and D from thymus (also have all lineages present) to have more power?
- incorporate batch F to include more cells (even though experimental design is not great)?
- try variance partitioning


## 2.4. LRT test on PBMCs in batches E-F-I, pairwise contrasts
Parameters:
- padj < 0.05
- logFC > 0.1
- for each lineage, gene upregulated compared to at least 3 other lineages (except NKT, at least upregulated in 2 contrasts, because very few genes DE)

Heatmap:
<img src="allPBMCs_batchE-F-I_lineagespecific_upregulatedgenes_heatmap_padj0_05.pdf" height="400"/>

