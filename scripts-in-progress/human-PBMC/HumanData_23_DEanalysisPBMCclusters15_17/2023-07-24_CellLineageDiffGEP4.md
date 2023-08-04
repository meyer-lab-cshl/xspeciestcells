Title: DE genes btw lineages clusters 15-17
Date: 2023-07-24
Category: Human project
Summary: DE analysis on PBMCs to detect lineage-specific signature in clusters 15-17
Tags: scRNAseq, effector programs

[TOC]

# 1. Context
From clustering, cNMF, DE analysis we observe that cell state transcriptomic signal is stronger than lineage.
Laurent's lab exposed GDT cells or MAIT cells to different combinations of cytokines:
- IL12+IL18
- IL15
- IL15+IL17
We would expect a small proportion of GDT, and most of MAITs to react (cells in GEP1 express receptors to these cytokines). However, only GDT proliferated in presence of these cytokines, so clearly there must be still some differences btw MAIT and GDTs in the "same" cell state (or what we think is the same cell state).

We have been investigating those differences in [2023-07-19_CellLineageDiffGEP1.md](../HumanData_23_DEanalysisPBMCclusters13_14/2023-07-19_CellLineageDiffGEP1.md), but now we want to repeat it in clusters 15-17 (GEP4).


<br/>

# 2. DE analysis

## 2.1. LRT test on all PBMCs in clusters 15-17, excluding batches G-H
First I performed a LRT test on all PBMCs present in clusters 15-17 and batches E-F-I (avoiding batches G-H where we don't know for sure what is CD8 and what is MAIT).

- grouping cells by batch/donor, lineage, and keeping only groups with **50 cells** or more (threshold at 100 cells min removes a MAIT_E group)
- full model: `batchdonor_id + lineage_id`
- reduced model: `~ batchdonor_id`
- threshold: padj < 0.05

Heatmap:
<img src="All_heatmap_clust15-17_batchEFI_padj0_05.pdf" height="400"/>

After a bit of thinking, I thought that it was more thorough to do this analysis on batches/donors E5 and I11, both batches where we have all 5 lineages sequenced. Because in the previous heatmap we could see XIST as a DE gene, which is just because MAIT_E and MAIT_F are female donors, while batch I contains mostly male donors.

<br/>

## 2.2. LRT test on all PBMCs in clusters 15-17, including batches E5 and I11
E5 and I11 are the only batch/donor pairs where all 5 lineages were sequenced.

- grouping cells by batch/donor, lineage, and keeping only groups with **10 cells** or more (threshold at 50 cells min removes several groups...)
- full model: `batchdonor_id + lineage_id`
- reduced model: `~ batchdonor_id`
- threshold: padj < 0.05
<img src="All_heatmap_clust15-17_batchEI_padj0_05.pdf" height="400"/>
