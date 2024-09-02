## RUN METANEIGHBOR
## author: Salomé Carcy
## date: 02 September 2024
## purpose: script to run on a HPC to run pyMetaNeighbor between our thymocytes dataset and the Park et al. Human Thymus Atlas

# %%
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import pymn

# %%
# Import data
print("-- IMPORT ANNDATA --")
adata = sc.read_h5ad('./input/human_thym_merged_with_park_240902.h5ad')   
print(adata)

# %%
# Make sure observation columns of interest are strings
adata.obs['study'] = adata.obs['study'].astype(str)
adata.obs['clusters_to_compare'] = adata.obs['clusters_to_compare'].astype(str)
print(adata.obs.head())


# %%
# Import the highly variable genes and add to anndata object
print("-- IMPORT HVG LIST --")
hvg = pd.read_csv(r'./input/human_thym_merged_with_park_list_hvg_240902.csv', index_col=0)
hvg.head()
# add HVGs to anndata object
adata.var = pd.merge(adata.var, hvg, on=['features'], how='inner')
adata.var.index.name = 'features'
adata.var.index = adata.var.index.astype(str) # overcome error message
print(adata.var.head())
print(adata.var.highly_variable.value_counts()) # should have 3104 HVGs

# %%
# Subset anndata object to only HVG
print("-- SUBSET ANNDATA TO ONLY HVGs --")
adatasub = adata[:,adata.var['highly_variable']]
print(adatasub) # 3104 features

# %%
# Run metaneighbor
print("-- RUN METANEIGHBOR --")
pymn.MetaNeighborUS(
    adatasub,
    study_col='study',
    ct_col='clusters_to_compare',
    fast_version=False
)

# %%
# quick look at output
print(adatasub.uns['MetaNeighborUS'].head())

# %%
# Save to CSV file
print("-- SAVE METANEIGHBOR TO CSV --")
adatasub.uns['MetaNeighborUS'].to_csv("./output/human_thym_merged_with_park_pymtn_slow_version_output_240902.csv")