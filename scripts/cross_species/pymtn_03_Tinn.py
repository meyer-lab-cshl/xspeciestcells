# %%
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import pymn

# %%
# Import data
print("-- IMPORT ANNDATA --")
adata = sc.read_h5ad('/grid/meyer/home/scarcy/HumanThymus/crosspecies-innateT-metaneighbor/input/innateT_ms_hu_full_seu_gene_names.h5ad')   
print(adata)

# %%
# Make sure observation columns of interest are strings
adata.obs['study'] = adata.obs['study'].astype(str)
adata.obs['cell_annot'] = adata.obs['cell_annot'].astype(str)
adata.obs['cell.ident'] = adata.obs['cell.ident'].astype(str)
adata.obs = adata.obs.rename(columns={'orig.ident':'orig_ident', 'cell.ident':'cell_ident'})
print(adata.obs.head())

# %%
# Import the highly variable genes and add to anndata object
print("-- IMPORT HVG LIST --")
hvg = pd.read_csv(r'/grid/meyer/home/scarcy/HumanThymus/crosspecies-innateT-metaneighbor/input/list_hvg.csv', index_col=0)
adata.var = pd.merge(adata.var, hvg, on=['features'], how='inner')
adata.var.index = adata.var.index.astype(str) # overcome error message
print(adata.var.head())
print(adata.var.highly_variable.value_counts()) # should have 4229 HVGs

# %%
# Subset anndata object to only HVG
print("-- SUBSET ANNDATA TO ONLY HVGs --")
adatasub = adata[:,adata.var['highly_variable']]
print(adatasub) # 4229 variable features

# %%
# Run metaneighbor
print("-- RUN METANEIGHBOR --")
pymn.MetaNeighborUS(adatasub,
                    study_col='study',
                    ct_col='cell_annot',
                    fast_version=False)

# %%
print(adatasub.uns['MetaNeighborUS'].head())

# %%
# Save to CSV file
print("-- SAVE METANEIGHBOR TO CSV --")
adatasub.uns['MetaNeighborUS'].to_csv("/grid/meyer/home/scarcy/HumanThymus/crosspecies-innateT-metaneighbor/output/pymtn_crosspecies_innateT_slowversion_2023-11-13.csv")


