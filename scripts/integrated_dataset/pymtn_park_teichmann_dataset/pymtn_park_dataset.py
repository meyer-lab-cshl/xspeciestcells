# %%
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import pymn

# %%
# Import data
print("-- IMPORT ANNDATA --")
adata = sc.read_h5ad('/grid/meyer/home/scarcy/HumanThymus/park-metaneighbor/input/seuratobj_convert.h5ad')   
print(adata)

# %%
# Make sure observation columns of interest are strings
adata.obs['study'] = adata.obs['study'].astype(str)
adata.obs['cell_type'] = adata.obs['cell_type'].astype(str)
# adata.obs[["study", "cell_type"]] # should have 114363 cells


# %%
# Import the highly variable genes and add to anndata object
print("-- IMPORT HVG LIST --")
hvg = pd.read_csv(r'/grid/meyer/home/scarcy/HumanThymus/park-metaneighbor/input/list_hvg.csv', index_col=0)
adata.var = pd.merge(adata.var, hvg, on=['features'], how='inner')
adata.var.index = adata.var.index.astype(str) # overcome error message
print(adata.var.head())

# %%
# Subset anndata object to only HVG
# adata.var['features'] = adata.var['features'].astype(str)
# hvg_list = adata.var.loc[adata.var['highly_variable']==True][["features"]].values
# len(hvg_list) # 3106 variable features
print("-- SUBSET ANNDATA TO ONLY HVGs --")
adatasub = adata[:,adata.var['highly_variable']]
print(adatasub) # 3106 variable features

# %%
# Personal test on very few cells
# adatasub = adatasub[adatasub.obs['cell_type'].isin(["γδT", "c6", "NKT", "c17"])]
# adatasub

# Make sure the obs are strings
adatasub.obs['study'] = adatasub.obs['study'].astype(str)
adatasub.obs['cell_type'] = adatasub.obs['cell_type'].astype(str)

# %%
# Run metaneighbor
print("-- RUN METANEIGHBOR --")
pymn.MetaNeighborUS(adatasub,
                    study_col='study',
                    ct_col='cell_type',
                    fast_version=False)

# %%
#pymn.plotMetaNeighborUS(adatasub, figsize=(10, 10), cmap='coolwarm', fontsize=10)
print(adatasub.uns['MetaNeighborUS'].head())

# %%
# Save to CSV file
print("-- SAVE METANEIGHBOR TO CSV --")
adatasub.uns['MetaNeighborUS'].to_csv("/grid/meyer/home/scarcy/HumanThymus/park-metaneighbor/output/pymtn_park_slowversion_2023-07-11.csv")


