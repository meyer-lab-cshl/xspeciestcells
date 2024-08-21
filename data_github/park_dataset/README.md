If you would like to reproduce Fig 7 or Supplementary Fig 4, you will need to download the [Thymus Atlas from Park et al.](https://www.science.org/doi/10.1126/science.aay3224). Here are the files you should download, and the corresponding links:

- `h5ad` file of Human Thymus Atlas: [c6e08ab6-ab3b-41dc-8058-8e6442e081ec.h5ad](https://cellxgene.cziscience.com/collections/de13e3e2-23b6-40ed-a413-e9e12d7d3910), download the one containing 255,901 cells.
- `h5ad` file of Human Thymus Atlas: [HTA08_v01_A05_Science_human_fig1.h5ad](https://app.cellatlas.io/thymus-development/dataset/8/scatterplot), from which we will extract the metadata with the correct cell annotations (the `h5ad` file from CellxGene doesn't contain cell annotation corresponding to the figures in the original paper).

<br/>

When running the script [download_park_thymic_atlas.Rmd](../../scripts/cross_species/download_park_thymic_atlas.Rmd), it will generate two seurat objects:

- `park_seurat_human.rds`
- `park_seurat_mouse.rds`

These seurat objects are necessary to reproduce Fig 7 and Supplementary Fig 4.