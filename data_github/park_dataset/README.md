If you would like to reproduce Fig 7 or Supplementary Fig 4, you will need to download the [Thymus Atlas from Park et al.](https://www.science.org/doi/10.1126/science.aay3224). Here are the files you should download (which we have gitignored), and the corresponding links:

- full data of Human Thymus Atlas: on [cellxgene website](https://cellxgene.cziscience.com/collections/de13e3e2-23b6-40ed-a413-e9e12d7d3910), download the `.h5ad` file containing 255,901 cells (named `c6e08ab6-ab3b-41dc-8058-8e6442e081ec.h5ad`).
- thymocytes only data of Human Thymus Atlas: on [cellxgene website](https://cellxgene.cziscience.com/collections/de13e3e2-23b6-40ed-a413-e9e12d7d3910), download the dataset containing 76,994 cells, download **both the `.h5ad` and the `.rds` files** (named `de8665a2-0476-4865-b4af-c7b8d3b1b87f.h5ad` and `de8665a2-0476-4865-b4af-c7b8d3b1b87f.rds`);
- `h5ad` file of Human Thymus Atlas: [HTA08_v01_A05_Science_human_fig1.h5ad](https://app.cellatlas.io/thymus-development/dataset/8/scatterplot), from which we will extract the metadata with the correct cell annotations (the `h5ad` file from CellxGene doesn't contain cell annotation corresponding to Figure 1C in the original paper).
- `h5ad` file of Human Thymus Atlas: [HTA08_v01_A06_Science_human_tcells.h5ad](https://app.cellatlas.io/thymus-development/dataset/10/scatterplot), from which we will extract the metadata with the correct cell annotations (the `h5ad` file from CellxGene doesn't contain cell annotation corresponding to Figure 2C in the original paper).
- `h5ad` file of Murine Thymus Atlas: [HTA08_v02_A04_Science_mouse_total.h5ad](https://zenodo.org/records/5500511), which we will directly use to look at the expression pattern of _Cd1d1, Slamf1, Slamf6_ in the murine thymus (see Fig 7 in our manuscript).

<br/>

When running the script [download_park_thymic_atlas.Rmd](../../scripts/cross_species/download_park_thymic_atlas.Rmd), it will generate three seurat objects:

- `park_seurat_human.rds`
- `park_seurat_mouse.rds`
- `park_seurat_human_thymocytes.rds`

These seurat objects are necessary to reproduce Fig 7 and Supplementary Fig 4.