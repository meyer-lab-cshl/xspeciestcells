This directory contains data related to analyses performed with cNMF.

More specifically:
- `limited_nonimputed_genes_per_gep_post_rank_threshold_k12.csv` contains the list of top genes associated with each GEP, this file is necessary to run the script [compare_geps_with_literature_signatures.Rmd](../../scripts/cNMF/compare_geps_with_literature_signatures.Rmd);
- `non_imputed_cNMF_allcells.usages.k_12.dt_0_02.consensus.txt` contains the usage of each GEP in each cell, this file is necessary to run the scripts [plot_gep12_batcheffect.Rmd](../../scripts/cNMF/plot_gep12_batcheffect.Rmd), [plot_GEPs_on_thymic_umaps.Rmd](../../scripts/cNMF/plot_GEPs_on_thymic_umaps.Rmd), and [ribbonplot_pbmc_Tlineages_GEPusage.Rmd](../../scripts/cNMF/ribbonplot_pbmc_Tlineages_GEPusage.Rmd). These GEP usages were also added to the integrated human seurat object deposited on GEO (GSE249684), and can be observed on our interactive [shinyapp](https://xspeciestcells.cshl.edu/).
