# Unraveling the Phenotypic States of Human innate-like T Cells: Comparative Insights with Conventional T Cells and Mouse Models

The thymus provides a unique and necessary environment for the proper development of T cells. While T cell development is deeply characterized in mouse models, there is still limited understanding of the development of T cells in the human thymus, more particularly for rarer T cell lineages such as iNKT, MAIT and &gamma;&delta; T cells, commonly referred to as innate-like T cells. Here, we used a combination of fluorescence-associated cell sorting, single-cell genomics and flow cytometry strategies in order to better understand commonalities and differences in the development of conventional and innate T cells in the human thymus, in addition to characterizing their molecular phenotype in the periphery.

This study is available as a preprint on [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.12.07.570707v1). To explore the data interactively, you can use our [ShinyCell app](http://xspeciestcells.cshl.edu/).

This repository contains all the scripts used to generate the figures contained in our manuscript. The repository to generate the interactive shiny app can be found [here](https://github.com/meyer-lab-cshl/xspeciestcells-shiny).


## Software requirements

Data analysis and visualization was done in R version 4.1.3 and python version 3.8. Package versions were detailed in the methods section of the preprint. When necessary, specific scripts were run on an HPC cluster.

## Directory structure

- [`scripts`](./scripts/) contains the final scripts that allow the (1) preprocessing of the data and (2) generation of figures in the manuscript;
- [`figures`](./figures) contains the final figures present in the manuscript as PDF files, upon peer-reviewed publication.
- [`data_geo`](./data_geo) was created for the user to put the data files deposited on GEO (ADD LINK HERE);
- [`data_github`](./data_github) contains small data files which were not uploaded on GEO (cNMF output files, gene signatures from the literature, etc.).
