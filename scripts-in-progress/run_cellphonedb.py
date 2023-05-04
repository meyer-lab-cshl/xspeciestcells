
import pandas as pd
import sys
import os

import anndata
from cellphonedb.src.core.methods import cpdb_degs_analysis_method

## required files
#cpdb_file_path: (mandatory) path to the database.
#meta_file_path: tsv with cell barcodes (1 col, barcode_sample) and cluster labels (2 col, cell_type)
#counts_file_path:  normalized counts file (not z-transformed),
#   either in text format or h5ad (recommended).
#degs_file_path: tsv with clusters (1 col, cluster) and their DEG/Markers (2 col, gene)

basedir = '/Users/hmeyer/projects/20220809_Thymic-iNKT-CrossSpecies'
out_path = os.path.join(basedir, 'results/cellphonedb_tec-inkt-dp')

cpdb_file_path = os.path.join(basedir, 'cellphonedb/db/v4.1.0/cellphonedb.zip')
meta_file_path = os.path.join(basedir, 'data/human-thymus/HumanData_16_ParkIntegration/metadata_Tstromal.tsv')
counts_file_path = os.path.join(basedir, 'data/human-thymus/HumanData_16_ParkIntegration/park_garpin_full_seu_gene_names_updated_anno.h5ad')
degs_file_path = os.path.join(basedir, 'data/human-thymus/HumanData_16_ParkIntegration/markers_Tstromal.tsv')


deconvoluted, means, relevant_interactions, significant_means = cpdb_degs_analysis_method.call(
    cpdb_file_path = cpdb_file_path,
    meta_file_path = meta_file_path,
    counts_file_path = counts_file_path,
    degs_file_path = degs_file_path,
    counts_data = 'hgnc_symbol',
    # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threshold = 0.1,
    # Sets the rounding for the mean values in significan_means.
    result_precision = 3,
    # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    separator = '|',
    output_path = out_path                                     # Path to save results
    )
