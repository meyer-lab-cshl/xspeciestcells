from pyscenic.binarization import binarize
import pandas
import seaborn

auc_mtx_inkt_loc = "/home/tonya/Downloads/thymus.nkt_02_28_23_raw_counts_multirun_03162023/full_object_thymus.nkt_02_28_23_raw_counts_matrix_consensus_10runs_auc_mtx.csv"
auc_mtx = pandas.read_csv(auc_mtx_inkt_loc, index_col= 0)
binarized = binarize(auc_mtx, num_workers = 2)
binarized[0].to_csv("/home/tonya/Downloads/thymus.nkt_02_28_23_raw_counts_multirun_03162023/full_object_thymus.nkt_02_28_23_raw_counts_matrix_consensus_10runs_auc_mtx_binarized_threshold_applied.csv")
seaborn.clustermap(binarized[0])

auc_mtx_mait_loc = "/home/tonya/Downloads/thymus.mait_03_02_23_raw_counts_multirun_03192023/thymus.mait_03_02_23_raw_counts_matrix_consensus_10runs_auc_mtx.csv"
auc_mtx = pandas.read_csv(auc_mtx_mait_loc, index_col= 0)
binarized = binarize(auc_mtx, num_workers = 2)
binarized[0].to_csv("/home/tonya/Downloads/thymus.mait_03_02_23_raw_counts_multirun_03192023/thymus.mait_03_02_23_raw_counts_matrix_consensus_10runs_auc_mtx_binarized_threshold_applied.csv")
seaborn.clustermap(binarized[0])

auc_mtx_full_loc = "/home/tonya/Downloads/seurat_filtered_harmony_02_15_23_raw_counts_multirun_03202023/seurat_filtered_harmony_02_15_23_raw_counts_matrix_consensus_10runs_auc_mtx.csv"
auc_mtx = pandas.read_csv(auc_mtx_full_loc, index_col= 0)
binarized = binarize(auc_mtx, num_workers = 2)
binarized[0].to_csv("/home/tonya/Downloads/seurat_filtered_harmony_02_15_23_raw_counts_multirun_03202023/seurat_filtered_harmony_02_15_23_raw_counts_matrix_consensus_10runs_auc_mtx_binarized_threshold_applied.csv")
seaborn.clustermap(binarized[0])