#!/bin/bash

#SBATCH --nodes=1 # use one node
#SBATCH --time=01:00:00 #10 minutes
#SBATCH --account=amc-general # normal, amc, long, mem (use mem when using the amem partition)
#SBATCH --partition=amilan # amilian, ami100, aa100, amem, amc
#SBATCH --ntasks=3 # total processes/threads
#SBATCH --job-name=aucell
#SBATCH --output=scenic_merged_10runs_full_object_aucell_03212023_%J.log
#SBATCH --mem=10G # suffix K,M,G,T can be used with default to M
#SBATCH --mail-user=tonya.brunetti@cuanschutz.edu
#SBATCH --mail-type=END
#SBATCH --error=scenic_merged_10runs_full_object_aucell_03212023_%J.err

module load singularity/3.7.4

singularityContainer="/projects/brunetti@xsede.org/singularity_containers/aertslab-pyscenic-0.12.1.sif"
inAndOutDir="/scratch/alpine/brunetti@xsede.org/scenic_gep_data_with_threshold_50target_70regulon/"
exprMatrix="full_obj_seurat_filtered_harmony_02_15_23_raw_counts_matrix_transposed.csv" # IMPORTANT! MAKE SURE TO GO BACK TO FULL MATRIX, NOT SUBSET SINCE NEED BACKGROUND GENES TO PROPERLY ESTIMATE BACKGROUND!!!
regulonsConsensus="GEP12_seurat_filtered_harmony_02_15_23_raw_counts_multirun_04192023_reg70-target50_regulons_pruned_final.gmt"
aucellOutFile="GEP12_seurat_filtered_harmony_02_15_23_raw_counts_multirun_04192023_reg70-target50_consensus_10runs_auc_mtx.csv"

echo "Running Step3: AUCELL on 10 merged reproducible grn and regulon runs"

/usr/bin/time singularity run -B ${inAndOutDir}:/data ${singularityContainer} pyscenic aucell /data/${exprMatrix} /data/${regulonsConsensus} -o /data/${aucellOutFile}
