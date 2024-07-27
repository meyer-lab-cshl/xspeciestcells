#!/bin/bash

#SBATCH --nodes=1 # use one node
#SBATCH --time=00:02:00 #2 minutes
#SBATCH --account=amc-general # normal, amc, long, mem (use mem when using the amem partition)
#SBATCH --partition=amilan # amilian, ami100, aa100, amem, amc
#SBATCH --ntasks=1 # total processes/threads
#SBATCH --job-name=submissionscript
#SBATCH --output=scenic_job_submission_09072023_%J.log
#SBATCH --mem=1M # suffix K,M,G,T can be used with default to M
#SBATCH --mail-user=tonya.brunetti@cuanschutz.edu
#SBATCH --mail-type=END
#SBATCH --error=scenic_job_submission_09072023_%J.err

seedArray=($(shuf -i 1000-9999 -n 10))
processName="seurat_full_object_all_cells"
dataOutdir="/scratch/alpine/brunetti@xsede.org/"
dataIndir="/projects/brunetti@xsede.org/GEPs_seurat_filtered_harmony_08_28_23_raw_counts_SCENIC/"
scenicDb1="/projects/brunetti@xsede.org/reference_data/scenic/resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/"
scenicDb2="/projects/brunetti@xsede.org/reference_data/scenic/resources.aertslab.org/cistarget/motif2tf/"
tfs="/projects/brunetti@xsede.org/reference_data/scenic/resources.aertslab.org/cistarget/tf_lists/"
singularityContainer="/projects/brunetti@xsede.org/singularity_containers/aertslab-pyscenic-0.12.1.sif"
dataset="GEP_3_seurat_filtered_harmony_08_28_23_raw_counts_matrix_transposed.csv"
prefix="GEP3_seurat_filtered_harmony_08_28_23_raw_counts_matrix"
cpus=30

echo ${seedArray[@]}


for i in ${seedArray[@]}
do
	sbatch --job-name=${i}"_fullobj" --ntasks=${cpus} --mem=120G --time=7:00:00 --nodes=1 --output=${dataOutdir}"seed_"${i}"_"${prefix}"_%J.log" --error=${dataOutdir}"seed_"${i}"_"${prefix}"_%J.err" --export=seed=${i},processName=${processName},outDir=${dataOutdir},inDir=${dataIndir},scenicDb1=${scenicDb1},scenicDb2=${scenicDb2},tfs=${tfs},sif=${singularityContainer},indat=${dataset},prefix=${prefix},cpus=${cpus} GEP3_scenic_subprocesses_09072023_full_object.sh
	echo "Submitted seed: $i"
done
