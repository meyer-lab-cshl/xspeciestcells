#!/bin/bash
  
#SBATCH --account=amc-general # normal, amc, long, mem (use mem when using the amem partition)
#SBATCH --partition=amilan # amilian, ami100, aa100, amem, amc

echo "Proccessing seed id: $seed for full_object_scenic_GEP3:"
echo "     seed: ${seed}"
echo "     processName: ${processName}"
echo "     outDir: ${outDir}"
echo "     inDir: ${inDir}"
echo "     scenicDb1: ${scenicDb1}"
echo "     scenicDb2: ${scenicDb2}"
echo "     tfs: ${tfs}"
echo "     sif: ${sif}"
echo "     indat: ${indat}"
echo "     prefix: ${prefix}"
echo "     cpus: ${cpus}"

echo "Making the following directory in the outDir path: seed_${seed}_scenic_${prefix}"

mkdir ${outDir}"seed_"${seed}"_scenic_"${prefix}

module load singularity/3.7.4

echo "Running Step1: GRN"

/usr/bin/time singularity run -B ${outDir}"seed_"${seed}"_scenic_"${prefix}:/data,${inDir}:/input,${scenicDb1}:/database1,${scenicDb2}:/database2,${tfs}:/tfsdir ${sif} pyscenic grn --num_workers ${cpus} -o /data/${seed}"_"${prefix}"_adjacencies.tsv" /input/${indat} /tfsdir/allTFs_hg38.txt

wait 

echo "Running Step2: CTX"

/usr/bin/time singularity run -B ${outDir}"seed_"${seed}"_scenic_"${prefix}:/data,${inDir}:/input,${scenicDb1}:/database1,${scenicDb2}:/database2,${tfs}:/tfsdir ${sif} pyscenic ctx /data/${seed}"_"${prefix}"_adjacencies.tsv" /database1/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /database1/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname /database2/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /input/${indat} --mode "custom_multiprocessing" --output /data/${seed}"_"${prefix}"_regulons.csv" --num_workers ${cpus}

wait

echo "Finished running!!"
