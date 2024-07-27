from pyscenic.cli.utils import load_signatures
from ctxcore.genesig import Regulon
import os
import collections
import pandas
import typing


#TODO: ADD PARSER AND HELP MESSAGING
grnboostAdjacenciesLoc = ''
regulonLoc = '/home/tonya/tmp_SCENIC/mouse_scenic_runs_10102023/full_innate_Tcell_integration_including_gd/'
outputFile = '/home/tonya/tmp_SCENIC/mouse_scenic_runs_10102023/full_innate_Tcell_integration_including_gd/nes2_maskdropouts_mouse_innate_gd_intergration_50_multirun_reg100-target80_regulons_pruned_final.gmt'
outputFileCSV = '/home/tonya/tmp_SCENIC/mouse_scenic_runs_10102023/full_innate_Tcell_integration_including_gd/nes2_maskdropouts_mouse_innate_gd_intergration_50_multirun_reg100-target80_regulons_pruned_final.csv'
reproducibiltyPercentOutFile = '/home/tonya/tmp_SCENIC/mouse_scenic_runs_10102023/full_innate_Tcell_integration_including_gd/nes2_maskdropouts_mouse_innate_gd_intergration_50_multirun_reg100-target80_regulons_pruned_final_reproducibility_metrics.txt'
targetGeneCounts='/home/tonya/tmp_SCENIC/mouse_scenic_runs_10102023/full_innate_Tcell_integration_including_gd/nes2_maskdropouts_mouse_innate_gd_intergration_50_multirun_reg100-target80_regulons_pruned_final_target_gene_counts.csv'
regulonThresh = 1
geneTargetThresh = 0.80


def reproducibility_across_runs() -> None:
    totalRuns = 0
    allRegulons = []
    geneTargets = {}

    for regulonsFiles in os.listdir(regulonLoc):
        if regulonsFiles.endswith('regulons.csv'):
            totalRuns += 1
            sig = load_signatures(os.path.join(regulonLoc, regulonsFiles)) 
            tfs = [x.transcription_factor for x in sig]
            for tfIdx in range(0, len(sig)):
                if tfs[tfIdx] in geneTargets:
                    geneTargets[tfs[tfIdx]].extend(list(sig[tfIdx].genes))
                else:
                    geneTargets[tfs[tfIdx]] = list(sig[tfIdx].genes)

            allRegulons.extend(tfs)

    reproducibility = collections.Counter(allRegulons)

    passingRegulons = []
    for key,value in reproducibility.items():
        if float(value)/float(totalRuns) >= regulonThresh:
            passingRegulons.append(key)
            
    print('{}/{} regulons pass the {} threshold, indicating {}/{} regulons occur in {}% of {} total runs'.format(len(passingRegulons), len(reproducibility), regulonThresh, len(passingRegulons), len(reproducibility), regulonThresh*100, totalRuns))

    prunedRegulons = {key:collections.Counter(value) for key,value in geneTargets.items() if key in passingRegulons}

    reproducibilityCounts = pandas.DataFrame.from_dict(prunedRegulons, orient='columns')
    reproducibilityCounts.to_csv(targetGeneCounts, index = True)

    finalRegulons = {}
    for key,values in prunedRegulons.items():
        passingGeneTargets = {gene:counts for gene,counts in values.items() if float(counts)/float(reproducibility[key]) >= geneTargetThresh}
        if len(passingGeneTargets) > 0:
            finalRegulons[key] = passingGeneTargets

    # Regulon(name='Aggregated regulon', gene2weight={'target1': 1.0, 'target2': 1.0}, transcription_factor='TF1', gene2occurrence={'target1': 1, 'target2': 1 })
    regulonObjs = []
    for key,value in finalRegulons.items():
        regulonObjs.append(Regulon(name='Aggregated regulon', gene2weight={genes: 1.0 for genes in finalRegulons[key].keys()}, transcription_factor=key, gene2occurrence=finalRegulons[key]))

    eof = 0 
    with open(outputFile, 'w') as gmt_file:
        for regulon in regulonObjs:
            eof += 1
            gmt_file.write('\t'.join([regulon.transcription_factor, regulon.name]))
            gmt_file.write('\t' + '\t'.join(list(regulon.gene2occurrence.keys())))
            if eof != len(regulonObjs):
                gmt_file.write('\n')


    eof = 0 
    with open(outputFileCSV, 'w') as gmt_file:
        for regulon in regulonObjs:
            eof += 1
            gmt_file.write(','.join([regulon.transcription_factor, regulon.name]))
            gmt_file.write(',' + ','.join(list(regulon.gene2occurrence.keys())))
            if eof != len(regulonObjs):
                gmt_file.write('\n')

    # write out text file of the number of times a regulon was reproducible across runs
    eof = 0
    with open(reproducibiltyPercentOutFile, 'w') as repFile:
        for regulons in list(reproducibility):
            eof += 1
            repFile.write('{}: {}/{}'.format(regulons, reproducibility[regulons], totalRuns))
            if eof != len(list(reproducibility)):
                repFile.write('\n')
            

def regulon_pruning(regulons, minGene) -> None:
    # regulons should be the .gmt file output of running reproducibility_across_runs() (each row is regulon, followed by "Aggregated Regulon", followed by tab-separated list of genes in regulon)
    # minGene is inclusive to being kept
    regulons_pruned_out = {}
    regulon_counts = 0
    regulons = "/home/tonya/tmp_SCENIC/scenic_50_runs_dewakss/dewakss_reg70-target50_regulons_pruned_final.gmt"
    passing_regulons = {}
    with open(regulons) as regulon_possibilities:
        for line in regulon_possibilities:
            regulon_counts += 1
            line = line.strip("\n").split("\t")
            if len(line[2:]) >= minGene:
                passing_regulons[line[0]] = line[2:]
            else:
                regulons_pruned_out[line[0]] = line[2:]


def target_pruning(intersectionCSV :  str, colsToUse: list, regulonFileGmt: str, colForPdIndex : int, outputFile:str, outputFileCSV:str) -> None:
    '''
    intersectionCSV: contains 1 or more columns of genes you want to keep in targets with a header
    regulonFileGmt: the output from reproducibility_across_runs() -> outputFile (tab-delimited)
    colsToUse: columns in intersectionCSV to use for keeping targets across regulons
    '''
    intersectionCSV = '/home/tonya/tmp_SCENIC/new_scenic_100_runs_09132023/full_object/nes2_maskdropouts/subset_targets_to_geps/limited_nonimputed_genes_per_gep_post_rank_threshold_k12.csv'
    colsToUse = ['GEP1', 'GEP2', 'GEP3', 'GEP4', 'GEP5', 'GEP6','GEP7', 'GEP8', 'GEP9', 'GEP10', 'GEP11']
    regulonFileGmt  = '/home/tonya/tmp_SCENIC/new_scenic_100_runs_09132023/full_object/nes2_maskdropouts/subset_targets_to_geps/reg100_target100/nes2_maskdropouts_full_obj_seurat_filtered_harmony_08_28_23_raw_counts_100_multirun_reg100-target100_regulons_pruned_final.gmt'
    outputFile = "/home/tonya/tmp_SCENIC/new_scenic_100_runs_09132023/full_object/nes2_maskdropouts/subset_targets_to_geps/reg100_target100/nes2_maskdropouts_full_obj_seurat_filtered_harmony_08_28_23_raw_counts_100_multirun_reg100-target100_regulons_pruned_final_GEP_INTERSECTION_TARGET_PRUNED.gmt"
    outputFileCSV  = "/home/tonya/tmp_SCENIC/new_scenic_100_runs_09132023/full_object/nes2_maskdropouts/subset_targets_to_geps/reg100_target100/nes2_maskdropouts_full_obj_seurat_filtered_harmony_08_28_23_raw_counts_100_multirun_reg100-target100_regulons_pruned_final_GEP_INTERSECTION_TARGET_PRUNED.csv"
    colForPdIndex = 0
    
    keepTargetDf = pandas.read_csv(intersectionCSV,dtype = str, index_col=colForPdIndex, usecols = colsToUse)
    keepTargetDf.columns.name = 'columns'
    # pandas columns becaome all one column in new DF (column -> value mapping so each row is 1 column to 1 value mapping; you get the number of rows per column is equal to the number of values in that column)
    transformedDf = keepTargetDf.stack()   
    transformedDf.name = 'target_genes' # name of new column after stacking
    transformedDf = transformedDf.reset_index()
    keepTargets = set(transformedDf['target_genes'].tolist())
    
    # update regulons with intersection of target genes in CSV vs target genes in regulonFileGmt
    regulonTargetGroups = {}
    with open(regulonFileGmt, 'r') as gmtFile:
        for line in gmtFile:
            sharedTargets = set(line.split('\t')[2:]) & set(keepTargets)
            if len(sharedTargets) != 0: # do not keep any regulons that have 0 target genes
                regulonTargetGroups[line.split('\t')[0]] = sharedTargets
    
    # write output to gmt format
    eof = 0 
    with open(outputFile, 'w') as gmt_file:
        for regulon, targets in regulonTargetGroups.items():
            print(regulon)
            eof += 1
            gmt_file.write('{}\t{}'.format(regulon, "Aggregated regulon"))
            gmt_file.write('\t' + '\t'.join(targets))
            if eof != len(regulonTargetGroups):
                gmt_file.write('\n')

    # write output as csv
    eof = 0 
    with open(outputFileCSV, 'w') as csv_file:
        for regulon, targets in regulonTargetGroups.items():
            eof += 1
            csv_file.write('{},{}'.format(regulon, "Aggregated regulon"))
            csv_file.write(',' + ','.join(targets))
            if eof != len(regulonTargetGroups):
                csv_file.write('\n')

'''
Regulons are derived from adjacencies based on three methods.

The first method to create the TF-modules is to select the best targets for each transcription factor:

1. Targets with importance > the 50th percentile.
2. Targets with importance > the 75th percentile.
3. Targets with importance > the 90th percentile.

The second method is to select the top targets for a given TF:

1. Top 50 targets (targets with highest weight)

'''
def explore_regulon_adjacencies(adjacenciesDir:stsr) -> None:
    from arboreto.algo import grnboost2
    from arboreto.utils import load_tf_names
    from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
    from pyscenic.utils import modules_from_adjacencies, load_motifs
    from pyscenic.prune import prune2df, df2regulons
    from pyscenic.aucell import aucell
    adjFiles = "~/1303_dewakss_imputed_10kcells_seurat_filtered_harmony_08_28_23_raw_counts_matrix_adjacencies.tsv"
    tfs = "/home/tonya/tmp_SCENIC/scenic_50_runs_dewakss/scenic_files/allTFs_hg38.txt" 
    db_names = "/home/tonya/tmp_SCENIC/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
    motif_annotation_file = "/home/tonya/tmp_SCENIC/scenic_50_runs_dewakss/scenic_files/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    # tf_names = load_tf_names(tfs) does same thing as below
    tf_names = [line.strip() for line in open(tfs, 'r')]
    dbs = [RankingDatabase(fname = db_names, name = db_names.split('/')[-1])]
    
    # expr_matrix columns as genes and rows as cell names
    expr_matrix = pandas.read_csv("/home/tonya/tmp_SCENIC/scenic_50_runs_dewakss/subsampled_dewakssimputed_cellsxgenes.csv", sep=',', header=0, index_col=0)
    for adjFiles in os.listdir(adjacenciesDir):
        if adjFiles.endswith('matrix_adjacencies.tsv'):
            grn_matrix = pandas.read_csv(adjFiles, sep = '\t', header=0)
            q90 = grn_matrix[grn_matrix.importance > grn_matrix.importance.quantile(q = 0.90)]
            # modules_from_adjacencies calculates Pearson correlations between TF and target on all cells (rho_mask_dropouts = False)
            modules = list(modules_from_adjacencies(grn_matrix, expr_matrix, rho_mask_dropouts=False)) 
            # modules_from_adjacencies calculates Pearson correlations between TF and target with cells that have a non-zero expression value (rho_mask_dropouts = True)
            modules_ignore_dropout = list(modules_from_adjacencies(grn_matrix, expr_matrix, rho_mask_dropouts=True)) 

            [modules[regulon] for regulon in range(0, len(modules)) if 'ZBTB16' in modules[regulon].name]
            [modules_ignore_dropout[regulon] for regulon in range(0, len(modules)) if 'ZBTB16' in modules[regulon].name]
            [{regulon:modules_ignore_dropout[regulon]} for regulon in range(0, len(modules)) if 'CEBPD' in modules[regulon].name]
            
            df = prune2df(dbs, modules, motif_annotation_file)
            regulonTable = df2regulons(df)
    
    
    feather_file = pandas.read_feather(db_names)
