## Add some information used by spatialLIBD
spe$key <- paste0(spe$sample_id, "_", colnames(spe))
spe$sum_umi <- colSums(counts(spe))
spe$sum_gene <- colSums(counts(spe) > 0)

## Download the Gencode v32 GTF file and cache it
bfc <- BiocFileCache::BiocFileCache()
gtf_cache <- BiocFileCache::bfcrpath(
  bfc,
  paste0(
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
    "release_32/gencode.v32.annotation.gtf.gz"
  )
)

## Import into R (takes ~1 min)
gtf <- rtracklayer::import(gtf_cache)

## Subset to genes only
gtf <- gtf[gtf$type == "gene"]

## Remove the .x part of the gene IDs
gtf$gene_id <- gsub("\\..*", "", gtf$gene_id)

## Set the names to be the gene IDs
names(gtf) <- gtf$gene_id

## Match the genes
match_genes <- match(rowData(spe)$gene_id, gtf$gene_id)
table(is.na(match_genes))


## Drop the few genes for which we don't have information
spe <- spe[!is.na(match_genes), ]
match_genes <- match_genes[!is.na(match_genes)]

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_name", "gene_type")]

## Save the "interest"ing columns from our original spe object
interesting <- rowData(spe)[, grepl("interest", colnames(rowData(spe)))]

## Add the gene info to our SPE object
rowRanges(spe) <- gtf[match_genes]

## Add back the "interest" coolumns
rowData(spe) <- cbind(rowData(spe), interesting)

## Inspect the gene annotation data we added
rowRanges(spe)

spatialLIBD::vis_gene(
  spe = spe,
  sampleid = "GSM6281322",
  geneid = "CD1D",
)