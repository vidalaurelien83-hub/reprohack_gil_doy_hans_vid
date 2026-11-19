#!/usr/bin/env nextflow
/*
 * RNA-seq pipeline — Transcriptomic Profiling of Mouse
 * Hu et al. (2020), workflow reconstruction for the 12 ENA runs in samples.csv.
 * Launch: nextflow run main.nf -resume
 */
nextflow.enable.dsl = 2

params.samplesheet = "${projectDir}/samples.csv"
params.outdir      = "${projectDir}/results"
params.datadir     = "${projectDir}/data"
params.cdna_url    = 'https://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz'
params.gtf_url     = 'https://ftp.ensembl.org/pub/release-100/gtf/mus_musculus/Mus_musculus.GRCm38.100.gtf.gz'
// The paper does not report these exact Trimmomatic options; they are explicit reconstruction choices.
params.trim_adapters = '/opt/conda/share/trimmomatic/adapters/TruSeq3-PE.fa'
params.trim_args     = "ILLUMINACLIP:${params.trim_adapters}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36"
// Salmon library type is automatically inferred; these bias/mapping options are recorded here.
params.salmon_extra  = '--validateMappings --seqBias --gcBias'
// Optional GO/KEGG enrichment with the clusterProfiler image.
params.run_enrichment = false

process DOWNLOAD_FASTQ {
    tag "$sample"
    container = "${projectDir}/images/wget_1.21.3.sif"
    publishDir "${params.datadir}/fastq", mode: 'symlink'

    input:
    tuple val(sample), val(url1), val(md51), val(url2), val(md52)

    output:
    tuple val(sample), path("${sample}_1.fastq.gz"), path("${sample}_2.fastq.gz"), emit: reads

    script:
    """
    set -euo pipefail
    wget --retry-connrefused --waitretry=5 --tries=10 -q -O ${sample}_1.fastq.gz '${url1}'
    wget --retry-connrefused --waitretry=5 --tries=10 -q -O ${sample}_2.fastq.gz '${url2}'
    echo '${md51}  ${sample}_1.fastq.gz' | md5sum -c -
    echo '${md52}  ${sample}_2.fastq.gz' | md5sum -c -
    """
}

process DOWNLOAD_REFERENCE {
    container = "${projectDir}/images/wget_1.21.3.sif"
    publishDir "${params.datadir}/reference", mode: 'copy'

    output:
    tuple path('transcripts.fa'), path('annotation.gtf'), emit: refs

    script:
    """
    set -euo pipefail
    wget --retry-connrefused --waitretry=5 --tries=10 -q -O cdna.fa.gz '${params.cdna_url}'
    wget --retry-connrefused --waitretry=5 --tries=10 -q -O annotation.gtf.gz '${params.gtf_url}'
    gzip -dc cdna.fa.gz > transcripts.fa
    gzip -dc annotation.gtf.gz > annotation.gtf
    """
}

process TRIMMOMATIC {
    tag "$sample"
    cpus 4
    container = "${projectDir}/images/trimmomatic_0.39.sif"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}_R1.paired.fastq.gz"), path("${sample}_R2.paired.fastq.gz"), emit: paired

    script:
    """
    set -euo pipefail
    trimmomatic PE -threads ${task.cpus} -phred33 \\
      ${read1} ${read2} \\
      ${sample}_R1.paired.fastq.gz ${sample}_R1.unpaired.fastq.gz \\
      ${sample}_R2.paired.fastq.gz ${sample}_R2.unpaired.fastq.gz \\
      ${params.trim_args}
    """
}

process SALMON_INDEX {
    cpus 4
    container = "${projectDir}/images/salmon_1.2.0.sif"
    publishDir "${params.datadir}/reference", mode: 'copy'

    input:
    path cdna

    output:
    path 'salmon_index', emit: index

    script:
    """
    salmon index -t ${cdna} -i salmon_index -k 31
    """
}

process SALMON_QUANT {
    tag "$sample"
    cpus 4
    container = "${projectDir}/images/salmon_1.2.0.sif"
    publishDir "${params.outdir}/quantification", mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)
    path index

    output:
    tuple val(sample), path("${sample}.quant.sf"), emit: quant

    script:
    """
    set -euo pipefail
    salmon quant -i ${index} -l A -1 ${read1} -2 ${read2} \\
      -p ${task.cpus} ${params.salmon_extra} -o ${sample}_quant
    cp ${sample}_quant/quant.sf ${sample}.quant.sf
    """
}

process DESEQ2 {
    cpus 4
    container = "${projectDir}/images/deseq2_1.26.0.sif"
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    path quant_files
    path samplesheet
    path gtf

    output:
    path 'results_dpi*.tsv', emit: results
    path 'deg_dpi*.tsv', emit: degs
    path 'deseq2_summary.txt', emit: summary

    script:
    """
    cat > deseq2.R <<'RSCRIPT'
    suppressPackageStartupMessages({
      library(tximport)
      library(DESeq2)
    })
    args <- commandArgs(trailingOnly=TRUE)
    sheet_file <- args[1]
    gtf_file <- args[2]
    quant_files <- args[-c(1,2)]
    if (length(quant_files) != 12) stop('Expected 12 Salmon quant.sf files; got ', length(quant_files))

    samples <- read.csv(sheet_file, stringsAsFactors=FALSE, check.names=FALSE)
    names(quant_files) <- sub('\\\\.quant\\\\.sf$', '', basename(quant_files))
    if (!setequal(names(quant_files), samples$sample)) stop('Sample names in quant files do not match samples.csv')
    samples <- samples[match(names(quant_files), samples$sample), ]
    if (anyNA(samples$sample)) stop('Some quantifications have no matching sample row')

    gtf <- read.delim(gtf_file, header=FALSE, sep='\\t', comment.char='#', quote='', stringsAsFactors=FALSE)
    attrs <- gtf[[9]]
    tx <- sub('.*transcript_id "([^"]+)".*', '\\\\1', attrs)
    gene <- sub('.*gene_id "([^"]+)".*', '\\\\1', attrs)
    tx2gene <- unique(data.frame(tx=tx[tx != attrs], gene=gene[tx != attrs], stringsAsFactors=FALSE))
    tx2gene <- tx2gene[!grepl('transcript_id', tx2gene$tx), ]
    if (nrow(tx2gene) == 0) stop('Could not extract transcript_id/gene_id mapping from GTF')

    rownames(samples) <- samples$sample
    for (dpi in c(11, 33)) {
      keep <- samples$dpi == dpi
      coldata <- samples[keep, , drop=FALSE]
      coldata$status <- factor(coldata$status, levels=c('control', 'infected'))
      if (nrow(coldata) != 6 || any(table(coldata$status) != 3)) stop('Expected 3 infected and 3 control samples at dpi ', dpi)
      files_dpi <- quant_files[coldata$sample]
      if (anyNA(files_dpi)) stop('Missing Salmon quantification at dpi ', dpi)
      txi <- tximport(files_dpi, type='salmon', tx2gene=tx2gene, ignoreTxVersion=TRUE)
      dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=~status)
      dds <- DESeq(dds, quiet=TRUE)
      res <- results(dds, contrast=c('status','infected','control'), alpha=0.05)
      out <- as.data.frame(res)
      out$gene_id <- rownames(out)
      out <- out[, c('gene_id','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')]
      out <- out[order(out$padj, na.last=TRUE), ]
      write.table(out, paste0('results_dpi', dpi, '.tsv'), sep='\\t', quote=FALSE, row.names=FALSE)
      # Save both paper-reported filters because Methods and Results specify different p-value rules.
      by_p <- out[!is.na(out$pvalue) & out$pvalue < 0.05 & abs(out$log2FoldChange) >= 1.5, ]
      by_padj <- out[!is.na(out$padj) & out$padj < 0.05 & abs(out$log2FoldChange) > 1.5, ]
      write.table(by_padj, paste0('deg_dpi', dpi, '.tsv'), sep='\\t', quote=FALSE, row.names=FALSE)
      write.table(by_p, paste0('deg_dpi', dpi, '_methods_p.tsv'), sep='\\t', quote=FALSE, row.names=FALSE)
      cat('dpi', dpi, ':', nrow(by_padj), 'DEGs (padj<0.05, |log2FC|>1.5);', nrow(by_p), 'using Methods p-value threshold\\n')
    }
    RSCRIPT
    set -euo pipefail
    Rscript deseq2.R '${samplesheet}' '${gtf}' *.quant.sf | tee deseq2_summary.txt
    """
}

process CLUSTERPROFILER {
    cpus 2
    container = "${projectDir}/images/cluster_profiler.sif"
    publishDir "${params.outdir}/enrichment", mode: 'copy'

    input:
    path deg_files
    path result_files

    output:
    path 'GO_dpi*.tsv', emit: go
    path 'KEGG_dpi*.tsv', emit: kegg

    script:
    """
    cat > enrich.R <<'RSCRIPT'
    suppressPackageStartupMessages({library(clusterProfiler); library(org.Mm.eg.db)})
    deg_files <- list.files(pattern='^deg_dpi[0-9]+\\\\.tsv$')
    result_files <- list.files(pattern='^results_dpi[0-9]+\\\\.tsv$')
    if (!length(deg_files)) stop('No adjusted-p-value DEG tables found')
    for (f in deg_files) {
      dpi <- sub('^deg_dpi([0-9]+)\\\\.tsv$', '\\\\1', f)
      full <- paste0('results_dpi', dpi, '.tsv')
      if (!(full %in% result_files)) stop('Missing full result table: ', full)
      deg <- read.delim(f, stringsAsFactors=FALSE)
      all <- read.delim(full, stringsAsFactors=FALSE)
      ids <- unique(deg$gene_id)
      universe_ids <- unique(all$gene_id[!is.na(all$gene_id)])
      mapped <- bitr(ids, fromType='ENSEMBL', toType='ENTREZID', OrgDb=org.Mm.eg.db)
      universe <- bitr(universe_ids, fromType='ENSEMBL', toType='ENTREZID', OrgDb=org.Mm.eg.db)
      go <- enrichGO(gene=unique(mapped$ENTREZID), universe=unique(universe$ENTREZID), OrgDb=org.Mm.eg.db, keyType='ENTREZID', ont='BP', pAdjustMethod='BH', readable=TRUE)
      kegg <- enrichKEGG(gene=unique(mapped$ENTREZID), universe=unique(universe$ENTREZID), organism='mmu', pAdjustMethod='BH')
      write.table(as.data.frame(go), paste0('GO_dpi', dpi, '.tsv'), sep='\\t', quote=FALSE, row.names=FALSE)
      write.table(as.data.frame(kegg), paste0('KEGG_dpi', dpi, '.tsv'), sep='\\t', quote=FALSE, row.names=FALSE)
    }
    RSCRIPT
    Rscript enrich.R deg_dpi*.tsv results_dpi*.tsv
    """
}

workflow {
    sample_rows = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            if (!(row.layout?.toUpperCase() == 'PAIRED')) error "Expected paired-end sample: ${row.sample}"
            tuple(row.sample, row.fastq_1_url, row.fastq_1_md5, row.fastq_2_url, row.fastq_2_md5)
        }

    reads = DOWNLOAD_FASTQ(sample_rows)
    refs = DOWNLOAD_REFERENCE()
    trimmed = TRIMMOMATIC(reads.reads)
    index = SALMON_INDEX(refs.refs.map { it[0] })
    quantified = SALMON_QUANT(trimmed.paired, index.index)

    all_quants = quantified.quant.map { sample, quant -> quant }.collect()
    sheet_for_deseq = Channel.fromPath(params.samplesheet)
    gtf_for_deseq = refs.refs.map { it[1] }
    differential = DESEQ2(all_quants, sheet_for_deseq, gtf_for_deseq)

    if (params.run_enrichment) {
        CLUSTERPROFILER(differential.degs, differential.results)
    }
}



