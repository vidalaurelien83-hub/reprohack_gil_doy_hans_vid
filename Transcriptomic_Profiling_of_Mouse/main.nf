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
// Les FASTQ de l'ENA sont déjà nettoyés : leur nombre de reads est exactement
// la colonne "Clean reads" du Tableau 1. --skip_trimming les donne
// directement à Salmon, sans repasser par Trimmomatic.
params.skip_trimming = false
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
    // Les FASTQ nettoyés (~5 Go par échantillon) restent dans work/ : on ne
    // publie que le résumé (reads conservés), à comparer aux >60 M de l'article
    publishDir "${params.outdir}/trimming", mode: 'copy', pattern: '*.summary.txt'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}_R1.paired.fastq.gz"), path("${sample}_R2.paired.fastq.gz"), emit: paired
    path "${sample}.trimmomatic.summary.txt", emit: summary

    script:
    """
    set -euo pipefail
    trimmomatic PE -threads ${task.cpus} -phred33 \\
      -summary ${sample}.trimmomatic.summary.txt \\
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

// Analyse différentielle DESeq2 1.26 (script R : bin/deseq2.R)
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
    path 'det_dpi*.tsv', emit: dets          // det_dpi{11,33}.tsv + variantes _padj
    path 'deseq2_summary.txt', emit: summary

    script:
    """
    set -euo pipefail
    Rscript "${projectDir}/bin/deseq2.R" '${samplesheet}' '${gtf}' *.quant.sf | tee deseq2_summary.txt
    """
}

// Figure 2B : Venn des DETs 11 vs 33 dpi (script R : bin/venn.R)
// Deux versions : seuil principal (p < 0,05) et variante padj, pour comparer
process VENN {
    container = "${projectDir}/images/deseq2_1.26.0.sif"
    publishDir "${params.outdir}/figures", mode: 'copy'

    input:
    path det_files

    output:
    path 'venn_2B*.pdf', emit: pdf
    path 'venn_2B*.tsv', emit: tables
    path 'venn_summary.txt', emit: summary

    script:
    """
    set -euo pipefail
    Rscript "${projectDir}/bin/venn.R" det_dpi11.tsv det_dpi33.tsv venn_2B > venn_summary.txt
    Rscript "${projectDir}/bin/venn.R" det_dpi11_padj.tsv det_dpi33_padj.tsv venn_2B_padj >> venn_summary.txt
    cat venn_summary.txt
    """
}

// Enrichissement GO/KEGG (script R : bin/enrichment.R)
process CLUSTERPROFILER {
    cpus 2
    container = "${projectDir}/images/cluster_profiler.sif"
    publishDir "${params.outdir}/enrichment", mode: 'copy'

    input:
    path det_files
    path result_files

    output:
    path 'GO_dpi*.tsv', emit: go
    path 'KEGG_dpi*.tsv', emit: kegg

    script:
    """
    Rscript "${projectDir}/bin/enrichment.R"
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
    if (params.skip_trimming) {
        reads_for_quant = reads.reads
    } else {
        reads_for_quant = TRIMMOMATIC(reads.reads).paired
    }
    index = SALMON_INDEX(refs.refs.map { it[0] })
    quantified = SALMON_QUANT(reads_for_quant, index.index)

    all_quants = quantified.quant.map { sample, quant -> quant }.collect()
    sheet_for_deseq = Channel.fromPath(params.samplesheet)
    gtf_for_deseq = refs.refs.map { it[1] }
    differential = DESEQ2(all_quants, sheet_for_deseq, gtf_for_deseq)
    VENN(differential.dets)

    if (params.run_enrichment) {
        CLUSTERPROFILER(differential.dets, differential.results)
    }
}



