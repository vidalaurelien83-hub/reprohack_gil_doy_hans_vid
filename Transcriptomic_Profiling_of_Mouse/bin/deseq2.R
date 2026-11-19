#!/usr/bin/env Rscript
# =====================================================================
# deseq2.R - Transcrits différentiellement exprimés (DETs)
#            infectés vs contrôles, séparément à 11 et 33 dpi
#            (R 3.6.1 / DESeq2 1.26.0 / tximport 1.14.0)
# Usage   : Rscript deseq2.R samples.csv annotation.gtf S1.quant.sf ... S12.quant.sf
# Sorties : results_dpi{11,33}.tsv     résultats complets (un transcrit par ligne)
#           det_dpi{11,33}.tsv         DETs : p < 0,05 et |log2FC| >= 1,5
#                                      (seuil des Méthodes / consigne du hackathon)
#           det_dpi{11,33}_padj.tsv    variante : padj < 0,05 et |log2FC| > 1,5
#                                      (formulation de la partie Résultats)
#           stdout                     nombre de DETs comparé à l'article
# =====================================================================
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
})

args        <- commandArgs(trailingOnly = TRUE)
sheet_file  <- args[1]
gtf_file    <- args[2]
quant_files <- args[-c(1, 2)]
if (length(quant_files) != 12) stop('Expected 12 Salmon quant.sf files; got ', length(quant_files))

# Valeurs publiées (Hu et al. 2020, Fig. 2) : total, surexprimés, sous-exprimés
paper <- list('11' = c(total = 936,  up = 789,  down = 147),
              '33' = c(total = 2081, up = 1908, down = 173))

# --- 1. Plan d'expérience (samples.csv) ------------------------------
samples <- read.csv(sheet_file, stringsAsFactors = FALSE, check.names = FALSE)
names(quant_files) <- sub('\\.quant\\.sf$', '', basename(quant_files))
if (!setequal(names(quant_files), samples$sample)) stop('Sample names in quant files do not match samples.csv')
samples <- samples[match(names(quant_files), samples$sample), ]
rownames(samples) <- samples$sample

# --- 2. Annotation transcrit -> gène (GTF), utilisée seulement pour
#        annoter les tables : l'analyse elle-même reste au niveau transcrit
gtf   <- read.delim(gtf_file, header = FALSE, sep = '\t', comment.char = '#',
                    quote = '', stringsAsFactors = FALSE)
gtf   <- gtf[gtf[[3]] == 'transcript', ]
attr_of <- function(key) sub(paste0('.*', key, ' "([^"]+)".*'), '\\1', gtf[[9]])
annot <- data.frame(tx        = attr_of('transcript_id'),
                    gene_id   = attr_of('gene_id'),
                    gene_name = attr_of('gene_name'),
                    biotype   = attr_of('transcript_biotype'),
                    stringsAsFactors = FALSE)
# gene_name absent de certaines lignes : sub() renvoie alors la ligne entière
annot$gene_name[grepl('"', annot$gene_name)] <- NA
annot$biotype[grepl('"', annot$biotype)]     <- NA
annot <- annot[!duplicated(annot$tx), ]
if (nrow(annot) == 0) stop('No transcript lines found in GTF')

# --- 3. Un DESeq2 séparé par temps (infected vs control) -------------
# Choix : deux tests indépendants (comme les deux comparaisons In vs WT
# de l'article) plutôt qu'un modèle d'interaction infection x temps.
for (dpi in c(11, 33)) {
  coldata <- samples[samples$dpi == dpi, , drop = FALSE]
  # "control" = référence : log2FC > 0 veut dire plus exprimé chez les infectés
  coldata$status <- factor(coldata$status, levels = c('control', 'infected'))
  if (nrow(coldata) != 6 || any(table(coldata$status) != 3)) stop('Expected 3 infected and 3 control samples at dpi ', dpi)

  # txOut = TRUE : comptages par transcrit (pas d'agrégation par gène),
  # car l'article compte des transcrits (DETs)
  txi <- tximport(quant_files[coldata$sample], type = 'salmon', txOut = TRUE)
  dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~status)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c('status', 'infected', 'control'), alpha = 0.05)

  out <- as.data.frame(res)
  out$transcript_id <- rownames(out)
  # quant.sf : ENSMUST...1 (versionné) ; GTF Ensembl : ENSMUST... (sans version)
  tx_nover <- sub('\\.[0-9]+$', '', out$transcript_id)
  m <- match(tx_nover, annot$tx)
  out$gene_id   <- annot$gene_id[m]
  out$gene_name <- annot$gene_name[m]
  out$biotype   <- annot$biotype[m]
  out <- out[, c('transcript_id', 'gene_id', 'gene_name', 'biotype', 'baseMean',
                 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj')]
  out <- out[order(out$pvalue, na.last = TRUE), ]
  write.table(out, paste0('results_dpi', dpi, '.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)

  # --- 4. Sélection des DETs (les deux formulations de l'article) -----
  det      <- out[!is.na(out$pvalue) & out$pvalue < 0.05 & abs(out$log2FoldChange) >= 1.5, ]
  det_padj <- out[!is.na(out$padj)   & out$padj   < 0.05 & abs(out$log2FoldChange) >  1.5, ]
  write.table(det,      paste0('det_dpi', dpi, '.tsv'),      sep = '\t', quote = FALSE, row.names = FALSE)
  write.table(det_padj, paste0('det_dpi', dpi, '_padj.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)

  count <- function(d) c(total = nrow(d), up = sum(d$log2FoldChange > 0), down = sum(d$log2FoldChange < 0))
  p <- paper[[as.character(dpi)]]
  cat(sprintf('dpi %d  [total / up / down]\n', dpi))
  cat(sprintf('  article                       : %5d / %5d / %5d\n', p['total'], p['up'], p['down']))
  x <- count(det)
  cat(sprintf('  p < 0.05,    |log2FC| >= 1.5  : %5d / %5d / %5d\n', x['total'], x['up'], x['down']))
  x <- count(det_padj)
  cat(sprintf('  padj < 0.05, |log2FC| >  1.5  : %5d / %5d / %5d\n', x['total'], x['up'], x['down']))
}
