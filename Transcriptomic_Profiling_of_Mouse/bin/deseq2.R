#!/usr/bin/env Rscript
# =====================================================================
# deseq2.R - Analyse différentielle infectés vs contrôles, à 11 et 33 dpi
#            (R 3.6.1 / DESeq2 1.26.0 / tximport 1.14.0)
# Usage   : Rscript deseq2.R samples.csv annotation.gtf S1.quant.sf ... S12.quant.sf
# Sorties : results_dpi{11,33}.tsv           résultats complets
#           deg_dpi{11,33}.tsv               padj < 0,05 et |log2FC| > 1,5 (Résultats)
#           deg_dpi{11,33}_methods_p.tsv     p < 0,05 et |log2FC| >= 1,5 (Méthodes)
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

# --- 1. Plan d'expérience (samples.csv) ------------------------------
samples <- read.csv(sheet_file, stringsAsFactors = FALSE, check.names = FALSE)
names(quant_files) <- sub('\\.quant\\.sf$', '', basename(quant_files))
if (!setequal(names(quant_files), samples$sample)) stop('Sample names in quant files do not match samples.csv')
samples <- samples[match(names(quant_files), samples$sample), ]
if (anyNA(samples$sample)) stop('Some quantifications have no matching sample row')

# --- 2. Table transcrit -> gène extraite du GTF ----------------------
gtf   <- read.delim(gtf_file, header = FALSE, sep = '\t', comment.char = '#',
                    quote = '', stringsAsFactors = FALSE)
attrs <- gtf[[9]]
tx    <- sub('.*transcript_id "([^"]+)".*', '\\1', attrs)
gene  <- sub('.*gene_id "([^"]+)".*', '\\1', attrs)
tx2gene <- unique(data.frame(tx = tx[tx != attrs], gene = gene[tx != attrs],
                             stringsAsFactors = FALSE))
tx2gene <- tx2gene[!grepl('transcript_id', tx2gene$tx), ]
if (nrow(tx2gene) == 0) stop('Could not extract transcript_id/gene_id mapping from GTF')

# --- 3. Un DESeq2 séparé par temps (infected vs control) -------------
rownames(samples) <- samples$sample
for (dpi in c(11, 33)) {
  keep    <- samples$dpi == dpi
  coldata <- samples[keep, , drop = FALSE]
  # "control" = référence : log2FC > 0 veut dire plus exprimé chez les infectés
  coldata$status <- factor(coldata$status, levels = c('control', 'infected'))
  if (nrow(coldata) != 6 || any(table(coldata$status) != 3)) stop('Expected 3 infected and 3 control samples at dpi ', dpi)
  files_dpi <- quant_files[coldata$sample]
  if (anyNA(files_dpi)) stop('Missing Salmon quantification at dpi ', dpi)

  txi <- tximport(files_dpi, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = TRUE)
  dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~status)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c('status', 'infected', 'control'), alpha = 0.05)

  out <- as.data.frame(res)
  out$gene_id <- rownames(out)
  out <- out[, c('gene_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj')]
  out <- out[order(out$padj, na.last = TRUE), ]
  write.table(out, paste0('results_dpi', dpi, '.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)

  # --- 4. Les deux filtres de l'article (Méthodes et Résultats diffèrent)
  by_p    <- out[!is.na(out$pvalue) & out$pvalue < 0.05 & abs(out$log2FoldChange) >= 1.5, ]
  by_padj <- out[!is.na(out$padj) & out$padj < 0.05 & abs(out$log2FoldChange) > 1.5, ]
  write.table(by_padj, paste0('deg_dpi', dpi, '.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  write.table(by_p, paste0('deg_dpi', dpi, '_methods_p.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  cat('dpi', dpi, ':', nrow(by_padj), 'DEGs (padj<0.05, |log2FC|>1.5);',
      nrow(by_p), 'using Methods p-value threshold\n')
}
