#!/usr/bin/env Rscript
# =====================================================================
# enrichment.R - Enrichissement GO (BP) et KEGG des gènes DE, par temps
#                (R 3.6.2 / clusterProfiler, Bioconductor 3.10)
# Usage   : Rscript enrichment.R   (lit deg_dpi*.tsv et results_dpi*.tsv
#           dans le dossier courant)
# Sorties : GO_dpi{11,33}.tsv, KEGG_dpi{11,33}.tsv
# =====================================================================
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

deg_files    <- list.files(pattern = '^deg_dpi[0-9]+\\.tsv$')
result_files <- list.files(pattern = '^results_dpi[0-9]+\\.tsv$')
if (!length(deg_files)) stop('No adjusted-p-value DEG tables found')

for (f in deg_files) {
  dpi  <- sub('^deg_dpi([0-9]+)\\.tsv$', '\\1', f)
  full <- paste0('results_dpi', dpi, '.tsv')
  if (!(full %in% result_files)) stop('Missing full result table: ', full)

  deg <- read.delim(f, stringsAsFactors = FALSE)
  all <- read.delim(full, stringsAsFactors = FALSE)
  ids          <- unique(deg$gene_id)
  universe_ids <- unique(all$gene_id[!is.na(all$gene_id)])

  # Ensembl -> Entrez (identifiants attendus par enrichGO / enrichKEGG)
  mapped   <- bitr(ids, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
  universe <- bitr(universe_ids, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)

  go <- enrichGO(gene = unique(mapped$ENTREZID), universe = unique(universe$ENTREZID),
                 OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = 'BP',
                 pAdjustMethod = 'BH', readable = TRUE)
  kegg <- enrichKEGG(gene = unique(mapped$ENTREZID), universe = unique(universe$ENTREZID),
                     organism = 'mmu', pAdjustMethod = 'BH')

  write.table(as.data.frame(go), paste0('GO_dpi', dpi, '.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  write.table(as.data.frame(kegg), paste0('KEGG_dpi', dpi, '.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
}
