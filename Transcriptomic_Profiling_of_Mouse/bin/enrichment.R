#!/usr/bin/env Rscript
# =====================================================================
# enrichment.R - Enrichissement GO (BP) et KEGG des DETs, par temps
#                (R 3.6.2 / clusterProfiler, Bioconductor 3.10)
# Usage   : Rscript enrichment.R   (lit det_dpi*.tsv et results_dpi*.tsv
#           dans le dossier courant, sorties de deseq2.R)
# Sorties : GO_dpi{11,33}.tsv, KEGG_dpi{11,33}.tsv
# Note    : GO/KEGG travaillent sur des gènes : chaque DET est ramené à
#           son gene_id (plusieurs transcrits d'un gène = un seul gène).
# =====================================================================
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

det_files    <- list.files(pattern = '^det_dpi[0-9]+\\.tsv$')
result_files <- list.files(pattern = '^results_dpi[0-9]+\\.tsv$')
if (!length(det_files)) stop('No DET tables (det_dpi*.tsv) found')

for (f in det_files) {
  dpi  <- sub('^det_dpi([0-9]+)\\.tsv$', '\\1', f)
  full <- paste0('results_dpi', dpi, '.tsv')
  if (!(full %in% result_files)) stop('Missing full result table: ', full)

  det <- read.delim(f, stringsAsFactors = FALSE)
  all <- read.delim(full, stringsAsFactors = FALSE)
  ids          <- unique(det$gene_id[!is.na(det$gene_id)])
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
