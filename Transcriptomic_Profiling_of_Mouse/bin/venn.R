#!/usr/bin/env Rscript
# =====================================================================
# venn.R - Figure 2B : DETs communs / spécifiques à 11 et 33 dpi,
#          séparément pour les transcrits surexprimés et sous-exprimés
#          (R base uniquement : pas de paquet supplémentaire requis)
# Usage   : Rscript venn.R det_dpi11.tsv det_dpi33.tsv prefixe
# Sorties : <prefixe>.pdf          les deux diagrammes (up, down)
#           <prefixe>_counts.tsv   effectifs de chaque zone + valeurs de l'article
#           <prefixe>_shared.tsv   liste des transcrits communs aux deux temps
#           stdout                 comparaison avec l'article
# =====================================================================
args   <- commandArgs(trailingOnly = TRUE)
f11    <- args[1]
f33    <- args[2]
prefix <- args[3]

d11 <- read.delim(f11, stringsAsFactors = FALSE)
d33 <- read.delim(f33, stringsAsFactors = FALSE)

# Ensembles de transcrits par sens de variation
sets <- list(
  up   = list(d11 = d11$transcript_id[d11$log2FoldChange > 0],
              d33 = d33$transcript_id[d33$log2FoldChange > 0]),
  down = list(d11 = d11$transcript_id[d11$log2FoldChange < 0],
              d33 = d33$transcript_id[d33$log2FoldChange < 0]))

# Valeurs publiées (Hu et al. 2020, Fig. 2A/2B) : 11 seul, commun, 33 seul
# (789 up à 11 dpi, 1908 up à 33 dpi, 618 communs ; 147 / 173 / 6 pour down)
paper <- list(up   = c(only11 = 789 - 618, shared = 618, only33 = 1908 - 618),
              down = c(only11 = 147 - 6,   shared = 6,   only33 = 173 - 6))

zones <- function(s) c(only11 = length(setdiff(s$d11, s$d33)),
                       shared = length(intersect(s$d11, s$d33)),
                       only33 = length(setdiff(s$d33, s$d11)))

# --- Dessin d'un Venn à deux ensembles en R base --------------------
draw_venn <- function(n, main, col) {
  plot.new()
  plot.window(xlim = c(-2.2, 2.2), ylim = c(-1.7, 1.7), asp = 1)
  t <- seq(0, 2 * pi, length.out = 200)
  for (cx in c(-0.6, 0.6)) {
    polygon(cx + 1.2 * cos(t), 1.2 * sin(t),
            col = adjustcolor(col, alpha.f = 0.25), border = col, lwd = 2)
  }
  text(c(-1.15, 0, 1.15), 0, labels = n, cex = 1.6, font = 2)
  text(c(-0.9, 0.9), 1.45, labels = c('11 dpi', '33 dpi'), cex = 1.2)
  title(main = main, cex.main = 1.3)
}

counts <- lapply(sets, zones)

pdf(paste0(prefix, '.pdf'), width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
# titres en ASCII : évite les problèmes d'accents si la locale du conteneur est "C"
draw_venn(counts$up,   'Up-regulated DETs',   '#C0392B')
draw_venn(counts$down, 'Down-regulated DETs', '#2471A3')
invisible(dev.off())

# --- Tables et comparaison avec l'article ---------------------------
tab <- do.call(rbind, lapply(names(counts), function(k)
  data.frame(direction = k, zone = names(counts[[k]]),
             ours = as.integer(counts[[k]]), paper = as.integer(paper[[k]]))))
tab$diff <- tab$ours - tab$paper
write.table(tab, paste0(prefix, '_counts.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)

shared <- do.call(rbind, lapply(names(sets), function(k) {
  ids <- intersect(sets[[k]]$d11, sets[[k]]$d33)
  if (!length(ids)) return(NULL)
  a <- d11[match(ids, d11$transcript_id), c('transcript_id', 'gene_id', 'gene_name')]
  data.frame(direction = k, a,
             log2FC_11 = d11$log2FoldChange[match(ids, d11$transcript_id)],
             log2FC_33 = d33$log2FoldChange[match(ids, d33$transcript_id)])
}))
if (is.null(shared)) shared <- data.frame(direction = character(), transcript_id = character())
write.table(shared, paste0(prefix, '_shared.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)

cat(prefix, ': zones [11 seul / commun / 33 seul]\n')
for (k in names(counts)) {
  cat(sprintf('  %-4s  nous : %5d / %5d / %5d   article : %5d / %5d / %5d\n', k,
              counts[[k]][1], counts[[k]][2], counts[[k]][3],
              paper[[k]][1], paper[[k]][2], paper[[k]][3]))
}
