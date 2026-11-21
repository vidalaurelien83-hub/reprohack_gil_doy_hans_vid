#!/usr/bin/env Rscript
# =====================================================================
# venn.R - Figure 2B : un seul diagramme de Venn à 4 ensembles, comme
#          dans l'article (Acute/Chronic x up/down mélangés)
#          (R base uniquement : pas de paquet supplémentaire requis)
# Usage   : Rscript venn.R det_dpi11.tsv det_dpi33.tsv prefixe
# Sorties : <prefixe>.pdf          diagramme à 4 ellipses
#           <prefixe>_counts.tsv   effectifs des 15 zones + valeurs de l'article
#           <prefixe>_shared.tsv   transcrits présents dans plusieurs ensembles
#           stdout                 comparaison avec l'article
# =====================================================================
args   <- commandArgs(trailingOnly = TRUE)
f11    <- args[1]
f33    <- args[2]
prefix <- args[3]

d11 <- read.delim(f11, stringsAsFactors = FALSE)
d33 <- read.delim(f33, stringsAsFactors = FALSE)

# --- Les 4 ensembles, dans l'ordre et avec les noms de l'article ----
sets <- list(
  Acute_upregulated     = d11$transcript_id[d11$log2FoldChange > 0],
  Acute_downregulated   = d11$transcript_id[d11$log2FoldChange < 0],
  Chronic_upregulated   = d33$transcript_id[d33$log2FoldChange > 0],
  Chronic_downregulated = d33$transcript_id[d33$log2FoldChange < 0])
short <- c('AU', 'AD', 'CU', 'CD')       # codes courts des zones

# Appartenance de chaque transcrit aux 4 ensembles -> code de zone ("AU&CU"...)
ids    <- unique(unlist(sets))
member <- sapply(sets, function(s) ids %in% s)
if (!is.matrix(member)) member <- matrix(member, ncol = 4)
code   <- apply(member, 1, function(m) paste(short[m], collapse = '&'))

# Les 15 zones possibles (toutes les combinaisons non vides)
all_zones <- unlist(lapply(1:4, function(k)
  apply(combn(short, k), 2, paste, collapse = '&')))
ours <- setNames(as.integer(table(factor(code, levels = all_zones))), all_zones)

# Valeurs publiées (Hu et al. 2020, Fig. 2B) ; zones absentes = 0
paper <- setNames(integer(length(all_zones)), all_zones)
paper[c('AU', 'AD', 'CU', 'CD', 'AU&CU', 'AD&CU', 'AU&CD', 'AD&CD')] <-
      c(169,  139, 1288,  165,     618,       2,       2,       6)

# --- Géométrie : 4 ellipses disposées comme VennDiagram::draw.quad.venn
ell <- data.frame(cx  = c(0.35, 0.50, 0.50, 0.65),
                  cy  = c(0.47, 0.57, 0.57, 0.47),
                  a   = c(0.35, 0.35, 0.35, 0.35),
                  b   = c(0.20, 0.15, 0.15, 0.20),
                  rot = c(-45, -45, 45, 45) * pi / 180)
cols <- c('#4FA8F0', '#F2F24A', '#F25C5C', '#5CE65C')   # bleu, jaune, rouge, vert

ell_poly <- function(i, n = 300) {
  t <- seq(0, 2 * pi, length.out = n)
  with(ell[i, ], cbind(x = cx + a * cos(t) * cos(rot) - b * sin(t) * sin(rot),
                       y = cy + a * cos(t) * sin(rot) + b * sin(t) * cos(rot)))
}
inside <- function(i, x, y) with(ell[i, ], {
  dx <- x - cx; dy <- y - cy
  u  <-  dx * cos(rot) + dy * sin(rot)
  v  <- -dx * sin(rot) + dy * cos(rot)
  (u / a)^2 + (v / b)^2 <= 1
})

# Position de chaque étiquette : point de la zone le plus éloigné des bords
# (calculée sur une grille, donc valable quelle que soit la forme de la zone)
g     <- expand.grid(x = seq(0, 1, length.out = 220), y = seq(0, 1, length.out = 220))
in_g  <- sapply(1:4, function(i) inside(i, g$x, g$y))
gcode <- apply(in_g, 1, function(m) paste(short[m], collapse = '&'))
dmin  <- rep(Inf, nrow(g))
for (i in 1:4) {
  p <- ell_poly(i, 200)
  for (j in seq_len(nrow(p)))
    dmin <- pmin(dmin, sqrt((g$x - p[j, 1])^2 + (g$y - p[j, 2])^2))
}
label_xy <- t(sapply(all_zones, function(z) {
  k <- which(gcode == z)
  if (!length(k)) return(c(NA, NA))
  k <- k[which.max(dmin[k])]
  c(g$x[k], g$y[k])
}))

# --- Dessin ----------------------------------------------------------
pdf(paste0(prefix, '.pdf'), width = 7, height = 5.9)
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0.12, 1), asp = 1)
for (i in 1:4) polygon(ell_poly(i), col = adjustcolor(cols[i], alpha.f = 0.5), border = NA)
for (i in 1:4) polygon(ell_poly(i), col = NA, border = 'grey30', lwd = 0.8)
text(label_xy[, 1], label_xy[, 2], labels = ours[all_zones], cex = 1.25)
# noms des ensembles, placés comme dans l'article
text(c(0.06, 0.30, 0.70, 0.94), c(0.83, 0.97, 0.97, 0.83),
     labels = names(sets), cex = 0.85, xpd = NA)
invisible(dev.off())

# --- Tables et comparaison avec l'article ---------------------------
tab <- data.frame(zone = all_zones, ours = ours, paper = paper,
                  diff = ours - paper, row.names = NULL)
write.table(tab, paste0(prefix, '_counts.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)

multi <- rowSums(member) > 1
if (any(multi)) {
  shared <- data.frame(zone = code[multi], transcript_id = ids[multi])
  ann <- rbind(d11, d33)
  ann <- ann[!duplicated(ann$transcript_id), ]
  shared <- cbind(shared, ann[match(shared$transcript_id, ann$transcript_id), c('gene_id', 'gene_name')],
                  log2FC_11 = d11$log2FoldChange[match(shared$transcript_id, d11$transcript_id)],
                  log2FC_33 = d33$log2FoldChange[match(shared$transcript_id, d33$transcript_id)])
  shared <- shared[order(shared$zone, shared$transcript_id), ]
} else {
  shared <- data.frame(zone = character(), transcript_id = character())
}
write.table(shared, paste0(prefix, '_shared.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)

cat(prefix, ': zones non vides (nous / article)  AU=Acute_up AD=Acute_down CU=Chronic_up CD=Chronic_down\n')
for (z in all_zones[ours > 0 | paper > 0])
  cat(sprintf('  %-12s %5d / %5d\n', z, ours[z], paper[z]))
