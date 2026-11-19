#!/bin/bash
# =====================================================================
# run.sh - Lance le workflow Nextflow (projet Mouse, Hu et al. 2020)
# Usage : bash run.sh                         (reprend là où il s'était arrêté grâce à -resume)
#         bash run.sh --run_enrichment true   (ajoute GO/KEGG avec clusterProfiler)
#         Les autres options sont transmises à Nextflow (ex. --outdir autre_dossier)
# Prérequis : setup.sh puis build_images.sh
# =====================================================================
set -euo pipefail

cd "$(dirname "$0")"     # se placer à la racine du projet

# Version de Nextflow figée (téléchargée automatiquement au 1er lancement)
export NXF_VER=25.04.8

# Images nécessaires ; cluster_profiler seulement si l'enrichissement est demandé
requises="wget_1.21.3 trimmomatic_0.39 salmon_1.2.0 deseq2_1.26.0"
if [[ " $* " == *" --run_enrichment true "* ]]; then
    requises="$requises cluster_profiler"
fi

manquantes=""
for name in $requises; do
    [ -f "images/$name.sif" ] || manquantes="$manquantes $name"
done
if [ -n "$manquantes" ]; then
    echo "ERREUR : images manquantes :$manquantes" >&2
    echo "Lancez d'abord : bash build_images.sh" >&2
    exit 1
fi

nextflow run main.nf -resume "$@"
