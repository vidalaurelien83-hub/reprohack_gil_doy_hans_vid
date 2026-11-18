#!/bin/bash
# =====================================================================
# run.sh - Lance le workflow Nextflow
# Usage : bash run.sh   (reprend là où il s'était arrêté grâce à -resume)
# Prérequis : setup.sh puis build_images.sh
# =====================================================================
set -euo pipefail

cd "$(dirname "$0")"     # se placer à la racine du projet

# Version de Nextflow figée (téléchargée automatiquement au 1er lancement)
export NXF_VER=25.04.8

# Vérifie que toutes les images Apptainer sont construites
manquantes=""
for def in container/*.def; do
    name=$(basename "$def" .def)
    [ -f "images/$name.sif" ] || manquantes="$manquantes $name"
done
if [ -n "$manquantes" ]; then
    echo "ERREUR : images manquantes :$manquantes" >&2
    echo "Lancez d'abord : bash build_images.sh" >&2
    exit 1
fi

nextflow run main.nf -resume "$@"
