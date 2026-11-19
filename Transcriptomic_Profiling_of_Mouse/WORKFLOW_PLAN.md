# Plan de workflow — Transcriptomic Profiling of Mouse

## Objectif

Reproduire le traitement RNA-seq de Hu et al. (2020) en gardant visibles les paramètres non précisés ou ambigus. Le papier décrit Trimmomatic 0.39, Salmon 1.2.0, la référence souris GRCm38.p6 / Ensembl 100 et DESeq2 avec R 3.6.1. Il rapporte deux temps (11 et 33 dpi), avec trois échantillons par condition à chaque temps.

## Structure

Garder les définitions Apptainer présentes dans `container/`. Reprendre du dossier voisin uniquement l’organisation Nextflow : `setup.sh`, `run.sh`, `nextflow.config`, `main.nf`, `samples.csv`, `bin/`, `docs/` et `results/` généré à l’exécution. Les commandes Bowtie, Cutadapt, featureCounts, les échantillons et les résultats de `Staphylococcus_aureus/` ne s’appliquent pas ici.

## Étapes du workflow

1. Charger et valider le manifeste des 12 runs SRA.
2. Télécharger les FASTQ appariés et vérifier leurs MD5.
3. Télécharger le cDNA et le GTF de souris depuis Ensembl release 100.
4. Nettoyer les lectures avec Trimmomatic 0.39.
5. Indexer le transcriptome puis quantifier avec Salmon 1.2.0.
6. Importer les quantifications avec tximport et lancer séparément les contrastes infecté/contrôle à 11 et 33 dpi avec DESeq2 1.26.0.
7. Conserver les p-values et p-values ajustées, puis comparer les résultats aux tables supplémentaires publiées.

## Conteneurs déjà disponibles

- `wget_1.21.3.def` : téléchargements
- `trimmomatic_0.39.def` : nettoyage
- `salmon_1.2.0.def` : indexation et quantification
- `deseq2_1.26.0.def` : R 3.6.1, DESeq2 1.26.0 et tximport 1.14.0

`build_images.sh` construit les images `.sif`; il ne lance pas l’analyse.

## Points reproductibilité à documenter

- Les 12 runs du BioProject PRJNA483261 sont consignés dans `samples.csv`, avec URLs, MD5 et tailles ENA. SRR7657873 est infecté à 33 dpi (T332) et SRR7657882 est contrôle à 33 dpi (C332).
- L’article ne fournit pas les paramètres complets de Trimmomatic. Les paramètres retenus doivent être des options configurables, consignées comme choix de reconstruction.
- Le type de librairie Salmon sera détecté (`-l A`) si aucun type explicite n’est documenté; enregistrer la version et la commande.
- Les Méthodes indiquent p < 0,05 et |log2FC| >= 1,5; les Résultats indiquent p ajusté < 0,05 et |log2FC| > 1,5. Garder les deux mesures et expliciter le filtre utilisé pour reproduire les nombres publiés.
- Les analyses GO/KEGG, réseaux de régulation et interactions hôte-parasite forment une phase aval distincte.

## Répartition du travail

Les personnes peuvent prendre les lots ci-dessous. Les fichiers et paramètres partagés sont à valider ensemble avant le lancement complet.

### Lot A — Données et références

- Contrôler les 12 lignes de `samples.csv` contre ENA, les paires R1/R2 et les MD5.
- Fixer les URL Ensembl 100 du cDNA et du GTF, puis enregistrer leurs checksums et l’annotation utilisée.
- Livrables : manifeste validé et fiche `docs/references.md`.
- Peut avancer en parallèle du lot B.

### Lot B — Orchestration Nextflow et quantification

- Écrire `setup.sh`, `run.sh`, `nextflow.config` et `main.nf` pour téléchargement, vérification MD5, Trimmomatic, index Salmon et quantification.
- Relier les processus aux `.sif` existants; garder les résultats intermédiaires dans `work/` et publier les sorties utiles sous `results/`.
- Documenter les paramètres Trimmomatic choisis comme hypothèses, sans les présenter comme paramètres rapportés par le papier.
- Dépend de l’interface du manifeste définie au lot A.

### Lot C — Analyse différentielle

- Écrire le script R `tximport` + DESeq2 à partir des `quant.sf` et du manifeste.
- Produire les contrastes à 11 et 33 dpi séparément, les tableaux complets avec p-value/padj/log2FC et un résumé des seuils.
- Définir avec le lot B les noms et chemins des quantifications avant intégration.
- Peut préparer le script sur une petite entrée synthétique pendant que le lot B construit le pipeline.

### Lot D — Validation et comparaison

- Récupérer les tables supplémentaires de l’article et confronter les sorties (nombre de gènes, sens des effets, seuils, identifiants).
- Consigner divergences, paramètres testés et limites d’interprétation dans `docs/validation.md`.
- Démarre quand des sorties du lot C sont disponibles; peut préparer en amont une grille de comparaison.

### Intégration commune

- Les lots A à C doivent converger sur le schéma du manifeste et les noms de fichiers avant la première exécution complète.
- La première exécution doit être un essai réduit sur un échantillon, puis l’exécution des 12 échantillons; garder logs et versions des outils.
- Les analyses aval et les pushs GitHub seront décidés après inspection des premiers résultats, conformément au plan de travail.
