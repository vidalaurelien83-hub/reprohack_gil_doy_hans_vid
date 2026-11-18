# ReproHackathon 2026 – *Staphylococcus aureus* persisters


Le workflow est écrit avec **Nextflow** et chaque étape tourne dans un conteneur **Apptainer**, avec les versions d'outils de l'article.

## Différentes étapes de l'exécution

### Récupération du projet et installation des outils sur la machine virtuelle 

git clone https://github.com/vidalaurelien83-hub/reprohack_gil_doy_hans_vid.git

cd reprohack_gil_doy_hans_vid/Staphylococcus_aureus

bash setup.sh

exec bash

cd ~/reprohack_gil_doy_hans_vid/Staphylococcus_aureus

### Construire les images Apptainer

bash build_images.sh

### Lancer le workflow 

bash run.sh

`run.sh` fixe la version de Nextflow (`NXF_VER=25.04.8`) : elle est téléchargée automatiquement au premier lancement.

Remarque : les annotations KEGG (BRITE `sao00001`) et AureoWiki utilisées pour la Fig. 3c sont téléchargées au moment de l'exécution (étape `DOWNLOAD_ANNOT`) ; ces bases évoluent, le nombre de gènes de la traduction peut donc légèrement varier selon la date d'exécution.

### Récupération des figures

Dans le terminal du PC : 

scp "ubuntu@<IP-DE-LA-VM>:~/reprohack_gil_doy_hans_vid/Staphylococcus_aureus/results/deseq2/*.pdf" .


## Résultats

Figures obtenues avec ce workflow (copie de `results/deseq2/`) :
- [figures/MA_plot_all.pdf](figures/MA_plot_all.pdf) 
- [figures/MA_plot_translation.pdf](figures/MA_plot_translation.pdf)

