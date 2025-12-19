# Analyse transcriptomique TCGA-UCEC

Projet M2 AIDA - Analyse différentielle de l'expression génique dans le cancer de l'endomètre (Uterine Corpus Endometrial Carcinoma)

**Repository** : [https://github.com/mylaila/transcriptomics_project](https://github.com/mylaila/transcriptomics_project)

## Vue d'ensemble

Ce projet analyse les données RNA-seq du projet TCGA-UCEC (The Cancer Genome Atlas - Uterine Corpus Endometrial Carcinoma) pour identifier les gènes différentiellement exprimés et caractériser les signatures moléculaires du cancer de l'endomètre.

**Cohorte à analyser** :
- Données TCGA-UCEC
- Échantillons tumoraux et normaux
- À définir selon les analyses QC

## Structure du projet

```
TCGA_UCEC_project/
├── data/                              # Données brutes et traitées
│   ├── exports/                       # Exports pour analyses
│   │   └── DESeq2_Results/
│   └── figures/                       # Figures générées
│
├── documentation/                     # Documentation du projet
│
├── Results_R_Analysis/                # Résultats de l'analyse différentielle
│
├── tmp_cache/                         # Fichiers temporaires
│
└── Notebooks et scripts principaux
    ├── 00_Exploratory_Analysis.ipynb  # Exploration initiale des données
    ├── 01_Reduction_QC_Cohort.ipynb   # QC et sélection de cohorte
    ├── 02_Single_Cell_Analysis.ipynb  # Analyse principale
    ├── 03_Pseudobulk_For_DE.ipynb     # Préparation pour analyse différentielle
    ├── 04_R_Analyse_DESeq2.R          # Analyse différentielle avec DESeq2
    └── 05_R_GO_Enrichissement.R       # Enrichissement fonctionnel
```

## Démarrage rapide

### 1. Activer l'environnement

**Windows (PowerShell/CMD)** :
```cmd
conda activate tcga_env
```

**VS Code** : Sélectionnez l'environnement Python approprié dans vos notebooks.

### 2. Vérifier l'installation

```python
python -c "import pandas, numpy; print('✅ Environnement prêt!')"
```

## Pipeline d'analyse

### Étape 1 : Exploration des données

**Script** : [00_Exploratory_Analysis.ipynb](00_Exploratory_Analysis.ipynb)

- Exploration initiale des données TCGA-UCEC
- Caractérisation de la cohorte
- Visualisations préliminaires

### Étape 2 : Contrôle qualité

**Script** : [01_Reduction_QC_Cohort.ipynb](01_Reduction_QC_Cohort.ipynb)

- Métriques de qualité
- Filtrage des échantillons
- Sélection de la cohorte finale

### Étape 3 : Analyse principale

**Script** : [02_Single_Cell_Analysis.ipynb](02_Single_Cell_Analysis.ipynb)

- Analyses exploratoires
- Identification de patterns
- Visualisations clés

### Étape 4 : Préparation pour analyse différentielle

**Script** : [03_Pseudobulk_For_DE.ipynb](03_Pseudobulk_For_DE.ipynb)

- Préparation des matrices de comptage
- Export pour analyse différentielle

### Étape 5 : Analyse différentielle avec DESeq2

**Script** : [04_R_Analyse_DESeq2.R](04_R_Analyse_DESeq2.R)

- Normalisation DESeq2
- Identification des gènes différentiellement exprimés
- Seuils : padj < 0.05, |log2FC| > 1

**Sortie** : `Results_R_Analysis/`

### Étape 6 : Enrichissement fonctionnel

**Script** : [05_R_GO_Enrichissement.R](05_R_GO_Enrichissement.R)

- Analyse d'enrichissement GO
- Visualisations des voies biologiques

## Environnement technique

**Python** (3.8+) :
- pandas
- numpy
- matplotlib
- seaborn
- scipy
- scikit-learn

**R** (4.0+) :
- DESeq2
- clusterProfiler
- org.Hs.eg.db
- ggplot2
- pheatmap

## Notes méthodologiques

### Source des données

Données RNA-seq du projet TCGA-UCEC (The Cancer Genome Atlas - Uterine Corpus Endometrial Carcinoma).

### Analyse différentielle

Utilisation de DESeq2 pour l'identification des gènes différentiellement exprimés entre échantillons tumoraux et normaux.

### Seuils statistiques

- **padj < 0.05** : Seuil de significativité (FDR)
- **|log2FC| > 1** : Fold-change minimum

## Structure des résultats

Les résultats sont organisés dans :
- `data/exports/DESeq2_Results/` : Résultats bruts des analyses
- `Results_R_Analysis/` : Tableaux de gènes différentiellement exprimés
- `data/figures/` : Visualisations générées

## Références

- **TCGA** : The Cancer Genome Atlas Program
- **DESeq2** : Love, Huber & Anders (2014)
- **clusterProfiler** : Yu et al. (2012)

---

**Auteur** : Projet M2 AIDA
**Date** : Décembre 2025
