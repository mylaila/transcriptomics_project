# Analyse transcriptomique du cortex préfrontal dans les maladies neurodégénératives

Projet M2 AIDA - Analyse différentielle de l'expression génique dans la maladie d'Alzheimer et la maladie de Parkinson

**Repository** : [https://github.com/mylaila/transcriptomics_project](https://github.com/mylaila/transcriptomics_project)

## Vue d'ensemble

Ce projet analyse les données single-nucleus RNA-seq du cortex préfrontal pour identifier les gènes différentiellement exprimés entre patients atteints d'Alzheimer (AD), Parkinson (PD) et contrôles sains.

**Cohorte analysée** :
- 17 donneurs (8 AD, 3 PD, 6 CTRL)
- 62 800 noyaux après QC
- 7 types cellulaires annotés

## Structure du projet

```
transcriptomics_project/
├── data/                              # Données brutes et traitées
│   ├── AD_PD_CTRL.h5ad               # Dataset initial
│   ├── adata_filtered.h5ad            # Post-QC
│   ├── adata_pp.h5ad                  # Normalisé + clustering
│   └── adata_annotated.h5ad           # Annotated
│
├── exports/                           # Données pseudobulk pour R
│   ├── pseudobulk_counts.csv
│   └── pseudobulk_metadata.csv
│
├── figures/                           # Figures générées (PNG 300dpi)
│   ├── Fig1_QC_ViolinPlot.png
│   ├── Fig2_UMAP_MultiPanel.png
│   ├── Fig3_DotPlot_Markers.png
│   ├── Fig4_Volcano_Strategic.png
│   └── Fig5_ORA_Barplot.png
│
├── Results_R_Analysis/                # Résultats de l'analyse différentielle
│   ├── Microglia_AD_vs_CTRL.csv
│   ├── Microglia_PD_vs_CTRL.csv
│   └── *_GO_results.csv
│
└── Notebooks et scripts principaux
    ├── reduction_of_dataset.ipynb     # S1 : QC et filtrage
    ├── single_cell_pipeline.ipynb     # S2 : Normalisation, clustering, UMAP
    ├── add_cell_type_annotation.ipynb # S3 : Annotation manuelle des clusters
    ├── translation_to_R.ipynb         # S4 : Agrégation pseudobulk
    ├── final_analysis_DGE.R           # S5 : Analyse différentielle (limma/voom)
    └── POST_R_FIGURES_4_ET_5.ipynb    # Visualisation finale
```

## Démarrage rapide

### 1. Activer l'environnement

**Windows (PowerShell/CMD)** :
```cmd
conda activate transcriptomics
```

Ou utilisez le script fourni :
```cmd
.\.vscode\activate.bat
```

**VS Code** : L'environnement est configuré automatiquement. Vérifiez que le kernel "Python (transcriptomics)" est sélectionné dans vos notebooks.

### 2. Vérifier l'installation

```python
python -c "import scanpy, anndata, pandas; print('✅ Environnement prêt!')"
```

## Pipeline d'analyse

### Étape 1 : Contrôle qualité et sélection de cohorte

**Script** : `reduction_of_dataset.ipynb`

- Chargement du dataset initial (6.28 GB, 26 pathologies)
- Métriques QC : n_genes, n_counts, pct_counts_mt
- Filtres appliqués :
  - Gènes : 1000-6000 par cellule
  - Mitochondrial : < 5%
  - Sélection cohorte : AD, PD, CTRL uniquement
  - Exclusion : Ascendance africaine, comorbidités

**Sortie** : `data/adata_filtered.h5ad` (62 800 cellules, 34 176 gènes)

### Étape 2 : Normalisation et clustering

**Script** : `single_cell_pipeline.ipynb`

- Normalisation (10 000 counts/cellule)
- Sélection HVG (2000 gènes)
- PCA (50 composantes)
- Clustering Leiden (23 clusters)
- UMAP
- Identification des marqueurs

**Sortie** : `data/adata_pp.h5ad`

### Étape 3 : Annotation des types cellulaires

**Script** : `add_cell_type_annotation.ipynb`

Types cellulaires identifiés :
- Neurones excitateurs (36.3%)
- Neurones inhibiteurs (20.7%)
- Oligodendrocytes/OPC (15.4%)
- Microglie (7.8%)
- Astrocytes homéostatiques (7.6%)
- Astrocytes réactifs (1.1%)
- Support/Vasculaire/Immunitaire (10.5%)

**Sortie** : `data/adata_annotated.h5ad`

### Étape 4 : Agrégation pseudobulk

**Script** : `translation_to_R.ipynb`

- Agrégation par (type cellulaire × donneur)
- Filtres : ≥20 cellules/échantillon, ≥5 donneurs/type
- Export pour analyse R

**Sortie** : `exports/pseudobulk_counts.csv`, `exports/pseudobulk_metadata.csv`

### Étape 5 : Analyse différentielle

**Script** : `final_analysis_DGE.R`

- Normalisation TMM
- Modèle linéaire (limma/voom)
- Contrastes : AD vs CTRL, PD vs CTRL
- Enrichissement GO (clusterProfiler)
- Seuils : |log2FC| > 0.5, FDR < 0.05

**Sortie** : `Results_R_Analysis/*.csv`

### Étape 6 : Visualisation finale

**Script** : `POST_R_FIGURES_4_ET_5.ipynb`

Génération des volcano plots et barplots d'enrichissement GO pour le rapport.

## Génération des figures

Les 5 figures obligatoires sont générées automatiquement :

| Figure | Description | Script |
|--------|-------------|--------|
| Fig 1 | QC Violin Plot (avant/après) | S1 |
| Fig 2 | UMAP Multi-Panel (4 colorations) | S2 |
| Fig 3 | Dot Plot des marqueurs | S2 |
| Fig 4 | Volcano Plot (AD + PD) | Python post-R |
| Fig 5 | ORA Barplot | Python post-R |

**Pour générer toutes les figures** :

1. Exécuter les notebooks S1 à S4 dans l'ordre
2. Exécuter le script R : `Rscript final_analysis_DGE.R`
3. Exécuter le notebook Python post-R

Toutes les figures sont sauvegardées en PNG 300 dpi dans `figures/`.

## Environnement technique

**Activation de l'environnement** :
```bash
conda activate transcriptomics
```

**Python** (3.8+) :
- scanpy
- anndata
- pandas
- numpy
- matplotlib
- seaborn

**R** (4.0+) :
- limma
- edgeR
- clusterProfiler
- org.Hs.eg.db
- ggplot2

## Notes méthodologiques

### Format des identifiants de gènes

Les gènes sont au format **ENSEMBL IDs** (ENSG...). Le script R a été configuré pour utiliser `keyType = "ENSEMBL"` dans les analyses d'enrichissement.

### Types cellulaires analysés

Focus sur la **Microglie** pour les Figures 4 & 5, type cellulaire particulièrement pertinent dans les pathologies neurodégénératives (réponse inflammatoire, phagocytose).

### Limitations

- Pas de correction de batch (petite cohorte, surtout PD)
- Analyse exploratoire (validation fonctionnelle nécessaire)
- Dataset limité au cortex préfrontal

## Références

- Dataset source : CZ CELLxGENE Discover
- Pipeline basé sur : scanpy best practices (Luecken & Theis, 2019)
- Analyse différentielle : limma/voom (Law et al., 2014)

---

**Auteur** : Projet M2 AIDA
**Date** : Décembre 2025
