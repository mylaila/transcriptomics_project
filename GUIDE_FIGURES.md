# Guide de génération des figures pour le rapport

## Vue d'ensemble

Ce guide documente la procédure pour générer les 5 figures obligatoires du projet d'analyse transcriptomique (Alzheimer / Parkinson).

## Liste des figures

| Figure | Description | Script source | Format |
|--------|-------------|---------------|--------|
| **Figure 1** | QC Violin Plot (avant/après filtrage) | `reduction_of_dataset.ipynb` | PNG 300 dpi |
| **Figure 2** | UMAP Multi-Panel (clusters, diagnostic, sexe, ascendance) | `single_cell_pipeline.ipynb` | PNG 300 dpi |
| **Figure 3** | Dot Plot des marqueurs par cluster | `single_cell_pipeline.ipynb` | PNG 300 dpi |
| **Figure 4** | Volcano Plot (gènes différentiellement exprimés) | `POST_R_FIGURES_4_ET_5.ipynb` | PNG 300 dpi |
| **Figure 5** | Barplot d'enrichissement GO | `POST_R_FIGURES_4_ET_5.ipynb` | PNG 300 dpi |

## Procédure de génération

### Étape 1 : Configuration du script R

Le script `final_analysis_DGE.R` doit être configuré pour :
- Utiliser les identifiants ENSEMBL (format des gènes dans le dataset)
- Exporter les résultats d'enrichissement GO en CSV

**Modifications nécessaires dans `final_analysis_DGE.R`** :

**Ligne ~178** - Format des identifiants de gènes :
```r
ego <- enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
```

**Après ligne ~184** - Export des résultats GO :
```r
if (!is.null(ego) && nrow(ego) > 0) {
  ego_df <- as.data.frame(ego)
  ora_filename <- paste0(filename_base, "_GO_results.csv")
  write.csv(ego_df, ora_filename, row.names = FALSE)
}
```

### Étape 2 : Exécution du pipeline d'analyse

**2.1 Analyse différentielle (R)**
```bash
Rscript final_analysis_DGE.R
```

Durée estimée : 5-10 minutes

**2.2 Génération des figures Python**

Exécuter les notebooks dans l'ordre :
1. `single_cell_pipeline.ipynb` → Génère Figures 2 et 3
2. `POST_R_FIGURES_4_ET_5.ipynb` → Génère Figures 4 et 5

### Étape 3 : Vérification des sorties

Toutes les figures doivent être présentes dans le dossier `figures/` :
- `Fig1_QC_ViolinPlot.png`
- `Fig2_UMAP_MultiPanel.png`
- `Fig3_DotPlot_Markers.png`
- `Fig4_Volcano_Strategic.png`
- `Fig5_ORA_Barplot.png`

## Notes techniques

### Format des données

- **Identifiants de gènes** : ENSEMBL IDs (format ENSG...)
- **Type cellulaire analysé** : Microglia (pertinent pour les pathologies neurodégénératives)
- **Seuils d'expression différentielle** : |log2FC| > 0.5, FDR < 0.05

### Chemins des fichiers

Les résultats R sont stockés dans `Results_R_Analysis/` :
- Résultats DGE : `Microglia_AD_vs_CTRL.csv`, `Microglia_PD_vs_CTRL.csv`
- Résultats GO : `Microglia_AD_vs_CTRL_GO_results.csv`, `Microglia_PD_vs_CTRL_GO_results.csv`

### Résolution des erreurs courantes

| Erreur | Solution |
|--------|----------|
| `KeyError` sur une colonne | Vérifier que le notebook précédent a bien été exécuté |
| Fichiers CSV introuvables | Re-exécuter le script R |
| Images de mauvaise qualité | Vérifier le paramètre `dpi=300` |

## Structure du projet

```
transcriptomics-code/
├── data/                          # Données d'entrée
├── figures/                       # Figures générées (PNG)
├── Results_R_Analysis/            # Résultats de l'analyse R
├── reduction_of_dataset.ipynb     # S1 : QC et filtrage
├── single_cell_pipeline.ipynb     # S2 : Clustering et UMAP
├── add_cell_type_annotation.ipynb # S3 : Annotation cellulaire
├── translation_to_R.ipynb         # S4 : Préparation pseudobulk
├── final_analysis_DGE.R           # S5 : Analyse différentielle
└── POST_R_FIGURES_4_ET_5.ipynb    # Visualisation finale
```
