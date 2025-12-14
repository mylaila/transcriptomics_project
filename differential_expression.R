# Charger les bibliothèques nécessaires
library(HDF5Array)
library(DESeq2)

# Définir le chemin d'accès au fichier .h5ad
file_path <- "C:/Z/M2_AIDA/transcriptomics_project/data/HBCC_SCZ_CTRL_postQC_UMAP.h5ad"

# Charger les données en utilisant HDF5Array
adata <- HDF5Array::loadHDF5(file_path)

# Vérifier les informations de base sur le fichier
print(adata)

# Extraire les données d'expression (comptages des cellules)
expression_data <- adata$X

# Extraire les métadonnées des cellules (disease, donor_id)
meta_data <- adata$obs

# Vérifier les premières lignes des métadonnées pour s'assurer que les informations sont correctes
head(meta_data)

# Assurer que la colonne 'disease' contient SCZ et CTRL, et la convertir en facteur
meta_data$disease <- factor(meta_data$disease, levels = c("CTRL", "SCZ"))

# Créer un DataFrame pour DESeq2
dds <- DESeqDataSetFromMatrix(countData = expression_data,
                              colData = meta_data,
                              design = ~ disease)  # Comparaison entre SCZ vs CTRL

# Vérifier le contenu du DESeqDataSet
print(dds)

# Effectuer l'analyse différentielle
dds <- DESeq(dds)

# Extraire les résultats de l'analyse
res <- results(dds)

# Résumer les résultats
summary(res)

# Visualiser les résultats avec un plot MA
plotMA(res, main = "Analyse MA")

# Afficher les 10 premiers gènes différentiellement exprimés (si p-value ajustée < 0.05)
DE_genes <- res[res$padj < 0.05, ]
print(head(DE_genes))

# Exporter les résultats dans un fichier CSV
write.csv(as.data.frame(res), "differential_expression_results.csv")
