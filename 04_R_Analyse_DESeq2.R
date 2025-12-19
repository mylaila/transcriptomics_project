# ==============================================================================
# ANALYSE DIFFERENTIELLE TCGA-UCEC - DESeq2
# Projet : M2 AIDA - Cancer de l'endomètre
# ==============================================================================

# 1. CHARGEMENT DES BIBLIOTHEQUES
# ==============================================================================
cat("Chargement des bibliotheques...\n")
library(DESeq2)
library(tidyverse)

cat("\n======================================================================\n")
cat("ANALYSE DIFFERENTIELLE TCGA-UCEC - TUMEUR vs NORMAL\n")
cat("======================================================================\n\n")

# 2. CONFIGURATION
# ==============================================================================
base_dir <- "C:/Z/M2_AIDA/TCGA_UCEC_project/data/exports"

counts_file <- file.path(base_dir, "expression_counts.csv")
metadata_file <- file.path(base_dir, "clinical_metadata.csv")

cat("Verification des fichiers...\n")
if (!file.exists(counts_file)) stop("Fichier counts introuvable")
if (!file.exists(metadata_file)) stop("Fichier metadata introuvable")
cat("Fichiers detectes\n\n")

# 3. CHARGEMENT DES DONNEES
# ==============================================================================
cat("Chargement des donnees...\n")

counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)

# NETTOYAGE DE LA COLONNE DISEASE
metadata$disease <- trimws(metadata$disease)
metadata$disease <- factor(metadata$disease, levels = c("normal", "schizophrenia"))

cat("Donnees chargees\n")
cat("   Matrice counts :", nrow(counts), "genes x", ncol(counts), "echantillons\n")
cat("   Metadonnees :", nrow(metadata), "echantillons\n")
cat("   Disease distribution :\n")
print(table(metadata$disease))

# Verifier la correspondance
if (!all(colnames(counts) == rownames(metadata))) {
    stop("ERREUR : Colonnes counts != lignes metadata")
}

# 4. ANALYSE DIFFERENTIELLE
# ==============================================================================
cat("\nPreparation de l'analyse differentielle...\n")

# TODO: Adapter selon la structure de vos métadonnées
# Exemple : variable 'sample_type' contenant 'Tumor' et 'Normal'

# Verifier les conditions
# condition_column <- "sample_type"  # À adapter
# levels <- unique(metadata[[condition_column]])
# cat("Conditions detectees :", paste(levels, collapse = ", "), "\n")

# n_tumor <- sum(metadata[[condition_column]] == "Tumor", na.rm = TRUE)
# n_normal <- sum(metadata[[condition_column]] == "Normal", na.rm = TRUE)
# cat("Echantillons : Tumor =", n_tumor, ", Normal =", n_normal, "\n")

# if (n_normal < 3 || n_tumor < 3) {
#     stop("Replicats insuffisants pour l'analyse DESeq2")
# }
    
    # Convertir en matrice d'entiers
    counts_matrix <- as.matrix(counts_sub)
    counts_matrix <- round(counts_matrix)
    
    # S'assurer que disease est un facteur
    meta_sub$disease <- factor(meta_sub$disease, levels = c("normal", "schizophrenia"))
    
    # Creer l'objet DESeq2
    tryCatch({
        dds <- DESeqDataSetFromMatrix(
            countData = counts_matrix,
            colData = meta_sub,
            design = ~ disease
        )
        
        # Filtrer genes peu exprimes
        keep <- rowSums(counts(dds) >= 10) >= 3
        dds <- dds[keep, ]
        
        cat("Genes retenus apres filtrage :", sum(keep), "/", length(keep), "\n")
        
        if (sum(keep) < 100) {
            cat("SKIP : Trop peu de genes exprimes\n")
            next
        }
        
        # Lancer DESeq2
        cat("Execution de DESeq2...\n")
        dds <- DESeq(dds, quiet = TRUE)
        
# TODO: Adapter le code ci-dessous selon vos données

# Exemple de workflow DESeq2 (à décommenter et adapter):
#
# # Convertir en matrice d'entiers
# counts_matrix <- as.matrix(counts)
# counts_matrix <- round(counts_matrix)
# 
# # S'assurer que la condition est un facteur
# metadata$sample_type <- factor(metadata$sample_type, levels = c("Normal", "Tumor"))
# 
# # Creer l'objet DESeq2
# tryCatch({
#     dds <- DESeqDataSetFromMatrix(
#         countData = counts_matrix,
#         colData = metadata,
#         design = ~ sample_type
#     )
#     
#     # Filtrer genes peu exprimes
#     keep <- rowSums(counts(dds) >= 10) >= 3
#     dds <- dds[keep, ]
#     cat("Genes retenus apres filtrage :", sum(keep), "/", length(keep), "\n")
#     
#     # Analyse DESeq2
#     dds <- DESeq(dds)
#     
#     # Extraire les resultats
#     res <- results(dds, contrast = c("sample_type", "Tumor", "Normal"))
#     res_ordered <- res[order(res$pvalue), ]
#     
#     # Statistiques
#     n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
#     n_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
#     n_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
#     
#     cat("Resultats :\n")
#     cat("   Genes testes :", nrow(res_ordered), "\n")
#     cat("   DEGs (padj < 0.05) :", n_sig, "\n")
#     cat("   Up-regules (Tumor > Normal) :", n_up, "\n")
#     cat("   Down-regules (Tumor < Normal) :", n_down, "\n")
#     
# }, error = function(e) {
#     cat("ERREUR :", e$message, "\n")
# })

# 5. SAUVEGARDER LES RESULTATS
# ==============================================================================
# TODO: Décommenter et adapter selon vos résultats

# cat("\n\nSauvegarde des resultats...\n")
# 
# output_dir <- "C:/Z/M2_AIDA/TCGA_UCEC_project/Results_R_Analysis"
# dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# 
# # Sauvegarder les résultats
# results_file <- file.path(output_dir, "TCGA_UCEC_DESeq2_results.csv")
# write.csv(as.data.frame(res_ordered), results_file, row.names = TRUE)
# cat("Resultats sauvegardes dans :", results_file, "\n")
# 
# # Sauvegarder l'environnement
# save.image(file = file.path(output_dir, "DESeq2_workspace.RData"))
# cat("\nWorkspace sauvegarde dans : DESeq2_workspace.RData\n")

cat("\n======================================================================\n")
cat("SCRIPT TEMPLATE - TCGA-UCEC DESeq2\n")
cat("======================================================================\n")
cat("Ce script doit etre adapte selon la structure de vos donnees.\n")
cat("Consultez la documentation DESeq2 pour plus de details.\n")
cat("======================================================================\n")
