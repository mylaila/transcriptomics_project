# ==============================================================================
# ANALYSE DIFFERENTIELLE SCZ vs CTRL - DESeq2 (VERSION CORRIGEE)
# Projet : M2 AIDA - Transcriptomique HBCC
# ==============================================================================

# 1. CHARGEMENT DES BIBLIOTHEQUES
# ==============================================================================
cat("Chargement des bibliotheques...\n")
library(DESeq2)
library(tidyverse)

cat("\n======================================================================\n")
cat("ANALYSE DIFFERENTIELLE SCZ vs CTRL\n")
cat("======================================================================\n\n")

# 2. CONFIGURATION
# ==============================================================================
base_dir <- "C:/Z/M2_AIDA/transcriptomics_project/data/exports"

counts_file <- file.path(base_dir, "pseudobulk_counts.csv")
metadata_file <- file.path(base_dir, "pseudobulk_metadata.csv")

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

# 4. ANALYSE PAR TYPE CELLULAIRE
# ==============================================================================
results_list <- list()
summary_table <- data.frame()

cell_types <- unique(metadata$cell_type)
cat("\nNombre de types cellulaires a analyser :", length(cell_types), "\n\n")

for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    
    cat("\n==================================================================\n")
    cat("[", i, "/", length(cell_types), "] Type cellulaire :", cell_type, "\n")
    cat("==================================================================\n")
    
    # Subset
    idx <- metadata$cell_type == cell_type
    meta_sub <- metadata[idx, ]
    counts_sub <- counts[, idx]
    
    # Verifier les replicats
    n_ctrl <- sum(meta_sub$disease == "normal", na.rm = TRUE)
    n_scz <- sum(meta_sub$disease == "schizophrenia", na.rm = TRUE)
    
    cat("Echantillons : CTRL =", n_ctrl, ", SCZ =", n_scz, "\n")
    
    if (n_ctrl < 3 || n_scz < 3) {
        cat("SKIP : Replicats insuffisants\n")
        next
    }
    
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
        
        # Extraire les resultats
        res <- results(dds, contrast = c("disease", "schizophrenia", "normal"))
        res_ordered <- res[order(res$pvalue), ]
        
        # Stocker
        results_list[[cell_type]] <- as.data.frame(res_ordered)
        
        # Statistiques
        n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
        n_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
        n_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
        
        cat("Resultats :\n")
        cat("   Genes testes :", nrow(res_ordered), "\n")
        cat("   DEGs (FDR < 0.05) :", n_sig, "\n")
        cat("   Up-regules (SCZ > CTRL) :", n_up, "\n")
        cat("   Down-regules (SCZ < CTRL) :", n_down, "\n")
        
        # Stocker dans le tableau recapitulatif
        summary_table <- rbind(summary_table, data.frame(
            cell_type = cell_type,
            n_samples_ctrl = n_ctrl,
            n_samples_scz = n_scz,
            n_genes_tested = nrow(res_ordered),
            n_DEGs = n_sig,
            n_up = n_up,
            n_down = n_down,
            stringsAsFactors = FALSE
        ))
        
    }, error = function(e) {
        cat("ERREUR pour", cell_type, ":", e$message, "\n")
    })
}

# 5. SAUVEGARDER LES RESULTATS
# ==============================================================================
cat("\n\nSauvegarde des resultats...\n")

output_dir <- file.path(base_dir, "DESeq2_Results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Sauvegarder chaque fichier
for (cell_type in names(results_list)) {
    # Nettoyer le nom pour eviter les problemes de chemin
    clean_name <- gsub("/", "_", cell_type)
    clean_name <- gsub(" ", "_", clean_name)
    
    filename <- file.path(output_dir, paste0(clean_name, "_DESeq2.csv"))
    write.csv(results_list[[cell_type]], filename, row.names = TRUE)
    cat("   ", basename(filename), "\n")
}

# Sauvegarder le tableau recapitulatif
summary_file <- file.path(output_dir, "00_SUMMARY.csv")
write.csv(summary_table, summary_file, row.names = FALSE)
cat("   00_SUMMARY.csv\n")

# 6. RESUME FINAL
# ==============================================================================
cat("\n======================================================================\n")
cat("ANALYSE TERMINEE\n")
cat("======================================================================\n")
cat("Resultats dans :", output_dir, "\n")
cat("Types cellulaires analyses :", nrow(summary_table), "\n")
cat("Total DEGs (tous types confondus) :", sum(summary_table$n_DEGs), "\n\n")

# Afficher le tableau recapitulatif
cat("TABLEAU RECAPITULATIF :\n")
print(summary_table)

cat("\n======================================================================\n")
cat("Script termine avec succes !\n")
cat("======================================================================\n")

# Sauvegarder aussi l'environnement pour analyse ulterieure
save(results_list, summary_table, file = file.path(output_dir, "DESeq2_workspace.RData"))
cat("\nWorkspace sauvegarde dans : DESeq2_workspace.RData\n")
