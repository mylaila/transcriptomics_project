# ==============================================================================
# ENRICHISSEMENT FONCTIONNEL GO - TCGA-UCEC
# Projet : M2 AIDA - Cancer de l'endomètre
# ==============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

cat("======================================================================\n")
cat("ENRICHISSEMENT GO - TCGA-UCEC\n")
cat("======================================================================\n\n")

# TODO: Adapter selon vos résultats DESeq2

# Exemple de workflow (à décommenter et adapter):
#
# # Charger les résultats DESeq2
# degs <- read.csv("C:/Z/M2_AIDA/TCGA_UCEC_project/Results_R_Analysis/TCGA_UCEC_DESeq2_results.csv", 
#                  row.names = 1)
# 
# # Filtrer les gènes significatifs
# sig_genes <- degs[degs$padj < 0.05 & abs(degs$log2FoldChange) > 1, ]
# 
# cat("Nombre de DEGs significatifs :", nrow(sig_genes), "\n\n")
# 
# # Enrichissement GO
# # Note: Adapter le format des IDs de gènes (ENSEMBL, SYMBOL, ENTREZ, etc.)
# go_results <- enrichGO(
#     gene = rownames(sig_genes),
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENSEMBL",  # ou "SYMBOL" selon vos données
#     ont = "BP",           # Biological Process
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05
# )
# 
# # Visualisation Top 10 termes
# p <- dotplot(
#   go_results,
#   showCategory = 10,
#   orderBy = "x",
#   title = "Enrichissement GO - TCGA-UCEC (Tumor vs Normal)"
# )
# 
# ggsave("C:/Z/M2_AIDA/TCGA_UCEC_project/data/figures/GO_enrichment.png", 
#        p, width = 10, height = 8, dpi = 300)

cat("\n======================================================================\n")
cat("SCRIPT TEMPLATE - Enrichissement GO\n")
cat("======================================================================\n")
cat("Ce script doit etre adapte selon vos resultats DESeq2.\n")
cat("Assurez-vous d'utiliser le bon format d'identifiants de genes.\n")
cat("======================================================================\n")
