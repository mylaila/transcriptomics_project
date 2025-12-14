library(clusterProfiler)
library(org.Hs.eg.db)

# Charger les DEGs VIP (374 gènes)
vip_degs <- read.csv("C:/Z/M2_AIDA/transcriptomics_project/data/exports/DESeq2_Results/VIP_GABAergic_cortical_interneuron_DESeq2.csv")
sig_genes <- vip_degs[vip_degs$padj < 0.05, ]

# Enrichissement GO
go_results <- enrichGO(
    gene = rownames(sig_genes),
    OrgDb = org.Hs.eg.db,
    ont = "BP",  # Biological Process
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
)

# Top 10 termes
dotplot(
  go_results,
  showCategory = 10,
  orderBy = "x",
  title = "Enrichissement des processus biologiques (GO) dans les interneurones GABAergiques VIP – SCH vs CTRL"
)
