# ==============================================================================
# SCRIPT: Differential Gene Expression (DGE) & Pathway Enrichment Analysis
# PROJECT: M2 AIDA - Transcriptomic Analysis of the Prefrontal Cortex
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. ENVIRONMENT SETUP & DATA LOADING
# ------------------------------------------------------------------------------

rm(list = ls())
graphics.off()

cat("======================================================================\n")
cat("🚀 STARTING R ANALYSIS PIPELINE: DGE & ENRICHMENT\n")
cat("======================================================================\n")

# --- A. Locate Export Directory ---
input_dir <- "exports"
if (!dir.exists(input_dir) && dir.exists("export")) { input_dir <- "export" }

if (!dir.exists(input_dir)) {
  if (file.exists("pseudobulk_counts.csv")) {
    input_dir <- "."
  } else {
    stop("CRITICAL ERROR: 'exports' directory not found.")
  }
}

print(paste("📂 Data Directory identified:", input_dir))

# --- B. Load Data ---
file_counts <- file.path(input_dir, "pseudobulk_counts.csv")
file_meta   <- file.path(input_dir, "pseudobulk_metadata.csv")

if (!file.exists(file_counts)) stop("ERROR: 'pseudobulk_counts.csv' is missing.")
if (!file.exists(file_meta))   stop("ERROR: 'pseudobulk_metadata.csv' is missing.")

counts <- read.csv(file_counts, row.names = 1, check.names = FALSE)
metadata <- read.csv(file_meta, row.names = 1)

cat("✅ Data loaded successfully.\n")

if(!all(colnames(counts) == rownames(metadata))) {
  stop("DATA INTEGRITY ERROR: Count matrix columns do not match metadata rows.")
}

# --- C. RECODING DISEASE LABELS (CRITICAL FIX) ---
# Mapping long names to short codes for analysis
cat("🔄 Recoding disease labels to standard codes (AD, PD, CTRL)...\n")

metadata$disease[metadata$disease == "dementia || Alzheimer disease"] <- "AD"
metadata$disease[metadata$disease == "dementia || Parkinson disease"] <- "PD"
metadata$disease[metadata$disease == "normal"] <- "CTRL"

# Verify recoding
print(table(metadata$disease))

# ------------------------------------------------------------------------------
# 2. LIBRARY LOADING
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 3. ANALYSIS CONFIGURATION
# ------------------------------------------------------------------------------

colnames(metadata) <- make.names(colnames(metadata))
output_dir <- "Results_R_Analysis"
dir.create(output_dir, showWarnings = FALSE)

cell_types <- unique(metadata$cell_type_annotation)
cat(paste("\n🔬 Cell Types to analyze:", length(cell_types), "\n"))

# ------------------------------------------------------------------------------
# 4. MAIN LOOP
# ------------------------------------------------------------------------------

for (ct in cell_types) {
  
  cat(paste0("\n--------------------------------------------------\n"))
  cat(paste0("Processing: ", ct, "\n"))
  
  # --- A. Subset ---
  samples_keep <- rownames(metadata)[metadata$cell_type_annotation == ct]
  
  if(length(samples_keep) < 6) {
    cat("⚠️  Skipping: Insufficient sample size (< 6 samples).\n")
    next
  }
  
  curr_counts <- counts[, samples_keep]
  curr_meta <- metadata[samples_keep, ]
  
  # --- B. Filter ---
  keep_genes <- rowSums(cpm(curr_counts) > 1) >= (length(samples_keep) * 0.3)
  curr_counts <- curr_counts[keep_genes, ]
  
  # --- C. Norm & Model ---
  dge <- DGEList(counts = curr_counts)
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Use the RECODED 'disease' column
  group <- factor(curr_meta$disease)
  
  # Ensure CTRL is reference
  if("CTRL" %in% levels(group)) {
    group <- relevel(group, ref = "CTRL")
  } else {
    cat("⚠️  Warning: No 'CTRL' group found for this cell type. Skipping.\n")
    next
  }
  
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  
  # Contrasts
  contrasts_list <- c()
  if("AD" %in% colnames(design) & "CTRL" %in% colnames(design)) contrasts_list <- c(contrasts_list, "AD - CTRL")
  if("PD" %in% colnames(design) & "CTRL" %in% colnames(design)) contrasts_list <- c(contrasts_list, "PD - CTRL")
  
  if(length(contrasts_list) == 0) { 
    cat("⚠️  Skipping: No valid contrasts (AD/PD vs CTRL).\n")
    next 
  }
  
  cm <- makeContrasts(contrasts = contrasts_list, levels = design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  
  # --- D. Results ---
  
  for (comp_raw in contrasts_list) {
    comp_name <- gsub(" - ", "_vs_", comp_raw)
    top_res <- topTable(fit2, coef = comp_raw, number = Inf, sort.by = "P")
    
    filename_base <- file.path(output_dir, paste0(make.names(ct), "_", comp_name))
    write.csv(top_res, paste0(filename_base, ".csv"))
    
    n_sig <- sum(top_res$adj.P.Val < 0.05 & abs(top_res$logFC) > 0.5)
    cat(paste("   -> comparison", comp_name, ":", n_sig, "DEGs found.\n"))
    
    # Volcano Plot
    try({
      df_plot <- top_res
      df_plot$Significance <- "NS"
      df_plot$Significance[df_plot$adj.P.Val < 0.05 & df_plot$logFC > 0.5] <- "UP"
      df_plot$Significance[df_plot$adj.P.Val < 0.05 & df_plot$logFC < -0.5] <- "DOWN"
      
      p_vol <- ggplot(df_plot, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
        geom_point(alpha = 0.6, size = 1.5) +
        scale_color_manual(values = c("DOWN" = "blue", "NS" = "grey", "UP" = "red")) +
        theme_minimal() +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        labs(title = paste(ct, comp_name),
             subtitle = paste("DEGs:", n_sig),
             x = "Log2 Fold Change", y = "-Log10 FDR") +
        theme(legend.position = "top")
      
      ggsave(paste0(filename_base, "_Volcano.pdf"), p_vol, width = 6, height = 5)
    }, silent = TRUE)
    
    # Enrichment
    sig_genes <- rownames(top_res)[top_res$adj.P.Val < 0.05 & abs(top_res$logFC) > 0.5]
    
    if (length(sig_genes) >= 5) {
      cat("      -> Running GO Enrichment...\n")
      try({
        ego <- enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                        ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
        
        if (!is.null(ego) && nrow(ego) > 0) {
          p_dot <- dotplot(ego, showCategory=15) + ggtitle(paste("GO:", ct, comp_name))
          ggsave(paste0(filename_base, "_GO_Dotplot.pdf"), p_dot, width = 10, height = 7)
        }
      }, silent = TRUE)
    }
  }
}

cat("\n======================================================================\n")
cat("✅ PIPELINE COMPLETE. Results saved in 'Results_R_Analysis' folder.\n")
cat("======================================================================\n")