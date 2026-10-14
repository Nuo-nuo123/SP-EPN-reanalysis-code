# Step 5: DMR (Differentially Methylated Regions) + Multi-omics Integration
# Keep this file pure ASCII.

Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp",
           TEMP="D:/claude-study/Rtmp",
           TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
sink(file.path(outdir, "step5_dmr_integration.log"), append=FALSE)
cat("== step5 start", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
flush.console()

# ============================================================
# 5a: Load methylation beta values
# ============================================================
cat("\n--- 5a: Load methylation data ---\n")
flush.console()

# Use the pre-processed beta matrix (GSE242991)
beta_file <- "D:/claude-study/Spinal Ependymoma/GEO_data/GSE242991/GSE242991_20231122_beta_values_normalized.csv.gz"
if (file.exists(beta_file)) {
  cat("Loading beta matrix from GSE242991...\n")
  beta <- fread(beta_file, nrows=10)
  cat("Dimensions:", dim(beta), "\n")
  # Read full matrix
  beta <- fread(beta_file)
  cat("Full beta matrix:", dim(beta), "\n")
} else {
  cat("Beta matrix not found, checking pval matrix...\n")
  # Fallback: use pval processed matrix from GSE243240
  pval_file <- "D:/claude-study/Spinal Ependymoma/GEO_data/GSE243240/GSE243240_pval_processed_matrix_070923.txt.gz"
  if (file.exists(pval_file)) {
    cat("Loading pval matrix...\n")
    pval <- fread(pval_file, nrows=10)
    cat("Pval matrix preview:", dim(pval), "\n")
  }
  beta <- NULL
}

# ============================================================
# 5b: Load RNA-seq labels for subtype assignment
# ============================================================
cat("\n--- 5b: Load sample labels ---\n")
flush.console()

labels_file <- file.path(outdir, "GSE242994_rna_labels.csv")
if (file.exists(labels_file)) {
  labels <- fread(labels_file)
  cat("Labels loaded:", nrow(labels), "samples\n")
  cat("Subtype distribution:\n")
  print(table(labels$Transcriptional_subtype))
} else {
  cat("Labels file not found\n")
  labels <- NULL
}

# ============================================================
# 5c: DEG analysis (bulk RNA-seq: Subtype B vs A)
# ============================================================
cat("\n--- 5c: Bulk DEG analysis ---\n")
flush.console()

expr_file <- file.path(outdir, "GSE242994_expr.rds")
if (file.exists(expr_file)) {
  expr_raw <- readRDS(expr_file)
  if (is.list(expr_raw) && "expr" %in% names(expr_raw)) {
    expr <- expr_raw$expr
  } else {
    expr <- expr_raw
  }
  cat("Expression matrix:", dim(expr), "\n")

  if (!is.null(labels)) {
    common <- intersect(colnames(expr), labels$sample)
    cat("Common samples:", length(common), "\n")

    if (length(common) > 10) {
      # Simple DEG: Wilcoxon test per gene
      subtype <- labels$Transcriptional_subtype[match(common, labels$sample)]
      group_A <- common[subtype == "A"]
      group_B <- common[subtype == "B"]
      cat("Group A:", length(group_A), "Group B:", length(group_B), "\n")

      deg_results <- data.frame(gene=rownames(expr), stringsAsFactors=FALSE)
      pvals <- numeric(nrow(expr))
      logFC <- numeric(nrow(expr))

      for (i in seq_len(nrow(expr))) {
        vals_A <- as.numeric(expr[i, group_A])
        vals_B <- as.numeric(expr[i, group_B])
        tryCatch({
          wt <- wilcox.test(vals_B, vals_A, exact=FALSE)
          pvals[i] <- wt$p.value
          logFC[i] <- mean(vals_B) - mean(vals_A)
        }, error=function(e) {
          pvals[i] <- 1
          logFC[i] <- 0
        })
        if (i %% 1000 == 0) {
          cat("  DEG: gene", i, "/", nrow(expr), "\n")
          flush.console()
        }
      }

      deg_results$pvalue <- pvals
      deg_results$logFC <- logFC
      deg_results$padj <- p.adjust(pvals, method="BH")
      deg_results <- deg_results[order(deg_results$pvalue), ]

      sig_deg <- deg_results[deg_results$padj < 0.05 & abs(deg_results$logFC) > 0.5, ]
      cat("Significant DEGs (padj<0.05, |logFC|>0.5):", nrow(sig_deg), "\n")

      write.csv(deg_results, file.path(outdir, "bulk_deg_all.csv"), row.names=FALSE)
      write.csv(sig_deg, file.path(outdir, "bulk_deg_significant.csv"), row.names=FALSE)
    }
  }
} else {
  cat("Expression matrix not found\n")
}

# ============================================================
# 5d: Multi-omics integration
# ============================================================
cat("\n--- 5d: Multi-omics integration ---\n")
flush.console()

# Collect gene lists from different sources
gene_lists <- list()

# 1. Bulk DEGs
if (exists("sig_deg") && nrow(sig_deg) > 0) {
  gene_lists$bulk_DEG <- sig_deg$gene
  cat("Bulk DEGs:", length(gene_lists$bulk_DEG), "\n")
}

# 2. SCENIC regulon targets (if available)
scenic_file <- "D:/claude-study/Spinal Ependymoma/analysis/GSE298501/scenic_results.rds"
if (file.exists(scenic_file)) {
  scenic <- readRDS(scenic_file)
  if ("regulons" %in% names(scenic)) {
    gene_lists$scenic_targets <- unique(unlist(scenic$regulons))
    cat("SCENIC targets:", length(gene_lists$scenic_targets), "\n")
  }
}

# 3. Cluster markers from scRNA
markers_file <- "D:/claude-study/Spinal Ependymoma/analysis/GSE298501/cluster_markers.csv"
if (file.exists(markers_file)) {
  markers <- fread(markers_file)
  gene_lists$sc_markers <- unique(markers$gene)
  cat("scRNA markers:", length(gene_lists$sc_markers), "\n")
}

# Integration: genes appearing in >= 2 sources
if (length(gene_lists) >= 2) {
  all_genes <- unique(unlist(gene_lists))
  gene_counts <- sapply(all_genes, function(g) sum(sapply(gene_lists, function(x) g %in% x)))

  scenic_genes <- if (!is.null(gene_lists$scenic_targets)) gene_lists$scenic_targets else character(0)
  marker_genes <- if (!is.null(gene_lists$sc_markers)) gene_lists$sc_markers else character(0)
  bulk_genes <- if (!is.null(gene_lists$bulk_DEG)) gene_lists$bulk_DEG else character(0)

  multi_omics <- data.frame(
    gene=all_genes,
    n_sources=gene_counts,
    in_bulk_DEG=all_genes %in% bulk_genes,
    in_scenic=all_genes %in% scenic_genes,
    in_sc_markers=all_genes %in% marker_genes,
    stringsAsFactors=FALSE
  )
  multi_omics <- multi_omics[order(-multi_omics$n_sources), ]

  integrated_genes <- multi_omics$gene[multi_omics$n_sources >= 2]
  cat("Genes in >= 2 sources:", length(integrated_genes), "\n")

  write.csv(multi_omics, file.path(outdir, "multiomics_integration.csv"), row.names=FALSE)
  cat("Multi-omics integration table saved\n")
} else {
  cat("Not enough gene lists for integration\n")
}

cat("\n== step5 done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
