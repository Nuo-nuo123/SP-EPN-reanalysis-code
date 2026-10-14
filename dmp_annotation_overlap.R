# DMP annotation + overlap with published DMPs + gene mapping
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({library(data.table)})
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "dmp_annotation_overlap.log"), append=FALSE)
cat("== DMP annotation start", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# ============================================================
# 1. Load my DMP results
# ============================================================
dmp <- read.csv(file.path(outdir, "dmp_results.csv"), stringsAsFactors=FALSE)
dmp_sig <- dmp[!is.na(dmp$padj) & dmp$padj < 0.05, ]
cat("My significant DMPs:", nrow(dmp_sig), "\n")
dmp_sig_strong <- dmp_sig[abs(dmp_sig$delta_beta) > 0.2, ]
cat("My strong DMPs (|dB|>0.2):", nrow(dmp_sig_strong), "\n")

# ============================================================
# 2. Load published DMP table (pre-extracted by Python)
# ============================================================
pub_dmp_csv <- file.path(outdir, "published_dmp.csv")
stopifnot(file.exists(pub_dmp_csv))
pub <- read.csv(pub_dmp_csv, stringsAsFactors=FALSE, check.names=FALSE)
cat("Published DMP table rows:", nrow(pub), "\n")
pub_probes <- pub[[1]]
cat("Published probe count:", length(pub_probes), "\n")

# Overlap
ov <- intersect(dmp_sig$cpg, pub_probes)
cat("\nOverlap of my sig DMPs with published:", length(ov),
    "of", nrow(dmp_sig), "(", round(100*length(ov)/nrow(dmp_sig),1), "% )\n")
ov_strong <- intersect(dmp_sig_strong$cpg, pub_probes)
cat("Overlap of strong DMPs:", length(ov_strong), "of", nrow(dmp_sig_strong),
    "(", round(100*length(ov_strong)/nrow(dmp_sig_strong),1), "% )\n")

# Fisher test of enrichment
n_total <- nrow(dmp)
pub_in_dmp <- length(intersect(pub_probes, dmp$cpg))
tbl <- matrix(c(length(ov_strong),
                nrow(dmp_sig_strong) - length(ov_strong),
                pub_in_dmp - length(ov_strong),
                n_total - nrow(dmp_sig_strong) - pub_in_dmp + length(ov_strong)),
              nrow=2)
ft <- fisher.test(tbl)
cat("Fisher enrichment of published DMPs among strong DMPs: p =",
    format.pval(ft$p.value), "OR =", round(ft$estimate,1), "\n")

# ============================================================
# 3. Annotate DMPs to genes (EPIC manifest, pre-extracted by Python)
# ============================================================
man_file <- file.path(outdir, "epic_manifest_genes.csv")
stopifnot(file.exists(man_file))
cat("\nLoading EPIC manifest gene map...\n")
man <- read.csv(man_file, stringsAsFactors=FALSE)
cat("Manifest rows:", nrow(man), "\n")

# Map CpG -> gene
man_map <- man$gene
names(man_map) <- man$IlmnID
dmp_sig$gene <- man_map[match(dmp_sig$cpg, names(man_map))]
dmp_sig$gene <- sub(";.*", "", dmp_sig$gene)  # first gene only
cat("Annotated DMPs:", sum(!is.na(dmp_sig$gene) & dmp_sig$gene != ""), "/", nrow(dmp_sig), "\n")

# Genes with multiple DMPs
gene_counts <- table(dmp_sig$gene[dmp_sig$gene != "" & !is.na(dmp_sig$gene)])
cat("Unique genes with DMPs:", length(gene_counts), "\n")
cat("Top DMP genes:\n")
print(head(sort(gene_counts, decreasing=TRUE), 20))

# ============================================================
# 4. Check 6-gene signature methylation status
# ============================================================
cat("\n--- Signature gene DMP status ---\n")
sig6 <- c("BATF", "BTG1", "CXCR4", "GPBP1", "SRSF2", "STAT6")
for (g in sig6) {
  idx <- which(dmp_sig$gene == g)
  if (length(idx) > 0) {
    cat(g, ": DMPs =", length(idx),
        "| delta beta:", paste(round(dmp_sig$delta_beta[idx], 3), collapse=","), "\n")
  } else {
    cat(g, ": no significant DMP\n")
  }
}

# ============================================================
# 5. Save annotated results
# ============================================================
write.csv(dmp_sig, file.path(outdir, "dmp_significant_annotated.csv"), row.names=FALSE)
cat("\nAnnotated DMPs saved\n")

cat("\n== DMP annotation done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
