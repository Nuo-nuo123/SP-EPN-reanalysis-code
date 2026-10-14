# GSE243240 validation: apply DMP panel to 287-sample independent cohort
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
figdir <- "D:/claude-study/Spinal Ependymoma/manuscript/figures"
sink(file.path(outdir, "gse243240_validation.log"), append=FALSE)
cat("== GSE243240 validation start", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# ============================================================
# 1. Load validation beta + sample map + training DMPs
# ============================================================
beta <- readRDS(file.path(outdir, "gse243240_beta_dmp.rds"))
cat("Validation beta:", nrow(beta), "probes x", ncol(beta), "samples\n")

smap <- read.csv(file.path(outdir, "gse243240_287_map.csv"), stringsAsFactors=FALSE)
rownames(smap) <- smap$gsm
# beta column names are full IDAT basenames; extract GSM prefix
col_gsm <- sub("^(GSM[0-9]+)_.*", "\\1", colnames(beta))
common <- intersect(col_gsm, smap$gsm)
cat("Samples matched to metadata:", length(common), "\n")

dmp_train <- read.csv(file.path(outdir, "dmp_results.csv"), stringsAsFactors=FALSE)
rownames(dmp_train) <- dmp_train$cpg
probes_use <- intersect(rownames(beta), dmp_train$cpg)
cat("DMP probes available:", length(probes_use), "\n")
db_dir <- dmp_train[probes_use, "delta_beta"]  # A - B direction from training

keep_col <- col_gsm %in% common
beta2 <- beta[probes_use, keep_col]
colnames(beta2) <- col_gsm[keep_col]
sm <- smap[common, ]
cat("Class distribution:\n"); print(table(sm$class_v12))

# ============================================================
# 2. Subtype score = correlation of sample beta profile with
#    training delta-beta vector (A-hypermethylated positive)
# ============================================================
score <- apply(beta2, 2, function(x) cor(x, db_dir, use="complete.obs"))
sm$score <- score[sm$gsm]

# Distribution by class
cat("\nScore by class:\n")
for (cl in sort(unique(sm$class_v12))) {
  idx <- sm$class_v12 == cl
  if (sum(idx) >= 5) {
    cat(sprintf("  %-20s n=%3d mean=%+.3f sd=%.3f\n",
                cl, sum(idx), mean(sm$score[idx]), sd(sm$score[idx])))
  }
}

# ============================================================
# 3. Focus on spinal classes: does EPN_SPINE split into A/B-like groups?
# ============================================================
cat("\n--- EPN_SPINE subtype structure ---\n")
sp_idx <- sm$class_v12 == "EPN_SPINE"
sp <- sm[sp_idx, ]
cat("EPN_SPINE n:", nrow(sp), "\n")

# bimodality test: k-means 2 on score
set.seed(42)
if (nrow(sp) >= 20) {
  km <- kmeans(sp$score, centers=2, nstart=25)
  sp$km2 <- km$cluster
  cat("K-means 2 on score:\n")
  print(table(sp$km2))
  cat("Cluster means:", round(tapply(sp$score, sp$km2, mean), 3), "\n")
  # gap between clusters
  cat("Cluster separation:", round(diff(sort(tapply(sp$score, sp$km2, mean))), 3), "\n")
} else {
  cat("Too few EPN_SPINE samples for k-means\n")
}

# ============================================================
# 4. Unsupervised hierarchical clustering on full DMP panel
# ============================================================
cat("\n--- Hierarchical clustering (top 500 DMPs by |dB|) ---\n")
top500 <- dmp_train$cpg[order(dmp_train$pvalue)][1:500]
top500 <- intersect(top500, probes_use)
b500 <- beta2[top500, ]
b500 <- t(scale(t(b500)))  # z-score per probe
dist_m <- dist(t(b500))
hc <- hclust(dist_m, method="ward.D2")
k5 <- cutree(hc, k=5)
tab <- table(cluster=k5, class=sm$class_v12)
print(tab)
write.csv(tab, file.path(outdir, "gse243240_cluster_vs_class.csv"))

# ============================================================
# 5. Spinal vs non-spinal discrimination by our DMP panel
# ============================================================
cat("\n--- Spinal vs non-spinal classes ---\n")
sm$is_spinal <- sm$class_v12 %in% c("EPN_SPINE", "EPN_SPINE_MYCN", "EPN_MPE_A", "EPN_MPE_B")
sp_idx2 <- sm$is_spinal
if (sum(sp_idx2) >= 10 && sum(!sp_idx2) >= 10) {
  wt <- wilcox.test(sm$score[sp_idx2], sm$score[!sp_idx2])
  cat(sprintf("Spinal (n=%d) vs others (n=%d): score %.3f vs %.3f, p=%.4g\n",
              sum(sp_idx2), sum(!sp_idx2),
              mean(sm$score[sp_idx2]), mean(sm$score[!sp_idx2]), wt$p.value))
  library(pROC)
  roc_sp <- roc(sp_idx2, sm$score)
  cat("AUC:", round(auc(roc_sp), 3), "\n")
}

# ============================================================
# 6. Survival analysis (limited data)
# ============================================================
cat("\n--- Survival analysis (samples with data) ---\n")
surv_ok <- !is.na(sm$pfs_months) & sm$pfs_months != "" & sm$pfs_months != "NA"
if (sum(surv_ok) >= 20) {
  sm$pfs_m <- as.numeric(sm$pfs_months)
  sm$pfs_e <- ifelse(sm$pfs_status == "1", 1, 0)
  sm$score_hi <- sm$score > median(sm$score, na.rm=TRUE)
  sf <- survfit(Surv(pfs_m, pfs_e) ~ score_hi, data=sm[surv_ok, ])
  pdf(file.path(figdir, "Fig3_pfs_gse243240.pdf"), width=6, height=6)
  plot(sf, col=c("blue", "red"), lwd=2, xlab="Months", ylab="PFS",
       main="GSE243240: PFS by subtype score (high vs low)")
  legend("topright", c("score low", "score high"), col=c("blue", "red"), lwd=2)
  dev.off()
  cx <- coxph(Surv(pfs_m, pfs_e) ~ score, data=sm[surv_ok, ])
  cat("Cox PFS ~ score: HR =", round(exp(coef(cx)), 2), "p =", round(summary(cx)$coef[5], 4), "\n")
  cat("Samples with PFS data:", sum(surv_ok), "events:", sum(sm$pfs_e[surv_ok]), "\n")
} else {
  cat("Insufficient survival data (n=", sum(surv_ok), ")\n")
}

# save scores
write.csv(sm[, c("gsm", "sample", "class_v12", "score", "os_months", "os_status",
                 "pfs_months", "pfs_status")],
          file.path(outdir, "gse243240_subtype_scores.csv"), row.names=FALSE)
cat("\nScores saved\n")

cat("\n== Validation done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
