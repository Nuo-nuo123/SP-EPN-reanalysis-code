# GSE243240 CONSISTENT re-analysis: one documented pipeline
# Primary: hierarchical clustering (Ward.D2, Euclidean, z-scored beta) of 94 EPN_SPINE
#          on all 1,000 discovery DMPs, cut at k=2.
# Sensitivity: k-means (k=2, nstart=25) on the subtype score.
# Age comparison with Wilcoxon test.
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
figdir <- "D:/claude-study/Spinal Ependymoma/manuscript/figures"
sink(file.path(outdir, "gse243240_consistent.log"), append=FALSE)
cat("== GSE243240 consistent pipeline", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

beta <- readRDS(file.path(outdir, "gse243240_beta_dmp.rds"))  # 1000 probes x 287 samples
smap <- read.csv(file.path(outdir, "gse243240_287_map.csv"), stringsAsFactors=FALSE)
dmp_train <- read.csv(file.path(outdir, "dmp_results.csv"), stringsAsFactors=FALSE)
col_gsm <- sub("^(GSM[0-9]+)_.*", "\\1", colnames(beta))
rownames(smap) <- smap$gsm

sp_gsm <- smap$gsm[smap$class_v12 == "EPN_SPINE"]
sp_beta <- beta[, col_gsm %in% sp_gsm]
colnames(sp_beta) <- col_gsm[col_gsm %in% sp_gsm]
sp <- smap[colnames(sp_beta), ]
cat("EPN_SPINE:", nrow(sp), "probes:", nrow(sp_beta), "\n")

# ============================================================
# PRIMARY: hierarchical clustering on z-scored beta (1000 DMPs), k=2
# ============================================================
set.seed(42)
sp_z <- t(scale(t(sp_beta)))
d <- dist(t(sp_z))
hc <- hclust(d, method="ward.D2")
k2 <- cutree(hc, k=2)
cat("\nPRIMARY (hierarchical, 1000 DMPs, k=2):\n")
print(table(k2))

# subtype score per sample (A-B direction from training)
db_dir <- dmp_train$delta_beta[match(rownames(sp_beta), dmp_train$cpg)]
score_sp <- apply(sp_beta, 2, function(x) cor(x, db_dir, use="complete.obs"))
cat("Score by cluster:\n")
print(tapply(score_sp, k2, mean))
cat("Positive-score proportion by cluster:\n")
for (cl in sort(unique(k2))) cat("  cluster", cl, ":", round(mean(score_sp[k2==cl] > 0)*100, 1), "%\n")

# ============================================================
# SENSITIVITY: k-means on score (k=2, nstart=25)
# ============================================================
set.seed(42)
km <- kmeans(score_sp, centers=2, nstart=25)
cat("\nSENSITIVITY (k-means on score, k=2, nstart=25):\n")
print(table(km$cluster))
cat("Cluster means:", round(tapply(score_sp, km$cluster, mean), 3), "\n")
# agreement between two methods
tryCatch({
  suppressPackageStartupMessages(library(mclust))
  cat("ARI (hierarchical vs k-means):", round(adjustedRandIndex(k2, km$cluster), 3), "\n")
}, error=function(e) cat("mclust unavailable for ARI\n"))

# ============================================================
# Age comparison (Wilcoxon) — review item 17
# ============================================================
sp$k2 <- k2[match(sp$gsm, names(k2))]
sp$score <- score_sp[sp$gsm]
sp$age_n <- suppressWarnings(as.numeric(sp$age))
for (cl in sort(unique(sp$k2))) {
  idx <- sp$k2 == cl
  ages <- sp$age_n[idx]
  cat(sprintf("\nCluster %d (n=%d): age n=%d mean=%.1f sd=%.1f\n",
              cl, sum(idx), sum(!is.na(ages)), mean(ages, na.rm=TRUE), sd(ages, na.rm=TRUE)))
}
wt_age <- wilcox.test(sp$age_n[sp$k2==1], sp$age_n[sp$k2==2])
cat("Age Wilcoxon (cluster 1 vs 2): p =", format.pval(wt_age$p.value), "\n")

# ============================================================
# Full-cohort survival (175 samples with PFS) — report as negative
# ============================================================
# recompute score for all 287 samples
beta_all <- beta[, col_gsm %in% smap$gsm]
colnames(beta_all) <- col_gsm[col_gsm %in% smap$gsm]
score_all <- apply(beta_all, 2, function(x) cor(x, db_dir[match(rownames(beta_all), dmp_train$cpg)], use="complete.obs"))
sm2 <- smap[names(score_all), ]
sm2$score <- score_all
surv_ok <- !is.na(sm2$pfs_months) & sm2$pfs_months != "" & sm2$pfs_months != "NA"
sm2$pfs_m <- suppressWarnings(as.numeric(sm2$pfs_months))
sm2$pfs_e <- ifelse(sm2$pfs_status == "1", 1, 0)
ok2 <- surv_ok & !is.na(sm2$pfs_m) & sm2$pfs_m > 0
cat("\nFull cohort with PFS:", sum(ok2), "events:", sum(sm2$pfs_e[ok2]), "\n")
if (sum(ok2) >= 20) {
  cx <- coxph(Surv(pfs_m, pfs_e) ~ score, data=sm2[ok2, ])
  cat("Cox PFS ~ score (all classes): HR =", round(exp(coef(cx)), 3),
      "p =", round(summary(cx)$coef[5], 4), "\n")
  # stratified by class
  cx2 <- tryCatch(coxph(Surv(pfs_m, pfs_e) ~ score + strata(class_v12), data=sm2[ok2, ]),
                  error=function(e) NULL)
  if (!is.null(cx2)) cat("Cox PFS ~ score + strata(class): HR =", round(exp(coef(cx2)), 3),
                         "p =", round(summary(cx2)$coef[5], 4), "\n")
  # within SPINE only
  sp_ok <- ok2 & sm2$class_v12 == "EPN_SPINE"
  cat("EPN_SPINE with PFS:", sum(sp_ok), "\n")
  if (sum(sp_ok) >= 10) {
    cx3 <- coxph(Surv(pfs_m, pfs_e) ~ score, data=sm2[sp_ok, ])
    cat("Cox PFS ~ score (SPINE only): HR =", round(exp(coef(cx3)), 3),
        "p =", round(summary(cx3)$coef[5], 4), "\n")
  }
}

# save consistent outputs
out <- data.frame(gsm=sp$gsm, sample=sp$sample, score=score_sp,
                  cluster_hier=k2, cluster_kmeans=km$cluster,
                  age=sp$age_n)
write.csv(out, file.path(outdir, "gse243240_spine_consistent.csv"), row.names=FALSE)
cat("\nSaved consistent results\n")

cat("\n== Consistent pipeline done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
