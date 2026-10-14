# GSE243240: SPINE subtype structure validation — cluster EPN_SPINE and correlate
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
sink(file.path(outdir, "gse243240_spine_subtypes.log"), append=FALSE)
cat("== SPINE subtype structure", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# Load data
beta <- readRDS(file.path(outdir, "gse243240_beta_dmp.rds"))
smap <- read.csv(file.path(outdir, "gse243240_287_map.csv"), stringsAsFactors=FALSE)
dmp_train <- read.csv(file.path(outdir, "dmp_results.csv"), stringsAsFactors=FALSE)

col_gsm <- sub("^(GSM[0-9]+)_.*", "\\1", colnames(beta))
rownames(smap) <- smap$gsm

# SPINE only
sp_gsm <- smap$gsm[smap$class_v12 == "EPN_SPINE"]
sp_beta <- beta[, col_gsm %in% sp_gsm]
colnames(sp_beta) <- col_gsm[col_gsm %in% sp_gsm]
sp <- smap[colnames(sp_beta), ]
cat("EPN_SPINE with beta:", nrow(sp), "\n")

# ============================================================
# 1. Unsupervised clustering of 94 SPINE on top 1000 DMPs
# ============================================================
set.seed(42)
sp_beta_z <- t(scale(t(sp_beta)))
dist_sp <- dist(t(sp_beta_z))
hc_sp <- hclust(dist_sp, method="ward.D2")
k2 <- cutree(hc_sp, k=2)
k3 <- cutree(hc_sp, k=3)
cat("\nK=2 clusters:\n"); print(table(k2))
cat("K=3 clusters:\n"); print(table(k3))

# Compare k=2 with training direction: compute mean delta-beta signature per cluster
db_dir <- dmp_train$delta_beta[match(rownames(sp_beta), dmp_train$cpg)]
score_sp <- apply(sp_beta, 2, function(x) cor(x, db_dir, use="complete.obs"))
cat("\nScore by k2 cluster:\n")
print(tapply(score_sp, k2, mean))

# Which clustering matches A/B direction better?
# A/B direction: subtype A had higher beta at positive-delta probes
# (delta_beta = A - B). So cluster with positive score = A-like
cat("\nA-like (positive score) proportion:\n")
for (cl in 1:2) cat("  cluster", cl, ":", round(mean(score_sp[k2==cl] > 0)*100, 1), "% positive\n")

# ============================================================
# 2. Clinical correlates of the 2 clusters
# ============================================================
sp$k2 <- k2[match(sp$gsm, names(k2))]
sp$score <- score_sp[sp$gsm]

cat("\n--- Clinical characteristics by cluster ---\n")
for (cl in sort(unique(sp$k2))) {
  idx <- sp$k2 == cl
  ages <- as.numeric(sp$age[idx])
  cat(sprintf("Cluster %d: n=%d, age mean=%.1f (n=%d), female=%d/%d\n",
              cl, sum(idx),
              mean(ages, na.rm=TRUE), sum(!is.na(ages)),
              sum(sp$sex[idx] %in% c("female", "F", "f"), na.rm=TRUE), sum(idx)))
}
cat("\nSex table:\n"); print(table(sp$k2, sp$sex))
cat("\nLocation table:\n"); print(table(sp$k2, sp$loc))

# ============================================================
# 3. Survival by cluster
# ============================================================
surv_ok <- !is.na(sp$pfs_months) & sp$pfs_months != "" & sp$pfs_months != "NA"
cat("\nPFS data available:", sum(surv_ok), "of", nrow(sp), "\n")
if (sum(surv_ok) >= 10) {
  sp$pfs_m <- as.numeric(sp$pfs_months)
  sp$pfs_e <- ifelse(sp$pfs_status == "1", 1, 0)
  sf <- survfit(Surv(pfs_m, pfs_e) ~ k2, data=sp[surv_ok, ])
  pdf(file.path(figdir, "Fig3_pfs_spine_clusters.pdf"), width=6, height=6)
  plot(sf, col=c("blue", "red"), lwd=2, xlab="Months", ylab="PFS probability",
       main="GSE243240 EPN_SPINE: PFS by DMP-defined cluster")
  legend("topright", paste("cluster", 1:2), col=c("blue", "red"), lwd=2)
  dev.off()
  sd_test <- survdiff(Surv(pfs_m, pfs_e) ~ k2, data=sp[surv_ok, ])
  cat("Log-rank p:", round(1 - pchisq(sd_test$chisq, df=1), 4), "\n")
} else {
  cat("Insufficient PFS data\n")
}

# ============================================================
# 4. Save results
# ============================================================
out <- data.frame(gsm=sp$gsm, sample=sp$sample, score=sp$score,
                  cluster_k2=sp$k2, cluster_k3=k3[match(sp$gsm, names(k3))],
                  age=sp$age, sex=sp$sex, loc=sp$loc,
                  pfs_months=sp$pfs_months, pfs_status=sp$pfs_status)
write.csv(out, file.path(outdir, "gse243240_spine_subtypes.csv"), row.names=FALSE)

# density plot of score showing bimodality
pdf(file.path(figdir, "Fig3_score_density.pdf"), width=6, height=5)
ggplot(data.frame(score=score_sp), aes(x=score)) +
  geom_histogram(bins=30, fill="steelblue", alpha=0.7) +
  geom_density(alpha=0.3, color="darkblue") +
  labs(title="EPN_SPINE (n=94): subtype score distribution",
       x="Correlation with training A-B DMP direction", y="Count") +
  theme_minimal(base_size=12)
dev.off()
cat("\nResults saved\n")

cat("\n== SPINE subtype analysis done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
