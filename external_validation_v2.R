# External validation v2: 6-gene consensus signature
# Train: GSE242994 (72 spinal ependymoma, A=43/B=29)
# Validate: GSE66787 (35 spinal), GSE64415 (209, molecular subgroups),
#           GSE125861 (69, primary/recurrent), GSE66354 (149, multi-tumor)
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(randomForest)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "external_validation_v2.log"), append=FALSE)
cat("== External validation v2 start", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# ============================================================
# Load training data
# ============================================================
expr_raw <- readRDS(file.path(outdir, "GSE242994_expr.rds"))
if (is.list(expr_raw) && "expr" %in% names(expr_raw)) { expr <- expr_raw$expr } else { expr <- expr_raw }
labels <- fread(file.path(outdir, "GSE242994_rna_labels.csv"))
common <- intersect(colnames(expr), labels$sample)
y <- as.factor(labels$Transcriptional_subtype[match(common, labels$sample)])
X <- t(as.matrix(expr[, common]))
cat("Training:", nrow(X), "samples, A:", sum(y=="A"), "B:", sum(y=="B"), "\n")

sig6 <- c("BATF", "BTG1", "CXCR4", "GPBP1", "SRSF2", "STAT6")
X6 <- X[, intersect(sig6, colnames(X))]
cat("6-gene signature training features:", ncol(X6), "\n")

# ============================================================
# Train models on 6 genes: RF + LR
# ============================================================
set.seed(42)
rf_fit <- randomForest(X6, y, ntree=2000, importance=TRUE)
cat("\nRF OOB error:", rf_fit$err.rate[nrow(rf_fit$err.rate), "OOB"], "\n")

# Internal CV for honest performance estimate
set.seed(42)
folds <- sample(rep(1:10, length.out=nrow(X6)))
cv_pred <- rep(NA, nrow(X6))
for (i in 1:10) {
  tr <- folds != i; te <- folds == i
  m <- randomForest(X6[tr,,drop=FALSE], y[tr], ntree=1000)
  cv_pred[te] <- predict(m, X6[te,,drop=FALSE], type="prob")[,2]
}
roc_cv <- roc(as.numeric(y)-1, cv_pred)
cat("6-gene CV AUC:", round(auc(roc_cv), 3), "\n")

# ============================================================
# Helper: apply signature to external cohort
# ============================================================
apply_signature <- function(cohort_name, expr_file, sample_file, group_col=NULL) {
  cat("\n==========", cohort_name, "==========\n")
  ce_raw <- readRDS(file.path(outdir, expr_file))
  if (is.list(ce_raw) && "expr" %in% names(ce_raw)) { ce <- ce_raw$expr } else { ce <- ce_raw }
  cat("Matrix:", nrow(ce), "genes x", ncol(ce), "samples\n")

  avail <- intersect(sig6, rownames(ce))
  cat("Signature genes available:", length(avail), "->", paste(avail, collapse=","), "\n")
  if (length(avail) < 4) { cat("SKIPPED: too few genes\n"); return(NULL) }

  # Build full 6-column matrix; pad missing genes with training median
  Xe_full <- matrix(NA, nrow=ncol(ce), ncol=length(sig6),
                    dimnames=list(colnames(ce), sig6))
  Xe_avail <- t(as.matrix(ce[avail, , drop=FALSE]))
  Xe_full[, avail] <- Xe_avail
  missing_genes <- setdiff(sig6, avail)
  for (mg in missing_genes) {
    Xe_full[, mg] <- median(X6[, mg])
    cat("  Padding missing gene", mg, "with training median\n")
  }
  pred <- predict(rf_fit, Xe_full, type="prob")[,2]

  # Load sample metadata for group labels
  meta <- fread(file.path(outdir, sample_file))
  meta$pred_B <- pred[match(meta$sample, colnames(ce))]

  # Extract group label from characteristics
  meta$group <- sub(".*?: ", "", meta$characteristics_ch1)

  cat("Prediction summary: mean", round(mean(pred),3), "sd", round(sd(pred),3), "\n")
  cat("Group distribution:\n")
  for (g in sort(unique(meta$group))) {
    idx <- meta$group == g
    if (sum(idx) >= 3) {
      cat(sprintf("  %-45s n=%3d  mean_pred=%5.3f\n", g, sum(idx), mean(meta$pred_B[idx], na.rm=TRUE)))
    }
  }

  # Statistical test if meaningful groups
  groups <- names(table(meta$group))
  if (length(groups) >= 2) {
    # use the two largest groups or all
    gtab <- sort(table(meta$group), decreasing=TRUE)
    if (length(gtab) >= 2) {
      g1 <- names(gtab)[1]; g2 <- names(gtab)[2]
      wt <- wilcox.test(meta$pred_B[meta$group==g1], meta$pred_B[meta$group==g2])
      cat(sprintf("Wilcoxon %s vs %s: p=%.4g\n", g1, g2, wt$p.value))
    }
  }
  return(meta)
}

res66787 <- apply_signature("GSE66787 spinal ependymoma",
                            "GSE66787_expr.rds", "GSE66787_samples.csv")
res64415 <- apply_signature("GSE64415 molecular subgroups",
                            "GSE64415_expr.rds", "GSE64415_samples.csv")
res125861 <- apply_signature("GSE125861 primary/recurrent",
                             "GSE125861_expr.rds", "GSE125861_samples.csv")
res66354 <- apply_signature("GSE66354 multi-tumor",
                            "GSE66354_expr.rds", "GSE66354_samples.csv")

# ============================================================
# GSE64415 focused: spinal vs non-spinal subgroups
# ============================================================
cat("\n=== GSE64415 focused: SP subgroups vs others ===\n")
if (!is.null(res64415)) {
  res64415$is_spinal <- grepl("SP_", res64415$group)
  cat("SP_EPN:", sum(res64415$group=="SP_EPN"), "SP_MPE:", sum(res64415$group=="SP_MPE"), "\n")
  sp_idx <- res64415$is_spinal
  others_idx <- !sp_idx
  if (sum(sp_idx) >= 3 && sum(others_idx) >= 3) {
    wt <- wilcox.test(res64415$pred_B[sp_idx], res64415$pred_B[others_idx])
    cat(sprintf("Spinal subgroups (n=%d) vs others (n=%d): mean %.3f vs %.3f, Wilcoxon p=%.4g\n",
                sum(sp_idx), sum(others_idx),
                mean(res64415$pred_B[sp_idx], na.rm=TRUE),
                mean(res64415$pred_B[others_idx], na.rm=TRUE),
                wt$p.value))
    roc_sp <- roc(sp_idx, res64415$pred_B)
    cat("AUC (spinal vs non-spinal):", round(auc(roc_sp), 3), "\n")
  }
}

# ============================================================
# GSE125861 focused: primary vs recurrent
# ============================================================
cat("\n=== GSE125861: primary vs recurrent ===\n")
if (!is.null(res125861)) {
  pri <- res125861$group == "primary"
  rec <- res125861$group == "recurrent"
  if (sum(pri) >= 3 && sum(rec) >= 3) {
    wt <- wilcox.test(res125861$pred_B[pri], res125861$pred_B[rec])
    cat(sprintf("Primary (n=%d) vs recurrent (n=%d): mean %.3f vs %.3f, Wilcoxon p=%.4g\n",
                sum(pri), sum(rec),
                mean(res125861$pred_B[pri], na.rm=TRUE),
                mean(res125861$pred_B[rec], na.rm=TRUE),
                wt$p.value))
    roc_pr <- roc(rec, res125861$pred_B)
    cat("AUC (recurrent vs primary):", round(auc(roc_pr), 3), "\n")
  }
}

# ============================================================
# GSE66354: ependymoma vs other tumors (specificity)
# ============================================================
cat("\n=== GSE66354: ependymoma vs other brain tumors ===\n")
if (!is.null(res66354)) {
  res66354$is_epen <- grepl("ependymoma", res66354$group)
  epen_idx <- res66354$is_epen
  other_idx <- !epen_idx
  if (sum(epen_idx) >= 3 && sum(other_idx) >= 3) {
    wt <- wilcox.test(res66354$pred_B[epen_idx], res66354$pred_B[other_idx])
    cat(sprintf("Ependymoma (n=%d) vs other tumors (n=%d): mean %.3f vs %.3f, Wilcoxon p=%.4g\n",
                sum(epen_idx), sum(other_idx),
                mean(res66354$pred_B[epen_idx], na.rm=TRUE),
                mean(res66354$pred_B[other_idx], na.rm=TRUE),
                wt$p.value))
  }
}

cat("\n== External validation v2 done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
