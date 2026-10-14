# Methylation classifier with NESTED cross-validation (no feature-selection leakage)
# Outer: 10-fold CV; Inner: select top-k DMPs within training fold only
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(data.table)
  library(randomForest)
  library(pROC)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "methylation_nested_cv.log"), append=FALSE)
cat("== Methylation nested CV start", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# Load 72-sample beta
beta_sub <- fread(file.path(outdir, "beta_72samples.csv"))
cpg_ids <- beta_sub[[1]]
beta_mat <- as.matrix(beta_sub[, -1, with=FALSE])
rownames(beta_mat) <- cpg_ids
rm(beta_sub); gc()

labels <- read.csv(file.path(outdir, "GSE242994_rna_labels.csv"), stringsAsFactors=FALSE)
smap <- read.csv(file.path(outdir, "SPEPN_sentrix_map.csv"), stringsAsFactors=FALSE)
labels$paper_id_clean <- sub("^\\s+", "", labels$paper_id)
smap$spid_clean <- sub("^\\s+", "", smap[[1]])
mapped <- merge(labels, smap, by.x="paper_id_clean", by.y="spid_clean", all.x=TRUE)
subtype <- mapped$Transcriptional_subtype[match(colnames(beta_mat), mapped$Sentrix_ID)]
cat("Samples:", ncol(beta_mat), "A:", sum(subtype=="A"), "B:", sum(subtype=="B"), "\n")

ym <- as.factor(subtype)

# ============================================================
# Nested CV: outer folds; feature selection INSIDE training folds
# ============================================================
nested_cv <- function(k_top, n_folds=10, n_rep=5) {
  aucs <- c()
  for (rep in 1:n_rep) {
    set.seed(100 + rep)
    folds <- sample(rep(1:n_folds, length.out=nrow(t(beta_mat))))
    preds <- rep(NA, ncol(beta_mat))
    for (f in 1:n_folds) {
      te <- which(folds == f); tr <- which(folds != f)
      # feature selection within training fold only
      pv <- apply(beta_mat[, tr], 1, function(x) {
        a <- x[subtype[tr]=="A"]; b <- x[subtype[tr]=="B"]
        if (sum(!is.na(a))<3 || sum(!is.na(b))<3) return(NA)
        tryCatch(t.test(a, b)$p.value, error=function(e) NA)
      })
      top <- names(sort(pv, na.last=TRUE))[1:min(k_top, sum(!is.na(pv)))]
      Xtr <- t(beta_mat[top, tr, drop=FALSE]); Xtr[is.na(Xtr)] <- 0
      Xte <- t(beta_mat[top, te, drop=FALSE]); Xte[is.na(Xte)] <- 0
      m <- randomForest(Xtr, ym[tr], ntree=1000)
      preds[te] <- predict(m, Xte, type="prob")[,2]
    }
    aucs <- c(aucs, auc(roc(as.numeric(ym)-1, preds)))
  }
  return(aucs)
}

cat("\nRunning nested CV (this takes a few minutes per fold)...\n")
for (k in c(100, 1000)) {
  cat("k =", k, "...\n")
  aucs <- nested_cv(k)
  cat(sprintf("Nested CV AUC (top %d, 5x10-fold): mean=%.3f sd=%.3f range=[%.3f, %.3f]\n",
              k, mean(aucs), sd(aucs), min(aucs), max(aucs)))
}

# ============================================================
# Comparison: naive CV (leaky) as reference
# ============================================================
cat("\nNaive (leaky) CV for comparison...\n")
pv_all <- apply(beta_mat, 1, function(x) {
  a <- x[subtype=="A"]; b <- x[subtype=="B"]
  tryCatch(t.test(a, b)$p.value, error=function(e) NA)
})
for (k in c(100, 1000)) {
  top <- names(sort(pv_all, na.last=TRUE))[1:k]
  Xm <- t(beta_mat[top, , drop=FALSE]); Xm[is.na(Xm)] <- 0
  aucs <- c()
  for (rep in 1:5) {
    set.seed(200 + rep)
    folds <- sample(rep(1:10, length.out=nrow(Xm)))
    preds <- rep(NA, nrow(Xm))
    for (f in 1:10) {
      tr <- folds != f; te <- folds == f
      m <- randomForest(Xm[tr,,drop=FALSE], ym[tr], ntree=1000)
      preds[te] <- predict(m, Xm[te,,drop=FALSE], type="prob")[,2]
    }
    aucs <- c(aucs, auc(roc(as.numeric(ym)-1, preds)))
  }
  cat(sprintf("Naive CV AUC (top %d): mean=%.3f sd=%.3f\n", k, mean(aucs), sd(aucs)))
}

cat("\n== Nested CV done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
