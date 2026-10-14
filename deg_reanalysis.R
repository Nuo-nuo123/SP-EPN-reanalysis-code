# DEG re-analysis: proper pipeline with low-expression filtering
# 1. Load raw counts (GSE242992) or existing expr
# 2. Filter low-expression genes (CPM>=1 in >=20% samples)
# 3. limma-voom style: log2(CPM) + t-test with BH
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(data.table)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "deg_reanalysis.log"), append=FALSE)
cat("== DEG re-analysis start", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# Load expression (count matrix from GSE242992 raw)
expr_raw <- readRDS(file.path(outdir, "GSE242994_expr.rds"))
expr <- if (is.list(expr_raw)) expr_raw$expr else expr_raw
cat("Expr dim:", nrow(expr), "x", ncol(expr), "\n")
cat("Value range:", range(expr, na.rm=TRUE), "\n")
cat("Mean per gene sample:", head(sort(rowMeans(expr), decreasing=TRUE), 3), "\n")

labels <- read.csv(file.path(outdir, "GSE242994_rna_labels.csv"), stringsAsFactors=FALSE)
common <- intersect(colnames(expr), labels$sample)
y <- labels$Transcriptional_subtype[match(common, labels$sample)]
cat("Samples:", length(common), "A:", sum(y=="A"), "B:", sum(y=="B"), "\n")

# ============================================================
# 1. Low-expression filtering: CPM >= 1 in at least 20% of samples
# ============================================================
lib_sizes <- colSums(expr[, common])
cpm <- sweep(expr[, common], 2, lib_sizes / 1e6, "/")
expressed <- rowSums(cpm >= 1) >= round(0.2 * length(common))
cat("Genes passing CPM>=1 in >=20% samples:", sum(expressed), "\n")
expr_f <- expr[expressed, common]
cat("Filtered:", nrow(expr_f), "genes\n")

# ============================================================
# 2. log2 CPM + Welch t-test + BH
# ============================================================
logcpm <- log2(cpm[expressed, , drop=FALSE] + 1)
grpA <- which(y == "A")
grpB <- which(y == "B")

cat("Running tests...\n")
n <- nrow(logcpm)
pvals <- numeric(n)
lfc <- numeric(n)
for (i in 1:n) {
  a <- as.numeric(logcpm[i, grpA]); b <- as.numeric(logcpm[i, grpB])
  tt <- tryCatch(t.test(a, b), error=function(e) NULL)
  if (is.null(tt)) { pvals[i] <- NA; lfc[i] <- NA; next }
  pvals[i] <- tt$p.value
  lfc[i] <- mean(a) - mean(b)
}
padj <- p.adjust(pvals, method="BH")
cat("Significant (padj<0.05):", sum(padj < 0.05, na.rm=TRUE), "\n")
cat("Significant (padj<0.05 & |logFC|>1):", sum(padj < 0.05 & abs(lfc) > 1, na.rm=TRUE), "\n")
cat("Min p:", min(pvals, na.rm=TRUE), "gene:", rownames(logcpm)[which.min(pvals)], "\n")
cat("p<0.01 count:", sum(pvals < 0.01, na.rm=TRUE), "\n")
cat("p<0.001 count:", sum(pvals < 0.001, na.rm=TRUE), "\n")
cat("|logFC|>1 count:", sum(abs(lfc) > 1, na.rm=TRUE), "\n")
cat("|logFC|>2 count:", sum(abs(lfc) > 2, na.rm=TRUE), "\n")

# ============================================================
# 3. Save results
# ============================================================
res <- data.frame(gene=rownames(logcpm), logFC=lfc, pvalue=pvals, padj=padj)
res <- res[order(res$pvalue), ]
write.csv(res, file.path(outdir, "deg_reanalysis_results.csv"), row.names=FALSE)

# Top 20 by p
cat("\nTop 20 by p:\n")
print(head(res[!is.na(res$pvalue), ], 20))

# ============================================================
# 4. Sensitivity: Molecular_subtype labels instead of Transcriptional
# ============================================================
cat("\n--- Sensitivity: Molecular_subtype labels ---\n")
y_mol <- labels$Molecular_subtype[match(common, labels$sample)]
cat("Molecular labels: A:", sum(y_mol=="A"), "B:", sum(y_mol=="B"),
    "NA:", sum(is.na(y_mol) | y_mol==""), "\n")
# 6 discordant: Transcriptional=A, Molecular=B
disc <- labels$sample[labels$Transcriptional_subtype=="A" & labels$Molecular_subtype=="B"]
cat("Discordant samples:", length(disc), "->", paste(disc, collapse=","), "\n")

# redo DEG with molecular labels
keep_mol <- y_mol %in% c("A", "B")
logcpm2 <- logcpm[, keep_mol]
y2 <- y_mol[keep_mol]
gA <- which(y2 == "A"); gB <- which(y2 == "B")
pvals2 <- numeric(nrow(logcpm2)); lfc2 <- numeric(nrow(logcpm2))
for (i in 1:nrow(logcpm2)) {
  a <- as.numeric(logcpm2[i, gA]); b <- as.numeric(logcpm2[i, gB])
  tt <- tryCatch(t.test(a, b), error=function(e) NULL)
  if (is.null(tt)) { pvals2[i] <- NA; lfc2[i] <- NA; next }
  pvals2[i] <- tt$p.value
  lfc2[i] <- mean(a) - mean(b)
}
padj2 <- p.adjust(pvals2, method="BH")
cat("Molecular-label DEG: sig (padj<0.05):", sum(padj2 < 0.05, na.rm=TRUE), "\n")
cat("Molecular-label DEG: min p:", min(pvals2, na.rm=TRUE), "\n")

res2 <- data.frame(gene=rownames(logcpm2), logFC=lfc2, pvalue=pvals2, padj=padj2)
res2 <- res2[order(res2$pvalue), ]
write.csv(res2, file.path(outdir, "deg_sensitivity_molecular_labels.csv"), row.names=FALSE)

# ============================================================
# 5. Multivariate check: RF on top variable genes (from logcpm)
# ============================================================
cat("\n--- Multivariate signal check ---\n")
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
# top 1000 most variable genes
gvar <- apply(logcpm, 1, var)
top_var <- names(sort(gvar, decreasing=TRUE))[1:1000]
Xv <- t(logcpm[top_var, ])
set.seed(42)
folds <- sample(rep(1:10, length.out=nrow(Xv)))
cv_pred <- rep(NA, nrow(Xv))
for (i in 1:10) {
  tr <- folds != i; te <- folds == i
  m <- randomForest(Xv[tr,,drop=FALSE], as.factor(y[tr]), ntree=1000)
  cv_pred[te] <- predict(m, Xv[te,,drop=FALSE], type="prob")[,2]
}
roc_v <- roc(as.numeric(as.factor(y))-1, cv_pred)
cat("RF CV AUC (top 1000 variable genes):", round(auc(roc_v), 3), "\n")
# top 4000 variable (matching the source study's 4000)
top4k <- names(sort(gvar, decreasing=TRUE))[1:min(4000, length(gvar))]
X4 <- t(logcpm[top4k, ])
set.seed(42)
folds4 <- sample(rep(1:10, length.out=nrow(X4)))
cv4 <- rep(NA, nrow(X4))
for (i in 1:10) {
  tr <- folds4 != i; te <- folds4 == i
  m <- randomForest(X4[tr,,drop=FALSE], as.factor(y[tr]), ntree=1000)
  cv4[te] <- predict(m, X4[te,,drop=FALSE], type="prob")[,2]
}
roc4 <- roc(as.numeric(as.factor(y))-1, cv4)
cat("RF CV AUC (top 4000 variable genes):", round(auc(roc4), 3), "\n")

saveRDS(list(deg=res, deg_mol=res2, roc_1000=roc_v, roc_4000=roc4),
        file.path(outdir, "deg_reanalysis.rds"))
cat("\n== DEG re-analysis done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
