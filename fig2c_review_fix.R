# Fig2C rebuilt for review: (A) volcano from the CORRECT filtered Welch table
# (18,556 genes, log2FC of log2(CPM+1), A-B); (B) ROC with representative-fold
# labels + mean/SD + non-truncated legend.
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages({library(data.table); library(randomForest); library(pROC); library(ggplot2); library(gridExtra)})
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
figdir <- "D:/claude-study/Spinal Ependymoma/figures"

# ---- volcano from the filtered Welch table ----
deg <- fread(file.path(outdir, "deg_reanalysis_results.csv"))  # 18,556 genes; logFC = A - B
deg$sig <- "Other"
deg$sig[deg$pvalue < 0.01] <- "p < 0.01 (140 genes)"
deg$sig[abs(deg$logFC) > 1] <- "|log2FC| > 1 (14 genes)"
p_vol <- ggplot(deg, aes(x=logFC, y=-log10(pvalue), color=sig)) +
  geom_point(size=0.4, alpha=0.6) +
  scale_color_manual(values=c("Other"="grey60", "p < 0.01 (140 genes)"="darkorange",
                              "|log2FC| > 1 (14 genes)"="#d62728"), name=NULL) +
  labs(title="Bulk DEG: subtype A vs B (18,556 genes; none at FDR < 0.05)",
       x="log2 fold change (A - B)", y=expression(-log[10](p))) +
  theme_minimal(base_size=10) +
  theme(legend.position="bottom", legend.text=element_text(size=7),
        legend.key.size=unit(0.3, "cm")) +
  guides(color=guide_legend(nrow=2))

# ---- ROC: representative folds with explicit labels ----
expr_raw <- readRDS(file.path(outdir, "GSE242994_expr.rds"))
expr <- if (is.list(expr_raw)) expr_raw$expr else expr_raw
labels <- read.csv(file.path(outdir, "GSE242994_rna_labels.csv"), stringsAsFactors=FALSE)
common <- intersect(colnames(expr), labels$sample)
y <- as.factor(labels$Transcriptional_subtype[match(common, labels$sample)])
lib_sizes <- colSums(expr[, common])
cpm <- sweep(expr[, common], 2, lib_sizes / 1e6, "/")
logcpm <- log2(cpm + 1)
expressed <- rowSums(cpm >= 1) >= round(0.2 * length(common))
Xt <- t(logcpm[expressed, ])
set.seed(501)
folds <- sample(rep(1:10, length.out=nrow(Xt)))
preds <- rep(NA, nrow(Xt))
for (f in 1:10) {
  te <- which(folds == f); tr <- which(folds != f)
  gvar <- apply(Xt[tr, , drop=FALSE], 2, var)
  top <- names(sort(gvar, decreasing=TRUE))[1:4000]
  m <- randomForest(Xt[tr, top, drop=FALSE], y[tr], ntree=1000)
  preds[te] <- predict(m, Xt[te, top, drop=FALSE], type="prob")[,2]
}
roc_nest <- roc(as.numeric(y)-1, preds)
X_raw <- t(as.matrix(expr[, common]))
gvar_all <- apply(X_raw, 2, var)
top50 <- names(sort(gvar_all, decreasing=TRUE))[1:50]
X50 <- X_raw[, top50]
set.seed(1001); folds50 <- sample(rep(1:10, length.out=nrow(X50)))
cv50 <- rep(NA, nrow(X50))
for (i in 1:10) {
  tr <- folds50 != i; te <- folds50 == i
  m <- randomForest(X50[tr,,drop=FALSE], y[tr], ntree=1000)
  cv50[te] <- predict(m, X50[te,,drop=FALSE], type="prob")[,2]
}
roc50 <- roc(as.numeric(y)-1, cv50)
sig6 <- c("BATF", "BTG1", "CXCR4", "GPBP1", "SRSF2", "STAT6")
X6 <- X_raw[, sig6]
set.seed(2001); folds6 <- sample(rep(1:10, length.out=nrow(X6)))
cv6 <- rep(NA, nrow(X6))
for (i in 1:10) {
  tr <- folds6 != i; te <- folds6 == i
  m <- randomForest(X6[tr,,drop=FALSE], y[tr], ntree=1000)
  cv6[te] <- predict(m, X6[te,,drop=FALSE], type="prob")[,2]
}
roc6 <- roc(as.numeric(y)-1, cv6)
cat("representative AUCs: nested", round(auc(roc_nest), 3),
    "| 50-gene", round(auc(roc50), 3), "| 6-gene", round(auc(roc6), 3), "\n")

p_roc_data <- data.frame(
  fpr=c(1-roc_nest$specificities, 1-roc50$specificities, 1-roc6$specificities),
  tpr=c(roc_nest$sensitivities, roc50$sensitivities, roc6$sensitivities),
  grp=factor(c(
    rep(sprintf("Nested CV top-4,000: representative fold %.3f (mean 0.530 +/- 0.027)", auc(roc_nest)),
        length(roc_nest$sensitivities)),
    rep(sprintf("50-gene ensemble, conventional CV: representative fold %.3f (mean 0.519 +/- 0.037)", auc(roc50)),
        length(roc50$sensitivities)),
    rep(sprintf("6-gene signature, conventional CV: representative fold %.3f (mean 0.537 +/- 0.036)", auc(roc6)),
        length(roc6$sensitivities)))))
p_roc <- ggplot(p_roc_data, aes(x=fpr, y=tpr, color=grp, linetype=grp)) +
  geom_line(linewidth=0.9) +
  geom_abline(slope=1, intercept=0, linetype=2, color="grey60") +
  labs(title="Transcriptomic classifiers (representative folds; means per design in legend)",
       x="False positive rate", y="True positive rate", color=NULL, linetype=NULL) +
  scale_color_manual(values=c("forestgreen", "darkorange", "steelblue")) +
  scale_linetype_manual(values=c(1, 2, 3)) +
  coord_fixed() +
  theme_minimal(base_size=9) +
  theme(legend.position="bottom", legend.text=element_text(size=6.5),
        legend.key.width=unit(1.4, "cm"), plot.margin=margin(5, 8, 5, 5)) +
  guides(color=guide_legend(ncol=1), linetype=guide_legend(ncol=1))

pdf(file.path(figdir, "Fig2C_volcano.pdf"), width=6.5, height=6)
print(p_vol)
dev.off()
pdf(file.path(figdir, "Fig2D_roc.pdf"), width=6.5, height=6)
print(p_roc)
dev.off()
cat("Fig2C volcano + Fig2D ROC saved\n")
