# Verify the transcriptome pooled OOF AUC: 3 reps x 10-fold nested CV,
# top-4,000 variance-selected in-fold, RF, pooled predictions + DeLong CI.
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "verify_transcriptome_pooled.log"), append=FALSE)
cat("== verify transcriptome pooled", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
suppressPackageStartupMessages({library(data.table); library(randomForest); library(pROC)})

expr_raw <- readRDS(file.path(outdir, "GSE242994_expr.rds"))
expr <- if (is.list(expr_raw)) expr_raw$expr else expr_raw
labels <- read.csv(file.path(outdir, "GSE242994_rna_labels.csv"), stringsAsFactors=FALSE)
common <- intersect(colnames(expr), labels$sample)
y <- factor(labels$Transcriptional_subtype[match(common, labels$sample)], levels=c("A","B"))
cpm <- sweep(expr[, common], 2, colSums(expr[, common])/1e6, "/")
logcpm <- log2(cpm + 1)
expressed <- rowSums(cpm >= 1) >= round(0.2 * length(common))
X <- t(logcpm[expressed, ])  # 72 x 18,556

n_reps <- 3
all_preds <- list()
for (rep in 1:n_reps) {
  set.seed(rep * 1000 + 501)
  folds <- sample(rep(1:10, length.out=nrow(X)))
  preds <- rep(NA_real_, nrow(X))
  for (f in 1:10) {
    te <- which(folds == f); tr <- which(folds != f)
    gvar <- apply(X[tr, , drop=FALSE], 2, var)
    top <- names(sort(gvar, decreasing=TRUE))[1:4000]
    m <- randomForest(X[tr, top, drop=FALSE], y[tr], ntree=1000)
    preds[te] <- predict(m, X[te, top, drop=FALSE], type="prob")[, "B"]
  }
  all_preds[[rep]] <- preds
  cat("rep", rep, "done\n")
}
pooled <- do.call(rbind, lapply(all_preds, function(p) data.frame(y=as.character(y), p=p)))
ro_pooled <- roc(pooled$y, pooled$p, levels=c("A","B"))
cat("pooled (3 reps) OOF AUC:", round(auc(ro_pooled), 3),
    "| DeLong CI:", paste(round(as.numeric(ci.auc(ro_pooled)), 3), collapse="-"), "\n")
cat("== done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
