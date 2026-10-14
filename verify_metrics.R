# Verify the reviewer-inserted classifier metrics: 5x10-fold nested CV on
# top-100 DMPs (in-fold vectorized Welch ranking), RF; pooled OOF AUC across
# reps + one representative rep's AUC and confusion matrix at 0.5 (B = positive).
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "verify_metrics.log"), append=FALSE)
cat("== verify metrics", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
suppressPackageStartupMessages({library(data.table); library(randomForest); library(pROC)})

beta <- fread(file.path(outdir, "beta_72samples.csv"))
cpg_ids <- beta[[1]]; B <- as.matrix(beta[, -1, with=FALSE]); rownames(B) <- cpg_ids
rm(beta); gc()
labels <- read.csv(file.path(outdir, "GSE242994_rna_labels.csv"), stringsAsFactors=FALSE)
smap <- read.csv(file.path(outdir, "SPEPN_sentrix_map.csv"), stringsAsFactors=FALSE)
labels$paper_id_clean <- sub("^\\s+", "", labels$paper_id)
smap$spid_clean <- sub("^\\s+", "", smap[[1]])
mapped <- merge(labels, smap, by.x="paper_id_clean", by.y="spid_clean", all.x=TRUE)
subtype <- mapped$Transcriptional_subtype[match(colnames(B), mapped$Sentrix_ID)]
y <- factor(subtype, levels=c("A","B"))
cat("samples: A", sum(y=="A"), "B", sum(y=="B"), "\n")
X <- t(B)  # samples x probes
rm(B); gc()

welch_p <- function(tr) {
  ya <- y[tr] == "A"; yb <- y[tr] == "B"
  ma <- colMeans(X[tr, , drop=FALSE][ya, , drop=FALSE])
  mb <- colMeans(X[tr, , drop=FALSE][yb, , drop=FALSE])
  va <- apply(X[tr, , drop=FALSE][ya, , drop=FALSE], 2, var)
  vb <- apply(X[tr, , drop=FALSE][yb, , drop=FALSE], 2, var)
  na <- sum(ya); nb <- sum(yb)
  se <- sqrt(va/na + vb/nb)
  tstat <- (ma - mb) / se
  df <- (va/na + vb/nb)^2 / ((va/na)^2/(na-1) + (vb/nb)^2/(nb-1))
  p <- 2 * pt(-abs(tstat), df)
  p[!is.finite(p)] <- 1
  p
}

n_reps <- 5
all_preds <- list()
for (rep in 1:n_reps) {
  set.seed(rep * 1000 + 42)
  folds <- sample(rep(1:10, length.out=nrow(X)))
  preds <- rep(NA_real_, nrow(X))
  for (f in 1:10) {
    te <- which(folds == f); tr <- which(folds != f)
    p <- welch_p(tr)
    top <- names(sort(p))[1:100]
    m <- randomForest(X[tr, top, drop=FALSE], y[tr], ntree=1000)
    preds[te] <- predict(m, X[te, top, drop=FALSE], type="prob")[, "B"]
  }
  all_preds[[rep]] <- preds
  cat("rep", rep, "done\n")
}
pooled <- do.call(rbind, lapply(all_preds, function(p) data.frame(y=as.character(y), p=p)))
ro_pooled <- roc(pooled$y, pooled$p, levels=c("A","B"))
cat("pooled (5 reps) OOF AUC:", round(auc(ro_pooled), 3),
    "| DeLong CI:", paste(round(as.numeric(ci.auc(ro_pooled)), 3), collapse="-"), "\n")
rep1 <- all_preds[[1]]
ro1 <- roc(as.character(y), rep1, levels=c("A","B"))
cat("rep-1 OOF AUC:", round(auc(ro1), 3), "\n")
tab <- table(true=as.character(y), pred=ifelse(rep1 >= 0.5, "B", "A"))
print(tab)
ta <- sum(y=="A"); tb <- sum(y=="B")
TP <- sum(y=="B" & rep1 >= 0.5); FN <- tb - TP
FP <- sum(y=="A" & rep1 >= 0.5); TN <- ta - FP
cat(sprintf("sens (B) %.3f | spec %.3f | PPV %.3f | NPV %.3f | acc %.3f\n",
    TP/tb, TN/ta, TP/(TP+FP), TN/(TN+FN), (TP+TN)/(ta+tb)))
cat("== done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
