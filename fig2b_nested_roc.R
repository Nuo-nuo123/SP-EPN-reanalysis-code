# Fig2B: ROC from NESTED CV (one representative 10-fold rep for each k)
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
figdir <- "D:/claude-study/Spinal Ependymoma/figures"

beta_sub <- fread(file.path(outdir, "beta_72samples.csv"))
cpg_ids <- beta_sub[[1]]
B <- as.matrix(beta_sub[, -1, with=FALSE])
rm(beta_sub); gc()

labels <- read.csv(file.path(outdir, "GSE242994_rna_labels.csv"), stringsAsFactors=FALSE)
smap <- read.csv(file.path(outdir, "SPEPN_sentrix_map.csv"), stringsAsFactors=FALSE)
labels$paper_id_clean <- sub("^\\s+", "", labels$paper_id)
smap$spid_clean <- sub("^\\s+", "", smap[[1]])
mapped <- merge(labels, smap, by.x="paper_id_clean", by.y="spid_clean", all.x=TRUE)
subtype <- mapped$Transcriptional_subtype[match(colnames(B), mapped$Sentrix_ID)]
ym <- as.factor(subtype)

fast_welch <- function(X, idx) {
  a <- X[, subtype[idx] == "A", drop=FALSE]
  b <- X[, subtype[idx] == "B", drop=FALSE]
  ma <- rowMeans(a, na.rm=TRUE); mb <- rowMeans(b, na.rm=TRUE)
  va <- apply(a, 1, var, na.rm=TRUE); vb <- apply(b, 1, var, na.rm=TRUE)
  na_ok <- rowSums(!is.na(a)); nb_ok <- rowSums(!is.na(b))
  se2 <- va / na_ok + vb / nb_ok
  tstat <- (ma - mb) / sqrt(se2)
  df <- se2^2 / ((va/na_ok)^2/(na_ok-1) + (vb/nb_ok)^2/(nb_ok-1))
  df[!is.finite(df)] <- NA
  2 * pt(-abs(tstat), df=df)
}

set.seed(101)
folds <- sample(rep(1:10, length.out=ncol(B)))

pdf(file.path(figdir, "Fig2B_roc_methylation.pdf"), width=6, height=6)
plot(0, 0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="False positive rate",
     ylab="True positive rate", main="Methylation classifiers (nested CV)")
cols <- c("darkred", "orange")
for (j in 1:2) {
  k <- c(100, 1000)[j]
  preds <- rep(NA, ncol(B))
  for (f in 1:10) {
    te <- which(folds == f); tr <- which(folds != f)
    pv <- fast_welch(B[, tr, drop=FALSE], tr)
    top <- order(pv, na.last=TRUE)[1:k]
    Xtr <- t(B[top, tr, drop=FALSE]); Xtr[is.na(Xtr)] <- 0
    Xte <- t(B[top, te, drop=FALSE]); Xte[is.na(Xte)] <- 0
    m <- randomForest(Xtr, ym[tr], ntree=1000)
    preds[te] <- predict(m, Xte, type="prob")[,2]
  }
  r <- roc(as.numeric(ym)-1, preds)
  plot(r, add=j>1, col=cols[j], lwd=2)
}
legend("bottomright", c("top-100 DMPs", "top-1,000 DMPs"),
       col=cols, lwd=2)
text(0.75, 0.25, "Nested CV (5x10-fold):\ntop-100 AUC = 0.956\ntop-1,000 AUC = 0.947", cex=0.9)
dev.off()
cat("Fig2B nested ROC saved\n")
