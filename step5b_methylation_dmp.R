# Step 5b: DNA methylation DMP analysis for SP-EPN subtypes A vs B
# Beta matrix: GSE242991 (404,256 CpG x 226 samples), 72 RNA-seq samples labeled A/B
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(data.table)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "step5b_methylation_dmp.log"), append=FALSE)
cat("== step5b methylation DMP start", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# ============================================================
# 1. Build sample mapping: paper_id -> Sentrix_ID (from Suppl_Table_1)
# ============================================================
map_csv <- file.path(outdir, "SPEPN_sentrix_map.csv")
stopifnot(file.exists(map_csv))  # pre-extracted by extract_suppl1.py
smap <- read.csv(map_csv, stringsAsFactors=FALSE, check.names=FALSE)
cat("Suppl_Table_1 rows:", nrow(smap), "\n")
cat("Sentrix_ID sample:", head(smap$Sentrix_ID, 3), "\n")

# ============================================================
# 2. Load label file (72 RNA samples with subtype A/B)
# ============================================================
labels <- read.csv(file.path(outdir, "GSE242994_rna_labels.csv"), stringsAsFactors=FALSE)
cat("Label rows:", nrow(labels), "\n")
# paper_id like "SP-EPN 153" -> match to smap$`SP-EPN ID`
labels$paper_id_clean <- sub("^\\s+", "", labels$paper_id)
smap$spid_clean <- sub("^\\s+", "", smap[[1]])
mapped <- merge(labels, smap, by.x="paper_id_clean", by.y="spid_clean", all.x=TRUE)
cat("Mapped to methylation:", sum(!is.na(mapped$Sentrix_ID)), "/", nrow(mapped), "\n")

# ============================================================
# 3. Load beta matrix columns for the 72 samples
# ============================================================
beta_file <- "D:/claude-study/Spinal Ependymoma/GEO_data/GSE242991/GSE242991_20231122_beta_values_normalized.csv.gz"
cat("Loading beta header...\n")
con <- gzfile(beta_file, "rt")
header_line <- readLines(con, n=1)
close(con)
beta_cols <- strsplit(header_line, ",")[[1]]
beta_cols <- gsub('"', '', beta_cols)
cat("Beta matrix columns:", length(beta_cols), "\n")

# Map sentrix -> beta column
sentrix_needed <- mapped$Sentrix_ID[!is.na(mapped$Sentrix_ID)]
beta_col_clean <- sub(".*_([0-9]+_R[0-9]+C[0-9]+)$", "\\1", beta_cols[-1])  # strip GSM prefix
names(beta_col_clean) <- beta_cols[-1]
found <- sentrix_needed %in% beta_col_clean
cat("Sentrix found in beta:", sum(found), "/", length(sentrix_needed), "\n")
if (sum(found) < 50) {
  # also try direct contains matching
  found2 <- sapply(sentrix_needed, function(s) {
    hits <- which(grepl(s, beta_cols[-1], fixed=TRUE))
    if (length(hits)==1) hits else NA
  })
  cat("Contains-match:", sum(!is.na(found2)), "\n")
}

# Build final column selection
sel_cols <- c()
sel_sentrix <- c()
for (s in sentrix_needed) {
  hit <- which(beta_col_clean == s)
  if (length(hit) == 1) {
    sel_cols <- c(sel_cols, beta_cols[-1][hit])
    sel_sentrix <- c(sel_sentrix, s)
  } else {
    hit2 <- which(grepl(s, beta_cols[-1], fixed=TRUE))
    if (length(hit2) == 1) {
      sel_cols <- c(sel_cols, beta_cols[-1][hit2])
      sel_sentrix <- c(sel_sentrix, s)
    }
  }
}
cat("Final selected columns:", length(sel_cols), "\n")

# ============================================================
# 4. Read pre-extracted 72-sample subset (extract_beta_subset.py)
# ============================================================
cat("Loading beta subset (pre-extracted by Python)...\n")
sub_file <- file.path(outdir, "beta_72samples.csv")
stopifnot(file.exists(sub_file))
beta_sub <- fread(sub_file, showProgress=TRUE)
cat("Beta subset dim:", nrow(beta_sub), "x", ncol(beta_sub), "\n")

# First column is CpG ID
cpg_ids <- beta_sub[[1]]
beta_mat <- as.matrix(beta_sub[, -1, with=FALSE])
rownames(beta_mat) <- cpg_ids
rm(beta_sub); gc(verbose=FALSE)

# Reorder to match labels
sample_order <- mapped$Sentrix_ID[match(colnames(beta_mat), mapped$Sentrix_ID)]
subtype <- mapped$Transcriptional_subtype[match(colnames(beta_mat), mapped$Sentrix_ID)]
cat("Samples in beta subset:", ncol(beta_mat), "A:", sum(subtype=="A"), "B:", sum(subtype=="B"), "\n")

# ============================================================
# 5. Quality: remove probes with missing values
# ============================================================
na_frac <- rowMeans(is.na(beta_mat))
cat("Probes with >20% NA:", sum(na_frac > 0.2), "\n")
keep <- na_frac <= 0.2
beta_mat <- beta_mat[keep, ]
cat("Probes kept:", nrow(beta_mat), "\n")

# ============================================================
# 6. DMP: limma-style t-test per CpG with FDR
# ============================================================
cat("Running DMP tests...\n")
grpA <- which(subtype == "A")
grpB <- which(subtype == "B")

# Welch t-test per probe (vectorized via matrixStats or manual row-wise)
n <- nrow(beta_mat)
pvals <- numeric(n)
dmeans <- numeric(n)
for (i in 1:n) {
  x <- beta_mat[i, ]
  a <- x[grpA]; b <- x[grpB]
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (length(a) < 3 || length(b) < 3) { pvals[i] <- NA; next }
  tt <- tryCatch(t.test(a, b), error=function(e) NULL)
  if (is.null(tt)) { pvals[i] <- NA; next }
  pvals[i] <- tt$p.value
  dmeans[i] <- mean(a) - mean(b)
}
padj <- p.adjust(pvals, method="BH")
cat("DMP done. Significant (padj<0.05):", sum(padj < 0.05, na.rm=TRUE), "\n")
cat("|delta beta|>0.2 among sig:", sum(padj < 0.05 & abs(dmeans) > 0.2, na.rm=TRUE), "\n")

dmp_res <- data.frame(
  cpg=rownames(beta_mat),
  delta_beta=dmeans,
  pvalue=pvals,
  padj=padj
)
dmp_res <- dmp_res[order(dmp_res$pvalue), ]
write.csv(dmp_res, file.path(outdir, "dmp_results.csv"), row.names=FALSE)
cat("DMP results saved:", nrow(dmp_res), "probes\n")

# Top 10 DMPs
cat("\nTop 10 DMPs:\n")
print(head(dmp_res[!is.na(dmp_res$padj), ], 10))

# ============================================================
# 7. Methylation-based classifier (RF on top 1000 DMPs)
# ============================================================
cat("\n--- Methylation classifier ---\n")
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
top_n <- min(1000, sum(!is.na(padj)))
top_probes <- dmp_res$cpg[order(dmp_res$pvalue)][1:top_n]
Xm <- t(beta_mat[top_probes, , drop=FALSE])
Xm[is.na(Xm)] <- 0
ym <- as.factor(subtype)

set.seed(42)
folds <- sample(rep(1:10, length.out=nrow(Xm)))
cv_pred <- rep(NA, nrow(Xm))
for (i in 1:10) {
  tr <- folds != i; te <- folds == i
  m <- randomForest(Xm[tr,,drop=FALSE], ym[tr], ntree=1000)
  cv_pred[te] <- predict(m, Xm[te,,drop=FALSE], type="prob")[,2]
}
roc_m <- roc(as.numeric(ym)-1, cv_pred)
cat("Methylation classifier CV AUC (1000 top probes):", round(auc(roc_m), 3), "\n")

# Also try top 100
top100 <- dmp_res$cpg[order(dmp_res$pvalue)][1:100]
Xm100 <- t(beta_mat[top100, , drop=FALSE])
Xm100[is.na(Xm100)] <- 0
set.seed(42)
folds2 <- sample(rep(1:10, length.out=nrow(Xm100)))
cv_pred2 <- rep(NA, nrow(Xm100))
for (i in 1:10) {
  tr <- folds2 != i; te <- folds2 == i
  m <- randomForest(Xm100[tr,,drop=FALSE], ym[tr], ntree=1000)
  cv_pred2[te] <- predict(m, Xm100[te,,drop=FALSE], type="prob")[,2]
}
roc_m100 <- roc(as.numeric(ym)-1, cv_pred2)
cat("Methylation classifier CV AUC (100 top probes):", round(auc(roc_m100), 3), "\n")

# 6-gene expression comparison AUC
cat("(For comparison, 6-gene expression RF CV AUC was 0.514)\n")

saveRDS(list(dmp=dmp_res, roc_full=roc_m, roc_100=roc_m100,
             subtype=subtype, samples=colnames(beta_mat)),
        file.path(outdir, "methylation_dmp_results.rds"))
cat("\n== step5b done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
