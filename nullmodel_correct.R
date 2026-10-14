# Corrected null models on the 20,000-probe random pool (noob): 100 draws of
# 1,000 probes -> (a) age association rate, (b) Ward-vs-kmeans ARI distribution.
# DMP-panel reference recomputed on the 80-sample subset.
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "nullmodel_correct.log"), append=FALSE)
cat("== corrected null models", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
suppressPackageStartupMessages({library(data.table); library(mclust)})

V20 <- readRDS(file.path(outdir, "gse243240_beta_random20k_noob.rds"))
Vdmp <- readRDS(file.path(outdir, "gse243240_beta_dmp_noob.rds"))
g80 <- read.csv(file.path(outdir, "gse243240_spine_noob_80.csv"), stringsAsFactors=FALSE)
smap2 <- read.csv(file.path(outdir, "gse243240_287_map.csv"), stringsAsFactors=FALSE)
col_gsm <- function(x) sub("^(GSM[0-9]+)_.*", "\\1", x)
ages <- suppressWarnings(as.numeric(smap2$age[match(col_gsm(colnames(V20)), smap2$gsm)]))
ok_age <- !is.na(ages)

beta <- fread(file.path(outdir, "beta_72samples.csv"))
cpg_ids <- beta[[1]]; B <- as.matrix(beta[, -1, with=FALSE]); rownames(B) <- cpg_ids
rm(beta); gc()
labels <- read.csv(file.path(outdir, "GSE242994_rna_labels.csv"), stringsAsFactors=FALSE)
smapd <- read.csv(file.path(outdir, "SPEPN_sentrix_map.csv"), stringsAsFactors=FALSE)
labels$paper_id_clean <- sub("^\\s+", "", labels$paper_id)
smapd$spid_clean <- sub("^\\s+", "", smapd[[1]])
mapped <- merge(labels, smapd, by.x="paper_id_clean", by.y="spid_clean", all.x=TRUE)
subtype <- mapped$Transcriptional_subtype[match(colnames(B), mapped$Sentrix_ID)]

idx80 <- which(col_gsm(colnames(V20)) %in% g80$gsm)
cat("80-sample columns in pool:", length(idx80), "\n")
ages80 <- ages[idx80]
ok_age <- !is.na(ages80)
cat("age annotations in the 80 samples:", sum(ok_age), "\n")

run_panel <- function(M) {
  # M: probes x 80 samples (already subset); returns list(age_p, ari)
  dM <- rowMeans(B[rownames(M), subtype=="A"], na.rm=TRUE) -
        rowMeans(B[rownames(M), subtype=="B"], na.rm=TRUE)
  z <- t(scale(t(M))); z[is.na(z)] <- 0
  k2 <- cutree(hclust(dist(t(z)), method="ward.D2"), k=2)
  sc <- apply(M, 2, function(v) cor(v, dM, use="pairwise.complete.obs"))
  set.seed(42)
  km <- kmeans(sc, centers=2, nstart=25)
  ari <- adjustedRandIndex(k2, km$cluster)
  wt <- suppressWarnings(wilcox.test(ages80[ok_age] ~ k2[ok_age])$p.value)
  list(age_p = wt, ari = ari, k2 = k2, sc = sc)
}

# DMP panel reference on the 80 samples
M_dmp <- Vdmp[, idx80, drop=FALSE]
ref <- run_panel(M_dmp)
cat("DMP panel (80 samples): Ward vs kmeans ARI =", round(ref$ari, 3),
    "| age p =", format.pval(ref$age_p), "\n")

# 100 random 1,000-probe draws from the 20,000 pool
set.seed(42)
age_ps <- numeric(100); aris <- numeric(100)
for (i in 1:100) {
  rr <- sample(nrow(V20), 1000)
  M <- V20[rr, idx80, drop=FALSE]
  r <- run_panel(M)
  age_ps[i] <- r$age_p; aris[i] <- r$ari
}
cat("\nrandom-probe age null: significant (p<0.05):", sum(age_ps < 0.05), "of 100 | median p:",
    format(median(age_ps), digits=3), "\n")
cat("random-probe ARI: median", round(median(aris), 3), "| mean", round(mean(aris), 3),
    "| max", round(max(aris), 3), "| >=0.9:", sum(aris >= 0.9), "of 100\n")
cat("== done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
