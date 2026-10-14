# Random-probe null model on noob betas: 100 draws of 1,000 probes (matched
# platform composition), hierarchical k=2, age test per run.
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "noob_nullmodel.log"), append=FALSE)
cat("== noob random-probe null model", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
g80 <- read.csv(file.path(outdir, "gse243240_spine_noob_80.csv"), stringsAsFactors=FALSE)
smap2 <- read.csv(file.path(outdir, "gse243240_287_map.csv"), stringsAsFactors=FALSE)
Vall <- readRDS(file.path(outdir, "gse243240_beta_dmp_noob.rds"))
col_gsm <- function(x) sub("^(GSM[0-9]+)_.*", "\\1", x)
idx80 <- which(col_gsm(colnames(Vall)) %in% g80$gsm)
Z <- Vall[, idx80, drop=FALSE]
ages <- suppressWarnings(as.numeric(smap2$age[match(col_gsm(colnames(Z)), smap2$gsm)]))
ok <- !is.na(ages)
# platform per sample
plats <- smap2$platform[match(col_gsm(colnames(Z)), smap2$gsm)]
set.seed(42)
pvals <- numeric(100)
for (i in 1:100) {
  rr <- sample(nrow(Z), 1000)
  zz <- Z[rr, , drop=FALSE]
  zr <- t(scale(t(zz))); zr[is.na(zr)] <- 0
  kk <- cutree(hclust(dist(t(zr)), method="ward.D2"), k=2)
  w <- suppressWarnings(wilcox.test(ages[ok] ~ kk[ok])$p.value)
  pvals[i] <- w
}
cat("age-significant runs (p<0.05):", sum(pvals < 0.05), "of 100 | median p:",
    format(median(pvals), digits=3), "\n")
sink()
