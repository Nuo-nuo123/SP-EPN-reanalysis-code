# 91- and 81-sample sensitivity clusterings on noob betas
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "noob_sens91.log"), append=FALSE)
suppressPackageStartupMessages({library(data.table)})
ov <- read.csv(file.path(outdir, "round3_sample_overlap_noob.csv"), stringsAsFactors=FALSE)
sp_ov <- ov[ov$class == "EPN_SPINE" & !is.na(ov$maxcor), ]
top3 <- sp_ov$gsm[order(-sp_ov$maxcor)][1:3]
cat("top-3 excluded:", top3, "\n")
extra1 <- sp_ov$gsm[sp_ov$maxcor > 0.90 & sp_ov$maxcor <= 0.95]
cat("extra 0.90-0.95 sample:", extra1, "\n")
smap2 <- read.csv(file.path(outdir, "gse243240_287_map.csv"), stringsAsFactors=FALSE)
Vall <- readRDS(file.path(outdir, "gse243240_beta_dmp_noob.rds"))
col_gsm <- function(x) sub("^(GSM[0-9]+)_.*", "\\1", x)
beta <- fread(file.path(outdir, "beta_72samples.csv"))
cpg_ids <- beta[[1]]; B <- as.matrix(beta[, -1, with=FALSE]); rownames(B) <- cpg_ids
rm(beta); gc()
vcom <- intersect(rownames(Vall), rownames(B))
Bcom <- B[vcom, , drop=FALSE]
Vc <- Vall[vcom, , drop=FALSE]
labels <- read.csv(file.path(outdir, "GSE242994_rna_labels.csv"), stringsAsFactors=FALSE)
smapd <- read.csv(file.path(outdir, "SPEPN_sentrix_map.csv"), stringsAsFactors=FALSE)
labels$paper_id_clean <- sub("^\\s+", "", labels$paper_id)
smapd$spid_clean <- sub("^\\s+", "", smapd[[1]])
mapped <- merge(labels, smapd, by.x="paper_id_clean", by.y="spid_clean", all.x=TRUE)
subtype <- mapped$Transcriptional_subtype[match(colnames(Bcom), mapped$Sentrix_ID)]
delta <- rowMeans(Bcom[, subtype=="A"], na.rm=TRUE) - rowMeans(Bcom[, subtype=="B"], na.rm=TRUE)

clus <- function(excl, tag) {
  keep <- which(col_gsm(colnames(Vc)) %in% smap2$gsm[smap2$class_v12 == "EPN_SPINE"] &
                !(col_gsm(colnames(Vc)) %in% excl))
  Z <- Vc[, keep, drop=FALSE]
  sc <- apply(Z, 2, function(v) cor(v, delta, use="pairwise.complete.obs"))
  z <- t(scale(t(Z))); z[is.na(z)] <- 0
  k2 <- cutree(hclust(dist(t(z)), method="ward.D2"), k=2)
  if (mean(sc[k2==1]) < mean(sc[k2==2])) k2 <- 3 - k2
  cat(tag, "n=", length(keep), ": split", paste(table(k2), collapse=" vs "),
      "| score means", paste(round(tapply(sc, k2, mean), 3), collapse=" vs "), "\n")
}
clus(top3, "91-sample (exclude top-3)")
clus(setdiff(sp_ov$gsm[sp_ov$maxcor > 0.95], extra1), "81-sample (exclude 14, include 0.90)")
sink()
