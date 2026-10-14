# Pseudo-bulk tumor-level A/B comparison + GSE243240 platform consistency
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/GSE298501"
outdir2 <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "pseudobulk_platform.log"), append=FALSE)
cat("== Pseudo-bulk + platform check", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# ============================================================
# 1. Pseudo-bulk: aggregate Ependymal-tumor cells per library,
#    compare 4 SP-EPN-A libs vs 4 SP-EPN-B libs
# ============================================================
cat("--- 1. Pseudo-bulk tumor-level A vs B ---\n")
obj <- readRDS(file.path(outdir, "GSE298501_integrated.rds"))
DefaultAssay(obj) <- "RNA"
ep_idx <- obj$celltype == "Ependymal-tumor"
obj_ep <- subset(obj, cells=colnames(obj)[ep_idx])
rm(obj); gc()

data_ep <- GetAssayData(obj_ep, layer="data")  # log-normalized
bc <- colnames(obj_ep)
lib <- sub("^((GSM[0-9]+)_(MPE|SPEPN)_[AB]).*", "\\1", bc)
lib_sub <- sub("^GSM[0-9]+_(MPE|SPEPN)_([AB]).*", "\\1-\\2", lib)
libs <- sort(unique(lib_sub))
libs_spa <- grep("SPEPN-A$", libs, value=TRUE)
libs_spb <- grep("SPEPN-B$", libs, value=TRUE)
cat("SP-EPN-A libraries:", paste(libs_spa, collapse=","), "\n")
cat("SP-EPN-B libraries:", paste(libs_spb, collapse=","), "\n")

# pseudo-bulk: mean log-normalized expression per library
pb <- sapply(c(libs_spa, libs_spb), function(L) {
  rowMeans(data_ep[, lib_sub == L, drop=FALSE])
})
colnames(pb) <- c(libs_spa, libs_spb)
grp <- ifelse(colnames(pb) %in% libs_spa, "A", "B")

# tumor-level Wilcoxon on 4 vs 4
pvals_pb <- apply(pb, 1, function(x) {
  a <- x[grp == "A"]; b <- x[grp == "B"]
  if (sd(x) == 0) return(NA)
  tryCatch(wilcox.test(a, b)$p.value, error=function(e) NA)
})
lfc_pb <- apply(pb, 1, function(x) mean(x[grp=="A"]) - mean(x[grp=="B"]))
padj_pb <- p.adjust(pvals_pb, method="BH")
cat("Genes tested:", sum(!is.na(pvals_pb)), "\n")
cat("Sig (padj<0.05, 4v4 pseudo-bulk):", sum(padj_pb < 0.05, na.rm=TRUE), "\n")
cat("Min p:", format(min(pvals_pb, na.rm=TRUE), digits=3), "\n")
cat("Max |mean diff|:", round(max(abs(lfc_pb), na.rm=TRUE), 3), "\n")
cat("Genes with |diff|>0.5:", sum(abs(lfc_pb) > 0.5, na.rm=TRUE), "\n")

# also pseudo-bulk MPE-A vs MPE-B for contrast (15 vs 3 libs)
libs_ma <- grep("MPE-A$", libs, value=TRUE)
libs_mb <- grep("MPE-B$", libs, value=TRUE)
pb_m <- sapply(c(libs_ma, libs_mb), function(L) rowMeans(data_ep[, lib_sub == L, drop=FALSE]))
grp_m <- ifelse(colnames(pb_m) %in% libs_ma, "A", "B")
pvals_m <- apply(pb_m, 1, function(x) {
  a <- x[grp_m=="A"]; b <- x[grp_m=="B"]
  if (sd(x) == 0) return(NA)
  tryCatch(wilcox.test(a, b)$p.value, error=function(e) NA)
})
padj_m <- p.adjust(pvals_m, method="BH")
cat("\nMPE-A vs MPE-B (15v3 pseudo-bulk): sig genes:",
    sum(padj_m < 0.05, na.rm=TRUE), "min p:", format(min(pvals_m, na.rm=TRUE), digits=3), "\n")

res_pb <- data.frame(gene=rownames(pb), pvalue=pvals_pb, padj=padj_pb, mean_diff=lfc_pb)
write.csv(res_pb, file.path(outdir, "pseudobulk_tumor_AB.csv"), row.names=FALSE)
cat("Saved pseudo-bulk results\n")

# ============================================================
# 2. GSE243240 platform consistency: cluster within 450K vs EPIC subsets
# ============================================================
cat("\n--- 2. Platform consistency check ---\n")
beta <- readRDS(file.path(outdir2, "gse243240_beta_dmp.rds"))
smap <- read.csv(file.path(outdir2, "gse243240_287_map.csv"), stringsAsFactors=FALSE)
col_gsm <- sub("^(GSM[0-9]+)_.*", "\\1", colnames(beta))
sp_gsm <- smap$gsm[smap$class_v12 == "EPN_SPINE"]
sp_beta <- beta[, col_gsm %in% sp_gsm]
colnames(sp_beta) <- col_gsm[col_gsm %in% sp_gsm]
platforms <- smap$platform[match(colnames(sp_beta), smap$gsm)]
cat("EPN_SPINE by platform:", paste(table(platforms), collapse=", "), "\n")

cluster_sub <- function(idx, label) {
  x <- sp_beta[, idx]
  xz <- t(scale(t(x))); xz[is.na(xz)] <- 0
  d <- dist(t(xz)); hc <- hclust(d, method="ward.D2")
  k2 <- cutree(hc, k=2)
  db_dir <- read.csv(file.path(outdir2, "dmp_results.csv"))
  db <- db_dir$delta_beta[match(rownames(x), db_dir$cpg)]
  sc <- apply(x, 2, function(v) cor(v, db, use="complete.obs"))
  cat(label, ": n =", sum(idx), ", k2 =", paste(table(k2), collapse="/"),
      ", cluster score means:", paste(round(tapply(sc, k2, mean), 3), collapse=" / "), "\n")
  return(list(k2=k2, score=sc))
}
res_450 <- cluster_sub(platforms == "GPL13534", "450K-only")
res_epic <- cluster_sub(platforms == "GPL21145", "EPIC-only")

cat("\n== done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
