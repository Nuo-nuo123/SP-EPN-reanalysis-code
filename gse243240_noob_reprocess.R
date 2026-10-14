# GSE243240 noob reprocessing: preprocessNoob + detectionP masking from IDATs,
# then redo (A) overlap screen vs discovery, (B) 94-sample EPN_SPINE clustering,
# (C) 80-sample overlap-excluded clustering + platform subsets + age test.
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages({library(minfi); library(data.table)})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
idat_dir <- "D:/claude-study/Spinal Ependymoma/GEO_data/GSE243240/IDAT"
sink(file.path(outdir, "gse243240_noob.log"), append=FALSE)
cat("== GSE243240 noob reprocessing", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# top-1000 discovery DMPs
dmp <- read.csv(file.path(outdir, "dmp_results.csv"), stringsAsFactors=FALSE)
top_dmps <- dmp$cpg[order(dmp$pvalue)][1:1000]
cat("Top DMPs:", length(top_dmps), "\n")
smap <- read.csv(file.path(outdir, "gse243240_287_map.csv"), stringsAsFactors=FALSE)
gsm_450k <- smap$gsm[smap$platform == "GPL13534"]
gsm_epic <- smap$gsm[smap$platform == "GPL21145"]
cat("450K:", length(gsm_450k), "EPIC:", length(gsm_epic), "\n")

process_platform_noob <- function(gsms, plat_label, out_rds) {
  cat("\n--- noob:", plat_label, "---\n")
  all_files <- list.files(idat_dir, pattern="_Grn.idat$", full.names=TRUE)
  base_names <- sub("_Grn\\.idat$", "", all_files)
  gsm_of_base <- sub(".*/(GSM[0-9]+)_.*", "\\1", base_names)
  bases <- base_names[gsm_of_base %in% gsms]
  cat("Found", length(bases), "of", length(gsms), "samples\n")
  beta_list <- list()
  chunk <- 40
  for (start in seq(1, length(bases), by=chunk)) {
    end <- min(start+chunk-1, length(bases))
    tb <- data.frame(Basename = bases[start:end])
    cat("  Reading", start, "-", end, "...\n")
    rg <- read.metharray.exp(base=NULL, targets=tb, verbose=FALSE, extended=FALSE, force=TRUE)
    cat("  detectionP + preprocessNoob ...\n")
    detP <- detectionP(rg)
    mset <- preprocessNoob(rg)
    beta <- getBeta(mset)
    rm(rg, mset); gc()
    beta[detP > 0.01] <- NA
    rm(detP); gc()
    avail <- intersect(top_dmps, rownames(beta))
    beta_list[[length(beta_list)+1]] <- beta[avail, , drop=FALSE]
    cat("  Chunk:", nrow(beta_list[[length(beta_list)]]), "probes,",
        round(100*mean(is.na(beta_list[[length(beta_list)]])), 2), "% NA\n")
  }
  beta_all <- do.call(cbind, beta_list)
  cat(plat_label, "final:", nrow(beta_all), "x", ncol(beta_all), "\n")
  saveRDS(beta_all, file.path(outdir, out_rds))
  beta_all
}
b450n <- process_platform_noob(gsm_450k, "450K", "gse243240_beta_450k_noob.rds")
bepicn <- process_platform_noob(gsm_epic, "EPIC", "gse243240_beta_epic_noob.rds")
vcom <- intersect(rownames(b450n), rownames(bepicn))
Vall <- cbind(b450n[vcom, , drop=FALSE], bepicn[vcom, , drop=FALSE])
saveRDS(Vall, file.path(outdir, "gse243240_beta_dmp_noob.rds"))
cat("\nCombined noob:", nrow(Vall), "probes x", ncol(Vall), "samples\n")

# ---- discovery betas on the same probes ----
beta72 <- fread(file.path(outdir, "beta_72samples.csv"))
cpg_ids <- beta72[[1]]; B <- as.matrix(beta72[, -1, with=FALSE]); rownames(B) <- cpg_ids
rm(beta72); gc()
Bcom <- B[vcom, , drop=FALSE]
col_gsm <- function(x) sub("^(GSM[0-9]+)_.*", "\\1", x)

# ============================================================
# A. overlap screen: max correlation with any discovery sample
# ============================================================
cat("\n--- overlap screen (noob) ---\n")
maxc <- apply(Vall, 2, function(v) max(cor(v, Bcom, use="pairwise.complete.obs"), na.rm=TRUE))
names(maxc) <- colnames(Vall)
gsm_all <- col_gsm(colnames(Vall))
df <- data.frame(gsm=gsm_all, maxcor=as.numeric(maxc), stringsAsFactors=FALSE)
df$class <- smap$class_v12[match(df$gsm, smap$gsm)]
cat(">0.95:", sum(df$maxcor > 0.95, na.rm=TRUE),
    "| EPN_SPINE >0.95:", sum(df$class == "EPN_SPINE" & df$maxcor > 0.95, na.rm=TRUE),
    "| >0.99:", sum(df$maxcor > 0.99, na.rm=TRUE), "\n")
write.csv(df, file.path(outdir, "round3_sample_overlap_noob.csv"), row.names=FALSE)

# ============================================================
# B/C. EPN_SPINE clustering, 94-sample and 80-sample overlap-excluded
# ============================================================
sp_gsm <- smap$gsm[smap$class_v12 == "EPN_SPINE"]
sp_beta <- Vall[, col_gsm(colnames(Vall)) %in% sp_gsm]
sp <- smap[match(col_gsm(colnames(sp_beta)), smap$gsm), ]
cat("\nEPN_SPINE:", nrow(sp), "probes:", nrow(sp_beta), "\n")

# discovery A-B direction (delta beta) on common probes
labels <- read.csv(file.path(outdir, "GSE242994_rna_labels.csv"), stringsAsFactors=FALSE)
smapd <- read.csv(file.path(outdir, "SPEPN_sentrix_map.csv"), stringsAsFactors=FALSE)
labels$paper_id_clean <- sub("^\\s+", "", labels$paper_id)
smapd$spid_clean <- sub("^\\s+", "", smapd[[1]])
mapped <- merge(labels, smapd, by.x="paper_id_clean", by.y="spid_clean", all.x=TRUE)
subtype <- mapped$Transcriptional_subtype[match(colnames(Bcom), mapped$Sentrix_ID)]
dA <- rowMeans(Bcom[, subtype == "A"], na.rm=TRUE)
dB <- rowMeans(Bcom[, subtype == "B"], na.rm=TRUE)
delta <- dA - dB

run_cluster <- function(Z, tag) {
  z <- t(scale(t(Z))); z[is.na(z)] <- 0
  hc <- hclust(dist(t(z)), method="ward.D2")
  k2 <- cutree(hc, k=2)
  sc <- apply(Z, 2, function(v) cor(v, delta, use="pairwise.complete.obs"))
  set.seed(42)
  km <- kmeans(sc, centers=2, nstart=25)
  suppressPackageStartupMessages(library(mclust))
  ari <- adjustedRandIndex(k2, km$cluster)
  # orient hierarchical clusters by score sign so cluster 1 = positive
  if (mean(sc[k2 == 1], na.rm=TRUE) < mean(sc[k2 == 2], na.rm=TRUE)) {
    k2 <- 3 - k2
  }
  cat(sprintf("\n%s: hier %s | score means %.3f vs %.3f | kmeans %s | ARI %.3f\n",
      tag, paste(table(k2), collapse="/"),
      mean(sc[k2 == 1]), mean(sc[k2 == 2]),
      paste(table(km$cluster), collapse="/"), ari))
  # contingency hier vs kmeans oriented: kmeans cluster with larger mean = positive
  kmp <- if (mean(sc[km$cluster == 1]) > mean(sc[km$cluster == 2])) 1 else 2
  kmo <- ifelse(km$cluster == kmp, 1, 2)
  cat("  contingency (rows hier 1/2, cols kmeans 1/2):\n")
  print(table(k2, kmo))
  cat("  agreement:", sum(k2 == kmo), "of", length(k2), "\n")
  ages <- suppressWarnings(as.numeric(sp$age[match(col_gsm(colnames(Z)), sp$gsm)]))
  ok <- !is.na(ages) & length(unique(k2[!is.na(ages)])) == 2
  if (sum(!is.na(ages)) >= 6) {
    cat("  age:", round(tapply(ages, k2, mean), 1),
        "p =", format(wilcox.test(ages ~ k2)$p.value, digits=3), "\n")
  }
  # platform subsets
  plats <- sp$platform[match(col_gsm(colnames(Z)), sp$gsm)]
  for (pp in unique(plats)) {
    idx <- which(plats == pp)
    if (length(idx) < 10) next
    zz <- z[, idx, drop=FALSE]
    hh <- cutree(hclust(dist(t(zz)), method="ward.D2"), k=2)
    ss <- sc[idx]
    if (mean(ss[hh == 1], na.rm=TRUE) < mean(ss[hh == 2], na.rm=TRUE)) hh <- 3 - hh
    cat(sprintf("  %s-only (n=%d): split %s | means %.3f vs %.3f\n",
        pp, length(idx), paste(table(hh), collapse="/"),
        mean(ss[hh == 1]), mean(ss[hh == 2])))
  }
  invisible(list(k2=k2, sc=sc, km=km))
}

cat("\n=== 94-sample EPN_SPINE (noob) ===")
r94 <- run_cluster(sp_beta, "94-sample")
excl <- df$gsm[df$class == "EPN_SPINE" & df$maxcor > 0.95 & !is.na(df$maxcor)]
keep <- !(col_gsm(colnames(sp_beta)) %in% excl)
cat("\nExcluded overlaps:", sum(!keep), "of", length(keep), "\n")
cat("=== 80-sample overlap-excluded (noob) ===")
sp80_beta <- sp_beta[, keep, drop=FALSE]
r80 <- run_cluster(sp80_beta, "80-sample")
out80 <- data.frame(gsm=col_gsm(colnames(sp80_beta)),
                    score=r80$sc, cluster_hier=r80$k2)
write.csv(out80, file.path(outdir, "gse243240_spine_noob_80.csv"), row.names=FALSE)

cat("\n== noob reprocessing done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
