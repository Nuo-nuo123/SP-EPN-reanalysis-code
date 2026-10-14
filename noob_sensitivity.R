# Noob-sensitivity numbers for the manuscript: overlap thresholds, 91-sample
# clustering, gap statistic, spinal-vs-other Wilcoxon, positive-score proportions,
# missing-age Fisher test, MYCN score overlap.
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "noob_sensitivity.log"), append=FALSE)
cat("== Noob sensitivity", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
suppressPackageStartupMessages({library(data.table); library(cluster); library(diptest)})

ov <- read.csv(file.path(outdir, "round3_sample_overlap_noob.csv"), stringsAsFactors=FALSE)
smap2 <- read.csv(file.path(outdir, "gse243240_287_map.csv"), stringsAsFactors=FALSE)
scdf <- read.csv(file.path(outdir, "gse243240_subtype_scores_noob.csv"), stringsAsFactors=FALSE)
sp <- scdf[scdf$class_v12 == "EPN_SPINE", ]
cat("EPN_SPINE maxcor thresholds: >0.99:", sum(sp$gsm %in% ov$gsm[ov$maxcor > 0.99]),
    "| >0.95:", sum(sp$gsm %in% ov$gsm[ov$maxcor > 0.95]),
    "| >0.90:", sum(sp$gsm %in% ov$gsm[ov$maxcor > 0.90]),
    "| >0.80:", sum(sp$gsm %in% ov$gsm[ov$maxcor > 0.80]), "\n")
cat("duplicates platform (EPIC):", sum(smap2$platform[smap2$gsm %in% ov$gsm[ov$maxcor > 0.95 & ov$class == "EPN_SPINE"]] == "GPL21145"),
    "of", sum(ov$class == "EPN_SPINE" & ov$maxcor > 0.95), "\n")

# 91-sample sensitivity: exclude only >0.99
Vall <- readRDS(file.path(outdir, "gse243240_beta_dmp_noob.rds"))
col_gsm <- function(x) sub("^(GSM[0-9]+)_.*", "\\1", x)
beta <- fread(file.path(outdir, "beta_72samples.csv"))
cpg_ids <- beta[[1]]; B <- as.matrix(beta[, -1, with=FALSE]); rownames(B) <- cpg_ids
rm(beta); gc()
vcom <- intersect(rownames(Vall), rownames(B))
Bcom <- B[vcom, , drop=FALSE]; Vc <- Vall[vcom, , drop=FALSE]
labels <- read.csv(file.path(outdir, "GSE242994_rna_labels.csv"), stringsAsFactors=FALSE)
smapd <- read.csv(file.path(outdir, "SPEPN_sentrix_map.csv"), stringsAsFactors=FALSE)
labels$paper_id_clean <- sub("^\\s+", "", labels$paper_id)
smapd$spid_clean <- sub("^\\s+", "", smapd[[1]])
mapped <- merge(labels, smapd, by.x="paper_id_clean", by.y="spid_clean", all.x=TRUE)
subtype <- mapped$Transcriptional_subtype[match(colnames(Bcom), mapped$Sentrix_ID)]
delta <- rowMeans(Bcom[, subtype=="A"], na.rm=TRUE) - rowMeans(Bcom[, subtype=="B"], na.rm=TRUE)
delta <- delta[match(rownames(Vc), names(delta))]

excl91 <- ov$gsm[ov$class == "EPN_SPINE" & ov$maxcor > 0.99 & !is.na(ov$maxcor)]
keep91 <- which(col_gsm(colnames(Vc)) %in% smap2$gsm[smap2$class_v12 == "EPN_SPINE"] &
                !(col_gsm(colnames(Vc)) %in% excl91))
Z91 <- Vc[, keep91, drop=FALSE]
sc91 <- apply(Z91, 2, function(v) cor(v, delta, use="pairwise.complete.obs"))
z91 <- t(scale(t(Z91))); z91[is.na(z91)] <- 0
hc91 <- cutree(hclust(dist(t(z91)), method="ward.D2"), k=2)
if (mean(sc91[hc91==1]) < mean(sc91[hc91==2])) hc91 <- 3 - hc91
cat("\n91-sample: split", paste(table(hc91), collapse=" vs "),
    "| score means", paste(round(tapply(sc91, hc91, mean), 3), collapse=" vs "), "\n")

# positive-score proportion per cluster (80-sample)
g80 <- read.csv(file.path(outdir, "gse243240_spine_noob_80.csv"), stringsAsFactors=FALSE)
cat("positive-score proportion by cluster:", paste(round(tapply(g80$score > 0, g80$cluster_hier, mean)*100, 1), collapse="% vs "), "%\n")

# gap statistic on the 80-sample noob z-scored matrix
idx80 <- which(col_gsm(colnames(Vc)) %in% g80$gsm)
Z <- Vc[, idx80, drop=FALSE]
z <- t(scale(t(Z))); z[is.na(z)] <- 0
cat("\ngap statistic (k=1..4, B=100)...\n")
set.seed(42)
gs <- clusGap(t(z), FUNcluster=function(x, k, ...) list(cluster=kmeans(x, k, nstart=10)$cluster),
              K.max=4, B=100)
print(round(gs$Tab, 3))
cat("optimal k:", maxSE(gs$Tab[,3], gs$Tab[,4], method="globalSEmax"), "\n")

# dip test (80-sample)
cat("dip p:", format(dip.test(g80$score)$p.value, digits=4), "\n")

# spinal vs non-spinal Wilcoxon
scdf$spinal <- ifelse(grepl("SPINE", scdf$class_v12), 1, 0)
wt <- wilcox.test(scdf$score[scdf$spinal==1], scdf$score[scdf$spinal==0])
cat("\nspinal vs non-spinal Wilcoxon p:", format.pval(wt$p.value), "\n")
suppressPackageStartupMessages(library(pROC))
cat("AUC:", round(auc(roc(scdf$spinal, scdf$score)), 3), "\n")

# missing-age Fisher (80-sample)
ages <- suppressWarnings(as.numeric(smap2$age[match(g80$gsm, smap2$gsm)]))
miss <- is.na(ages)
cat("\nage annotations:", sum(!miss), "of", length(ages), "\n")
ft <- fisher.test(table(miss, g80$cluster_hier))
cat("missing-age vs cluster Fisher p:", format.pval(ft$p.value), "\n")

# MYCN score overlap
my <- scdf[scdf$class_v12 == "EPN_SPINE_MYCN", ]
cat("\nMYCN (n=", nrow(my), "): score range", paste(round(range(my$score), 3), collapse=" to "),
    "| mean", round(mean(my$score), 3), "\n")
cat("80-sample cluster score ranges: A-like",
    paste(round(range(g80$score[g80$cluster_hier==1]), 3), collapse=" to "),
    "| B-like", paste(round(range(g80$score[g80$cluster_hier==2]), 3), collapse=" to "), "\n")

cat("\n== sensitivity done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
