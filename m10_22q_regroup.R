# M10 corrected: 22q/NF2-status regrouping (bi-allelic proxy = NF2-mutated +
# heterozygous 22q loss; monoallelic proxy = NF2-wt + heterozygous 22q loss).
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "m10_22q_regroup.log"), append=FALSE)
suppressPackageStartupMessages({library(data.table)})
expr_raw <- readRDS(file.path(outdir, "GSE242994_expr.rds"))
expr <- if (is.list(expr_raw)) expr_raw$expr else expr_raw
labels <- read.csv(file.path(outdir, "GSE242994_rna_labels.csv"), stringsAsFactors=FALSE)
common <- intersect(colnames(expr), labels$sample)
labs72 <- labels[match(common, labels$sample), ]
mut <- labs72$NF2_status %in% c("mutation in tumor", "NF2-related schwannomatosis")
het <- labs72$X22q_loss %in% c("heterozygous 22q loss")
g22 <- ifelse(mut & het, "bi-allelic", ifelse(!mut & het, "monoallelic", NA))
cat("bi-allelic (NF2-mut + het 22q):", sum(g22 == "bi-allelic", na.rm=TRUE),
    "| monoallelic (NF2-wt + het 22q):", sum(g22 == "monoallelic", na.rm=TRUE),
    "| other/NA:", sum(is.na(g22)), "\n")
ok22 <- !is.na(g22)
cpm <- sweep(expr[, common], 2, colSums(expr[, common])/1e6, "/")
logcpm <- log2(cpm + 1)
kept <- rowSums(cpm >= 1) >= round(0.2 * length(common))
p22 <- apply(logcpm[kept, ], 1, function(v) {
  a <- v[ok22 & g22 == "bi-allelic"]; b <- v[ok22 & g22 == "monoallelic"]
  tryCatch(t.test(a, b)$p.value, error=function(e) NA)
})
padj22 <- p.adjust(p22, "BH")
cat("DEG (22q regrouping): tested", length(p22), "| padj<0.05:", sum(padj22 < 0.05, na.rm=TRUE),
    "| min padj:", format(min(padj22, na.rm=TRUE), digits=3), "\n")
# DMP
beta <- fread(file.path(outdir, "beta_72samples.csv"))
cpg_ids <- beta[[1]]; B <- as.matrix(beta[, -1, with=FALSE]); rownames(B) <- cpg_ids
rm(beta); gc()
smap <- read.csv(file.path(outdir, "SPEPN_sentrix_map.csv"), stringsAsFactors=FALSE)
labs72$paper_id_clean <- sub("^\\s+", "", labs72$paper_id)
smap$spid_clean <- sub("^\\s+", "", smap[[1]])
mapped <- merge(labs72, smap, by.x="paper_id_clean", by.y="spid_clean", all.x=TRUE)
g22_col <- mapped$X22q_loss[match(colnames(B), mapped$Sentrix_ID)]
mut_col <- mapped$NF2_status[match(colnames(B), mapped$Sentrix_ID)]
mut_c <- mut_col %in% c("mutation in tumor", "NF2-related schwannomatosis")
het_c <- g22_col %in% c("heterozygous 22q loss")
g22c <- ifelse(mut_c & het_c, "bi-allelic", ifelse(!mut_c & het_c, "monoallelic", NA))
okc <- !is.na(g22c)
p22b <- apply(B[, okc, drop=FALSE], 1, function(v) {
  a <- v[g22c[okc] == "bi-allelic"]; b <- v[g22c[okc] == "monoallelic"]
  tryCatch(t.test(a, b)$p.value, error=function(e) NA)
})
padj22b <- p.adjust(p22b, "BH")
cat("DMP (22q regrouping): tested", length(p22b), "| padj<0.05:", sum(padj22b < 0.05, na.rm=TRUE),
    "| min padj:", format(min(padj22b, na.rm=TRUE), digits=3), "\n")
dmp_ab <- read.csv(file.path(outdir, "dmp_significant_annotated.csv"), stringsAsFactors=FALSE)
ov <- sum(padj22b < 0.05 & !is.na(padj22b) & rownames(B) %in% dmp_ab$cpg)
cat("22q-DMPs overlapping the 1,722 A/B DMPs:", ov, "\n")
# how many A/B DMPs are nominally (p<0.05) 22q-DMPs?
nom22 <- p22b < 0.05 & !is.na(p22b)
cat("A/B DMPs with nominal 22q-group p<0.05:", sum(nom22 & rownames(B) %in% dmp_ab$cpg), "of", nrow(dmp_ab), "\n")
sink()
