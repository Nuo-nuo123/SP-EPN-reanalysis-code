# Review 2026-10-10 analyses:
# (M3) detection rates/CPM for the 8 edgeR + 13 DESeq2 genes;
# (M4) doublet-excluded purified-tumor-nucleus DEG/CytoTRACE sensitivity;
# (M10) 22q-status regrouping sensitivity (bi-allelic vs monoallelic) for DEG + DMP.
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "review1010_analyses.log"), append=FALSE)
cat("== review1010 analyses", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# ================= M3: detection rates of the edgeR/DESeq2 genes =================
expr_raw <- readRDS(file.path(outdir, "GSE242994_expr.rds"))
expr <- if (is.list(expr_raw)) expr_raw$expr else expr_raw
labels <- read.csv(file.path(outdir, "GSE242994_rna_labels.csv"), stringsAsFactors=FALSE)
common <- intersect(colnames(expr), labels$sample)
y <- labels$Transcriptional_subtype[match(common, labels$sample)]
cpm <- sweep(expr[, common], 2, colSums(expr[, common])/1e6, "/")
genes8 <- c("CWH43", "GPR101", "INMT", "RASGRF1", "UNC13A", "CCDC85A", "HAP1", "BASP1")
# DESeq2 13 genes from SM3
sm3 <- read.csv(file.path(outdir, "supply_material_3_full.csv"), stringsAsFactors=FALSE)
genes13 <- sm3$gene[sm3$padj_DESeq2 < 0.05 & !is.na(sm3$padj_DESeq2)]
cat("\n--- M3 detection rates ---\n")
for (g in unique(c(genes8, genes13))) {
  v <- cpm[g, ]
  cat(sprintf("%-10s nonzero %%A %.0f | %%B %.0f | median CPM overall %.2f | max %.1f\n",
      g, 100*mean(v[y=="A"] > 0), 100*mean(v[y=="B"] > 0),
      median(v), max(v)))
}
cat("DESeq2 genes:", paste(genes13, collapse=", "), "\n")

# ================= M4: doublet-excluded purified-tumor-nucleus sensitivity =================
suppressPackageStartupMessages({library(Seurat); library(data.table)})
obj <- readRDS("D:/claude-study/Spinal Ependymoma/analysis/GSE298501/GSE298501_integrated.rds")
db <- read.csv("D:/claude-study/Spinal Ependymoma/analysis/GSE298501/scdblfinder_results.csv",
               stringsAsFactors=FALSE)
rownames(db) <- db$barcode
meta <- obj@meta.data
meta$doublet <- db[rownames(meta), "doublet"]
cts <- read.csv("D:/claude-study/Spinal Ependymoma/analysis/GSE298501/cytotrace_scores.csv",
                stringsAsFactors=FALSE)
cts$doublet <- db[cts$cell, "doublet"]
sub_grp <- sub("^GSM[0-9]+_(SPEPN_A|SPEPN_B)_.*", "\\1", cts$cell)
# CytoTRACE within Ependymal-tumor nuclei — which cells are tumor nuclei? use the integrated meta celltype
meta$cell_bc <- rownames(meta)
tumor_bc <- rownames(meta)[meta$celltype == "Ependymal-tumor"]
cts$celltype <- meta$celltype[match(cts$cell, meta$cell_bc)]
tum <- cts[cts$celltype == "Ependymal-tumor" & sub_grp %in% c("SPEPN_A", "SPEPN_B"), ]
cat("\n--- M4 CytoTRACE within tumor nuclei ---\n")
cat("all: A", sum(tum$subtype_grp <- sub("^GSM[0-9]+_(SPEPN_A|SPEPN_B)_.*", "\\1", tum$cell) == "SPEPN_A"),
    "B", sum(sub("^GSM[0-9]+_(SPEPN_A|SPEPN_B)_.*", "\\1", tum$cell) == "SPEPN_B"), "\n")
tum$grp <- sub("^GSM[0-9]+_(SPEPN_A|SPEPN_B)_.*", "\\1", tum$cell)
tum_clean <- tum[!tum$doublet, ]
w_all <- wilcox.test(CytoTRACE_score ~ grp, data=tum)
w_clean <- wilcox.test(CytoTRACE_score ~ grp, data=tum_clean)
cat(sprintf("CytoTRACE A vs B: all n=%d p=%.3f | doublet-excluded n=%d p=%.3f\n",
    nrow(tum), w_all$p.value, nrow(tum_clean), w_clean$p.value))
cat(sprintf("means: all %.3f/%.3f | clean %.3f/%.3f\n",
    mean(tum$CytoTRACE_score[tum$grp=="SPEPN_A"]), mean(tum$CytoTRACE_score[tum$grp=="SPEPN_B"]),
    mean(tum_clean$CytoTRACE_score[tum_clean$grp=="SPEPN_A"]), mean(tum_clean$CytoTRACE_score[tum_clean$grp=="SPEPN_B"])))

# ================= M10: 22q-status regrouping sensitivity =================
cat("\n--- M10 22q-status regrouping ---\n")
labs72 <- labels[match(common, labels$sample), ]
labs72$group22 <- ifelse(labs72$X22q_loss %in% c("Homozygous", "Homozygous deletion", "bi-allelic", "Biallelic"),
                         "bi-allelic",
                  ifelse(labs72$X22q_loss %in% c("Heterozygous", "monoallelic", "Monoallelic"), "monoallelic", NA))
cat("22q annotation values:", paste(unique(labs72$X22q_loss), collapse=" | "), "\n")
cat("bi-allelic:", sum(labs72$group22 == "bi-allelic", na.rm=TRUE),
    "monoallelic:", sum(labs72$group22 == "monoallelic", na.rm=TRUE), "\n")
# DEG: Welch on log2CPM between the 22q groups
g22 <- labs72$group22
ok22 <- !is.na(g22)
logcpm <- log2(cpm + 1)
kept <- rowSums(cpm >= 1) >= round(0.2 * length(common))
p22 <- apply(logcpm[kept, ], 1, function(v) {
  a <- v[ok22 & g22 == "bi-allelic"]; b <- v[ok22 & g22 == "monoallelic"]
  tryCatch(t.test(a, b)$p.value, error=function(e) NA)
})
padj22 <- p.adjust(p22, "BH")
cat("22q-regrouped DEG: tested", length(p22), "| padj<0.05:", sum(padj22 < 0.05, na.rm=TRUE),
    "| min padj:", format(min(padj22, na.rm=TRUE), digits=3), "\n")
# DMP: Welch on beta between the 22q groups
beta <- fread(file.path(outdir, "beta_72samples.csv"))
cpg_ids <- beta[[1]]; B <- as.matrix(beta[, -1, with=FALSE]); rownames(B) <- cpg_ids
rm(beta); gc()
smap <- read.csv(file.path(outdir, "SPEPN_sentrix_map.csv"), stringsAsFactors=FALSE)
labs72$paper_id_clean <- sub("^\\s+", "", labs72$paper_id)
smap$spid_clean <- sub("^\\s+", "", smap[[1]])
mapped <- merge(labs72, smap, by.x="paper_id_clean", by.y="spid_clean", all.x=TRUE)
sentrix <- mapped$Sentrix_ID[match(colnames(B), mapped$Sentrix_ID)]
# need 22q per column: map via paper id in the same order as B columns
g22_col <- mapped$group22[match(colnames(B), mapped$Sentrix_ID)]
ok22c <- !is.na(g22_col)
p22b <- apply(B[, ok22c, drop=FALSE], 1, function(v) {
  a <- v[g22_col[ok22c] == "bi-allelic"]; b <- v[g22_col[ok22c] == "monoallelic"]
  tryCatch(t.test(a, b)$p.value, error=function(e) NA)
})
padj22b <- p.adjust(p22b, "BH")
cat("22q-regrouped DMPs: tested", length(p22b), "| padj<0.05:", sum(padj22b < 0.05, na.rm=TRUE),
    "| min padj:", format(min(padj22b, na.rm=TRUE), digits=3), "\n")
# overlap with the A/B DMPs
dmp_ab <- read.csv(file.path(outdir, "dmp_significant_annotated.csv"), stringsAsFactors=FALSE)
ov <- sum(padj22b < 0.05 & !is.na(padj22b) & rownames(B) %in% dmp_ab$cpg)
cat("22q-DMPs overlapping the 1,722 A/B DMPs:", ov, "\n")
cat("\n== review1010 analyses done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
