# NF2 in the deposited data vs the source DEG table (A vs B direction)
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "nf2_replication.log"), append=FALSE)
suppressPackageStartupMessages({library(DESeq2)})
expr_raw <- readRDS(file.path(outdir, "GSE242994_expr.rds"))
expr <- if (is.list(expr_raw)) expr_raw$expr else expr_raw
labels <- read.csv(file.path(outdir, "GSE242994_rna_labels.csv"), stringsAsFactors=FALSE)
common <- intersect(colnames(expr), labels$sample)
y <- labels$Transcriptional_subtype[match(common, labels$sample)]
cnt <- expr[, common]; cnt <- cnt[rowSums(cnt) > 0, ]
dds <- DESeqDataSetFromMatrix(countData=round(cnt), colData=data.frame(subtype=factor(y)), design=~subtype)
dds <- DESeq(dds)
# A vs B to match the source DEG table direction (negative = lower in A)
res <- results(dds, contrast=c("subtype", "A", "B"))
g <- "NF2"
cat("NF2 (A vs B): log2FC", round(res[g, "log2FoldChange"], 3),
    "| padj", format(res[g, "padj"], digits=3),
    "| p", format(res[g, "pvalue"], digits=3), "\n")
# Welch on log2CPM as in the manuscript
cpm <- sweep(cnt, 2, colSums(cnt)/1e6, "/")
lcpm <- log2(cpm + 1)
wt <- t.test(lcpm[g, y == "A"], lcpm[g, y == "B"])
cat("Welch log2CPM: A", round(mean(lcpm[g, y == "A"]), 3),
    "B", round(mean(lcpm[g, y == "B"]), 3),
    "p =", format.pval(wt$p.value), "\n")
sink()
