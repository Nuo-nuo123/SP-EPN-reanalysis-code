# 构建 GSE242994 RNA-seq（72 例 total RNA）表达矩阵
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6", TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(), "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
library(data.table)
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
gz <- "D:/claude-study/Spinal Ependymoma/GEO_data/GSE242994/GSE242994-GPL24676_series_matrix.txt.gz"

txt <- readLines(gzfile(gz), n=20000)
beg <- which(txt == "!series_matrix_table_begin")
dat <- read.table(gzfile(gz), skip=beg, header=TRUE, sep="\t", quote="\"",
                  comment.char="", check.names=FALSE, stringsAsFactors=FALSE)
cat("raw:", nrow(dat), "genes x", ncol(dat)-1, "samples\n")
expr <- as.matrix(dat[, -1]); rownames(expr) <- dat[[1]]
suppressWarnings(storage.mode(expr) <- "double")

# 样本注释
meta <- txt[1:(beg-1)]
g <- grepl("^!Sample_", meta)
m2 <- strsplit(gsub('"', '', meta[g]), "\t")
snames <- sapply(m2, function(z) z[1])
vals <- lapply(m2, function(z) z[-1])
sm <- data.frame(sample=colnames(expr), stringsAsFactors=FALSE)
for (i in seq_along(vals)) {
  if (length(vals[[i]]) == ncol(expr)) sm[[gsub("!Sample_", "", snames[i])]] <- vals[[i]]
}
mol <- sm[["molecule_ch1"]]
cat("molecules:", paste(names(table(mol)), table(mol), sep="=", collapse=" "), "\n")

keep <- mol == "total RNA"
expr <- expr[, keep, drop=FALSE]
sm <- sm[keep, , drop=FALSE]
cat("RNA samples:", ncol(expr), "\n")

saveRDS(list(expr=expr, samples=sm), file.path(outdir, "GSE242994_expr.rds"))
write.csv(sm, file.path(outdir, "GSE242994_samples.csv"), row.names=FALSE)
cat("== GSE242994 RNA build done ==\n")
