# Step1: GSE298501 snRNA 预处理 —— 加载26文库 + QC + 合并
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6", TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(), "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))

rawdir <- "D:/claude-study/Spinal Ependymoma/GEO_data/GSE298501/RAW"
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/GSE298501"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
logf <- file.path(outdir, "step1_qc.log")
sink(logf, split=TRUE); cat("== step1 start", format(Sys.time()), "==\n")

fs <- list.files(rawdir, pattern="_matrix[.]mtx[.]gz$", full.names=TRUE)
objs <- list(); meta <- list()
for (f in fs) {
  prefix <- sub("_matrix[.]mtx[.]gz$", "", basename(f))
  gsm <- sub("_.*", "", prefix)
  grp <- ifelse(grepl("MPE", prefix), "MPE", "SPEPN")
  cat("loading", prefix, "...\n")
  mat <- ReadMtx(mtx=f, cells=file.path(rawdir, paste0(prefix, "_barcodes.tsv.gz")),
                 features=file.path(rawdir, paste0(prefix, "_features.tsv.gz")))
  colnames(mat) <- paste0(prefix, "_", colnames(mat))
  obj <- CreateSeuratObject(counts=mat, project=prefix, min.cells=3, min.features=200)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern="^MT-")
  obj$lib <- prefix; obj$gsm <- gsm; obj$group <- grp
  objs[[prefix]] <- obj
  meta[[prefix]] <- data.frame(prefix=prefix, gsm=gsm, group=grp,
                               cells=ncol(obj), genes=nrow(obj),
                               median_nFeature=median(obj$nFeature_RNA),
                               median_pctMT=median(obj$percent.mt))
  cat("  ->", ncol(obj), "cells\n")
}
merged <- merge(objs[[1]], y=objs[-1])
cat("merged:", ncol(merged), "cells,", nrow(merged), "genes\n")

# QC 过滤（FFPE snRNA 宽松阈值；记录过滤前后）
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
cat("after QC:", ncol(merged), "cells\n")

saveRDS(merged, file.path(outdir, "GSE298501_merged_raw.rds"))
write.csv(do.call(rbind, meta), file.path(outdir, "lib_summary.csv"), row.names=FALSE)
cat("== step1 done", format(Sys.time()), "==\n")
sink()
