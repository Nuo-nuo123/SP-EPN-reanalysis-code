# scDblFinder per-library doublet detection (background)
# Seurat v5 object with split layers: counts.<lib> per library -> access layers directly
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/GSE298501"
sink(file.path(outdir, "scdblfinder.log"), append=FALSE)
cat("== scDblFinder run", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
suppressPackageStartupMessages({library(scDblFinder); library(Seurat)})

obj <- readRDS(file.path(outdir, "GSE298501_merged_raw.rds"))
cat("object:", ncol(obj), "nuclei,", nrow(obj), "genes\n")
assay_obj <- obj[[DefaultAssay(obj)]]
lyr_names <- names(assay_obj@layers)
cat("layers:", length(lyr_names), "\n")
bcs <- colnames(obj)
libs <- unique(sub("^([^_]+_[^_]+_[^_]+_[^_]+)_.*", "\\1", bcs))
cat("libraries:", length(libs), "\n")
all_res <- list()
for (li in libs) {
  idx <- which(sub("^([^_]+_[^_]+_[^_]+_[^_]+)_.*", "\\1", bcs) == li)
  cat(li, "nuclei:", length(idx), "\n")
  lyr <- paste0("counts.", li)
  cnt <- assay_obj@layers[[lyr]]
  if (is.null(cnt)) {
    cat("  !! layer not found, skipping\n")
    next
  }
  sce <- SingleCellExperiment(list(counts=as(cnt, "sparseMatrix")))
  set.seed(123)
  res <- tryCatch(scDblFinder(sce, clusters=FALSE), error=function(e) {
    cat("  !! scDblFinder error:", conditionMessage(e), "\n"); NULL
  })
  if (!is.null(res)) {
    all_res[[li]] <- data.frame(lib=li, barcode=bcs[idx],
                                doublet=as.logical(res$scDblFinder.class == "doublet"),
                                score=res$scDblFinder.score)
  }
}
out <- do.call(rbind, all_res)
write.csv(out, file.path(outdir, "scdblfinder_results.csv"), row.names=FALSE)
cat("overall doublet rate:", sprintf("%.2f%%", 100*mean(out$doublet)), "\n")
out$entity <- sub("^GSM[0-9]+_(MPE|SPEPN)_.*", "\\1", out$barcode)
print(tapply(out$doublet, out$entity, mean))
cat("== done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
