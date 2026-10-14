Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("D:/claude-study/Spinal Ependymoma/analysis/GSE298501/GSE298501_integrated.rds")
meta <- obj@meta.data
lib <- as.character(meta$lib)
ent <- sub("^GSM[0-9]+_(MPE|SPEPN)_.*", "\\1", lib)
tab <- table(data.frame(lib=lib, ent=ent, celltype=as.character(meta$celltype)))
n_by_lib <- apply(tab, 1, sum)
props <- sweep(tab, 1, n_by_lib, "/")
libs <- dimnames(props)$lib
lib_ent <- sub("^GSM[0-9]+_(MPE|SPEPN)_.*", "\\1", libs)
cat("ent dimnames:", paste(dimnames(props)$ent, collapse=","), "\n")
for (e in c("MPE", "SPEPN")) {
  v <- props[lib_ent == e, e, "Ependymal-tumor"]
  cat(e, "tumor-cell proportion: median", round(median(v)*100, 1), "% | mean", round(mean(v)*100, 1),
      "% | range", paste(round(range(v)*100, 1), collapse="-"), "% | n_lib", length(v), "\n")
}
