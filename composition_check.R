# Per-library cell-type composition in the single-nucleus data:
# tumor-cell and astrocyte proportions, SP-EPN-A vs SP-EPN-B (4 vs 4 libraries)
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/GSE298501"
sink(file.path(outdir, "composition_check.log"), append=FALSE)
suppressPackageStartupMessages(library(Seurat))
obj <- readRDS(file.path(outdir, "GSE298501_integrated.rds"))
meta <- obj@meta.data
ct <- as.character(meta$celltype)
lib <- as.character(meta$lib)
sub_grp <- sub("^GSM[0-9]+_(MPE_A|MPE_B|SPEPN_A|SPEPN_B)_.*", "\\1", lib)
tab <- table(data.frame(lib=lib, group=sub_grp, celltype=ct))  # 3D: lib x group x celltype
n_by_lib <- apply(tab, 1, sum)
props <- sweep(tab, 1, n_by_lib, "/")  # keep zeros
libs <- dimnames(tab)$lib
lib_group <- sub("^GSM[0-9]+_(MPE_A|MPE_B|SPEPN_A|SPEPN_B)_.*", "\\1", libs)
groups <- c("MPE_A", "MPE_B", "SPEPN_A", "SPEPN_B")

for (g in groups) {
  lidx <- lib_group == g
  cat(g, "libraries:", sum(lidx), "\n")
  for (cname in c("Ependymal-tumor", "Astrocyte")) {
    v <- props[lidx, g, cname]
    cat("  ", cname, ": mean", round(mean(v), 3),
        "range", paste(round(range(v), 3), collapse = "-"), "\n")
  }
}
cat("\nAll cell types by group (mean prop):\n")
for (g in groups) {
  lidx <- lib_group == g
  m <- sapply(dimnames(tab)$celltype, function(cn) mean(props[lidx, g, cn]))
  cat(" ", g, ":", paste(sprintf("%s %.3f", names(m), m), collapse = " | "), "\n")
}
cat("\nSP-EPN-A vs SP-EPN-B (4 vs 4 libraries):\n")
for (cname in dimnames(tab)$celltype) {
  wa <- props[lib_group == "SPEPN_A", "SPEPN_A", cname]
  wb <- props[lib_group == "SPEPN_B", "SPEPN_B", cname]
  wt <- suppressWarnings(wilcox.test(wa, wb))
  cat(sprintf("  %-18s A %.3f vs B %.3f | p = %s\n",
      cname, mean(wa), mean(wb), format.pval(wt$p.value)))
}
sink()
