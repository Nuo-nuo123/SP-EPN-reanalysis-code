# Exploratory single-cell comparison: SP-EPN-A vs SP-EPN-B (library-level labels)
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(data.table)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/GSE298501"
sink(file.path(outdir, "singlecell_AB_comparison.log"), append=FALSE)
cat("== Single-cell A/B comparison", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# ============================================================
# 1. CytoTRACE by subtype (library suffix: SPEPN_A_* vs SPEPN_B_*)
# ============================================================
ct <- read.csv(file.path(outdir, "cytotrace_scores.csv"))
ct$entity <- sub("^GSM[0-9]+_(MPE|SPEPN)_.*", "\\1", ct$cell)
ct$subtype <- ifelse(grepl("_SPEPN_A_", ct$cell), "SP-EPN-A",
              ifelse(grepl("_SPEPN_B_", ct$cell), "SP-EPN-B",
              ifelse(grepl("_MPE_A_", ct$cell), "MPE-A",
              ifelse(grepl("_MPE_B_", ct$cell), "MPE-B", "other"))))
cat("CytoTRACE by subtype:\n")
print(table(ct$subtype))
for (s in c("SP-EPN-A", "SP-EPN-B", "MPE-A", "MPE-B")) {
  idx <- ct$subtype == s
  if (sum(idx) > 10) {
    cat(sprintf("%s: n=%d mean=%.3f sd=%.3f\n", s, sum(idx),
                mean(ct$CytoTRACE_score[idx]), sd(ct$CytoTRACE_score[idx])))
  }
}
# A vs B within SP-EPN
a <- ct$CytoTRACE_score[ct$subtype == "SP-EPN-A"]
b <- ct$CytoTRACE_score[ct$subtype == "SP-EPN-B"]
if (length(a) > 10 && length(b) > 10) {
  wt <- wilcox.test(a, b)
  cat(sprintf("SP-EPN-A vs SP-EPN-B: p=%.4g\n", wt$p.value))
}
# MPE-A vs MPE-B
a2 <- ct$CytoTRACE_score[ct$subtype == "MPE-A"]
b2 <- ct$CytoTRACE_score[ct$subtype == "MPE-B"]
if (length(a2) > 10 && length(b2) > 10) {
  wt2 <- wilcox.test(a2, b2)
  cat(sprintf("MPE-A vs MPE-B: p=%.4g\n", wt2$p.value))
}

# ============================================================
# 2. Library-level counts by subtype (from lib_summary)
# ============================================================
lib <- read.csv(file.path(outdir, "lib_summary.csv"))
lib$subtype <- sub("^GSM[0-9]+_(MPE_A|MPE_B|SPEPN_A|SPEPN_B).*", "\\1", lib$prefix)
lib$subtype <- gsub("_", "-", lib$subtype)
cat("\nLibrary summary by subtype:\n")
print(aggregate(cells ~ subtype, data=lib, FUN=function(x) c(n=length(x), sum=sum(x))))

cat("\n== Single-cell A/B done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
