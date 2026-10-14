Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/GSE298501"
sink(file.path(outdir, "cytotrace_by_group.log"), append=FALSE)
cat("== CytoTRACE group comparison", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# Load CytoTRACE scores
ct <- read.csv(file.path(outdir, "cytotrace_scores.csv"))
cat("CytoTRACE rows:", nrow(ct), "\n")
# cell barcodes contain library prefix like GSM9015502_MPE_A_10_
# group = MPE or SPEPN
ct$group <- ifelse(grepl("_MPE_", ct$cell), "MPE",
            ifelse(grepl("_SPEPN_", ct$cell), "SPEPN", "unknown"))
ct$group2 <- sub("^(GSM[0-9]+)_(MPE|SPEPN)_.*", "\\2", ct$cell)
cat("Group counts:\n"); print(table(ct$group2))

for (g in c("MPE", "SPEPN")) {
  idx <- ct$group2 == g
  cat(sprintf("%s: n=%d, mean=%.3f, sd=%.3f, median=%.3f\n",
              g, sum(idx), mean(ct$CytoTRACE_score[idx]),
              sd(ct$CytoTRACE_score[idx]), median(ct$CytoTRACE_score[idx])))
}
wt <- wilcox.test(ct$CytoTRACE_score[ct$group2=="MPE"],
                  ct$CytoTRACE_score[ct$group2=="SPEPN"])
cat("Wilcoxon MPE vs SPEPN p:", format.pval(wt$p.value), "\n")

# Effect size
library(effsize)
cd <- cohen.d(ct$CytoTRACE_score[ct$group2=="MPE"],
              ct$CytoTRACE_score[ct$group2=="SPEPN"])
cat("Cohen's d:", round(cd$estimate, 3), "\n")

cat("\n== done ==")
sink()
