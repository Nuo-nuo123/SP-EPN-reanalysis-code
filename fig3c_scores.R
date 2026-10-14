# Fig3C: subtype score by methylation class (287 samples)
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(ggplot2)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
figdir <- "D:/claude-study/Spinal Ependymoma/figures"

sc <- read.csv(file.path(outdir, "gse243240_subtype_scores.csv"), stringsAsFactors=FALSE)
# focus on classes with n>=8
keep <- names(table(sc$class_v12))[table(sc$class_v12) >= 8]
sc2 <- sc[sc$class_v12 %in% keep, ]
sc2$class_v12 <- factor(sc2$class_v12, levels=rev(sort(unique(sc2$class_v12))))

pdf(file.path(figdir, "Fig3C_scores_by_class.pdf"), width=8, height=6)
p <- ggplot(sc2, aes(x=score, y=class_v12, fill=class_v12)) +
  geom_boxplot(outlier.size=0.5, alpha=0.7) +
  labs(title="Subtype score by MNP v12.5 class (GSE243240, n=287)",
       x="Correlation with training A-B DMP direction", y="") +
  theme_minimal(base_size=11) +
  theme(legend.position="none")
print(p)
dev.off()
cat("Fig3C saved\n")
