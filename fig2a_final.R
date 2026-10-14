Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages({library(ggplot2)})
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
figdir <- "D:/claude-study/Spinal Ependymoma/figures"
dmp <- read.csv(file.path(outdir, "dmp_results.csv"))
dmp$sig <- ifelse(!is.na(dmp$padj) & dmp$padj < 0.05 & abs(dmp$delta_beta) > 0.2, "Strong DMP",
                  ifelse(!is.na(dmp$padj) & dmp$padj < 0.05, "DMP", "NS"))
sig_idx <- which(dmp$sig != "NS")
ns_idx <- which(dmp$sig == "NS")
set.seed(42)
ns_keep <- sample(ns_idx, 50000)
dmp2 <- dmp[c(sig_idx, ns_keep), ]
dmp2$neglogp <- -log10(pmax(dmp2$padj, 1e-300))
pdf(file.path(figdir, "Fig2A_volcano.pdf"), width=6.5, height=6)
p <- ggplot(dmp2, aes(x=delta_beta, y=neglogp, color=sig)) +
  geom_point(size=0.4, alpha=0.6) +
  scale_color_manual(values=c("Strong DMP"="#b30000", "DMP"="orange", "NS"="grey70"),
                     name=NULL,
                     labels=c("Strong DMP"="Strong DMP (adjusted p < 0.05 AND |delta beta| > 0.2)",
                              "DMP"="DMP (adjusted p < 0.05)",
                              "NS"="Non-significant (adjusted p >= 0.05)")) +
  geom_vline(xintercept=c(-0.2, 0.2), linetype="dashed", color="grey40") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
  xlim(-0.5, 0.5) +
  labs(title="DMPs: subtype A vs B (404,256 CpGs)",
       x="Delta beta (A - B; positive = higher methylation in subtype A)",
       y="-log10(adjusted p)") +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom", legend.text=element_text(size=7.5),
        legend.key.size=unit(0.35, "cm"), legend.margin=margin(t=0, b=0))
print(p)
dev.off()
cat("Fig2A final done\n")
