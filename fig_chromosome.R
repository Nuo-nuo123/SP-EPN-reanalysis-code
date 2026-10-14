# Chromosome distribution bar chart of DMPs (supplementary figure for genomic distribution)
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages({library(ggplot2)})
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
figdir <- "D:/claude-study/Spinal Ependymoma/figures"

dmp_ann <- read.csv(file.path(outdir, "dmp_genomic_annotation_full.csv"), stringsAsFactors=FALSE)
man <- read.csv(file.path(outdir, "epic_manifest_full.csv"), stringsAsFactors=FALSE)

# density: DMPs per 10k array probes per chromosome
chr_dmp <- table(dmp_ann$CHR)
chr_bg <- table(man$CHR)
chrs <- intersect(names(chr_dmp), names(chr_bg))
density <- sapply(chrs, function(ch) 10000 * as.numeric(chr_dmp[ch]) / as.numeric(chr_bg[ch]))
dens_df <- data.frame(chr=chrs, density=as.numeric(density),
                      n_dmp=as.numeric(chr_dmp[chrs]), stringsAsFactors=FALSE)
dens_df <- dens_df[order(dens_df$density, decreasing=TRUE), ]
dens_df$chr <- factor(dens_df$chr, levels=dens_df$chr)
dens_df$sig <- ifelse(dens_df$chr == "22", "chr22 (depleted, OR=0.52)", "Other")

pdf(file.path(figdir, "FigS_chromosome_density.pdf"), width=9, height=5)
p <- ggplot(dens_df, aes(x=chr, y=density, fill=sig)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("chr22 (depleted, OR=0.52)"="#d62728", "Other"="steelblue"), name=NULL) +
  labs(title="DMP density per chromosome (DMPs per 10,000 array probes)",
       x="Chromosome", y="DMP density (per 10k probes)") +
  theme_minimal(base_size=11) +
  theme(legend.position="bottom")
print(p)
dev.off()
cat("chromosome figure saved\n")
