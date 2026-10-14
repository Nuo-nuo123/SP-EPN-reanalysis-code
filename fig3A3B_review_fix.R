# Fig3A (title truncation fix) and Fig3B (annotation bars: cluster/age/platform)
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(pheatmap)})
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
figdir <- "D:/claude-study/Spinal Ependymoma/figures"

g80 <- read.csv(file.path(outdir, "gse243240_spine_noob_80.csv"), stringsAsFactors=FALSE)
smap2 <- read.csv(file.path(outdir, "gse243240_287_map.csv"), stringsAsFactors=FALSE)
col_gsm <- function(x) sub("^(GSM[0-9]+)_.*", "\\1", x)
ages <- suppressWarnings(as.numeric(smap2$age[match(g80$gsm, smap2$gsm)]))
plats <- smap2$platform[match(g80$gsm, smap2$gsm)]

# ---- Fig3A: score density, shorter title, extra right margin ----
sc <- g80$score
set.seed(42)
km <- kmeans(sc, centers=2, nstart=25)
boundary <- mean(sort(tapply(sc, km$cluster, mean)))
pdf(file.path(figdir, "Fig3_score_density.pdf"), width=6.2, height=4.6)
p <- ggplot(data.frame(score=sc), aes(x=score)) +
  geom_histogram(aes(y=after_stat(density)), bins=30, fill="steelblue", alpha=0.6) +
  geom_density(alpha=0.3, color="darkblue") +
  geom_vline(xintercept=boundary, linetype="dashed", color="red") +
  annotate("text", x=boundary+0.08, y=0.75, label=paste("k-means boundary:", round(boundary, 2)),
           color="red", hjust=0, size=3.5) +
  labs(title="EPN_SPINE (n = 80, overlap-excluded): subtype scores",
       x="Subtype score (Pearson correlation with A-B DMP direction)", y="Density") +
  theme_minimal(base_size=11) + theme(plot.margin=margin(8, 12, 8, 8))
print(p)
dev.off()
cat("Fig3A saved\n")

# ---- Fig3B: heatmap with annotation bars ----
Vall <- readRDS(file.path(outdir, "gse243240_beta_dmp_noob.rds"))
idx80 <- which(col_gsm(colnames(Vall)) %in% g80$gsm)
Z <- Vall[, idx80, drop=FALSE]
# restrict the displayed rows to the 100 most significant of the 1,000 DMPs
dmp_res <- read.csv(file.path(outdir, "dmp_results.csv"), stringsAsFactors=FALSE)
top100 <- dmp_res$cpg[order(dmp_res$pvalue)][1:100]
top100 <- intersect(top100, rownames(Z))
Z <- Z[top100, , drop=FALSE]
z <- t(scale(t(Z))); z[is.na(z)] <- 0
# annotation uses the manuscript's full-panel (1,000-DMP) cluster assignment
gsm80 <- sub("^(GSM[0-9]+)_.*", "\\1", colnames(Z))
k2 <- as.character(g80$cluster_hier[match(gsm80, g80$gsm)])
ann <- data.frame(
  cluster = ifelse(k2 == "1", "A-like (35)", "B-like (45)"),
  age = ifelse(is.na(ages), "Unknown",
        ifelse(ages < median(ages, na.rm=TRUE), "Younger", "Older")),
  platform = plats,
  row.names = colnames(Z))
ann_colors <- list(
  cluster = c("A-like (35)"="#d62728", "B-like (45)"="#1f77b4"),
  age = c("Younger"="#fff7bc", "Older"="#a6bddb", "Unknown"="grey85"),
  platform = c("GPL13534"="#7fbf7b", "GPL21145"="#af8dc3"))
pdf(file.path(figdir, "Fig3B_spine_heatmap.pdf"), width=11, height=7.5)
pheatmap(z, cluster_rows=TRUE, cluster_cols=TRUE,
         show_rownames=FALSE, show_colnames=FALSE,
         annotation_col=ann, annotation_colors=ann_colors,
         annotation_legend=TRUE, legend=TRUE,
         main="EPN_SPINE (n = 80, overlap-excluded): top 1,000 discovery DMPs, z-scored beta\nWard hierarchical clustering, k = 2 (35 vs 45)",
         color=colorRampPalette(c("navy", "white", "firebrick3"))(100),
         fontsize=9, border_color=NA, annotation_names_col=TRUE)
dev.off()
cat("Fig3B saved\n")
