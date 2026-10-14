# 通过计数矩阵把 GEO 的 DKFZ 样本匹配到论文 SP-EPN ID
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6", TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(), "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages({library(readxl); library(data.table); library(org.Hs.eg.db)})

xlsx <- "D:/claude-study/Spinal Ependymoma/GEO_data/GSE242994/401_2023_2668_MOESM2_ESM.xlsx"
t3 <- read_excel(xlsx, sheet="Suppl_Table_3_RNASeq_counts")
pg <- as.character(t3[[1]])
Mp <- as.matrix(t3[, -1]); rownames(Mp) <- pg
cat("paper matrix:", nrow(Mp), "x", ncol(Mp), "\n")

geo <- fread("D:/claude-study/Spinal Ependymoma/GEO_data/GSE242994/GSE242992_SP-EPN_raw_RNA_countMatrix_72samples.csv.gz", quote="\"")
eg <- as.character(geo[[1]])
Gmat <- as.matrix(geo[, -1, with=FALSE]); rownames(Gmat) <- eg
cat("geo matrix:", nrow(Gmat), "x", ncol(Gmat), "\n")

# 1) 先按 ENSEMBL 精确交集
ens_common <- intersect(rownames(Mp), eg)
cat("ENSEMBL direct overlap:", length(ens_common), "\n")

sym <- suppressMessages(mapIds(org.Hs.eg.db, keys=eg, keytype="ENSEMBL", column="SYMBOL", multiVals="first"))
sym <- as.character(sym)
ok <- !is.na(sym)
sym2ens <- split(eg[ok], sym[ok])   # symbol -> ensembl 向量
cat("geo ENSEMBL with symbol:", sum(ok), " unique symbols:", length(sym2ens), "\n")
common <- intersect(pg, names(sym2ens))
cat("symbol overlap:", length(common), "\n")

# GEO 端按 symbol 汇总（同 symbol 多 ENSEMBL 求和）
Mg <- t(vapply(common, function(s) {
  colSums(Gmat[sym2ens[[s]], , drop=FALSE])
}, numeric(ncol(Gmat))))
rownames(Mg) <- common
Mp <- Mp[common, , drop=FALSE]
cat("matrices built\n")

# 2) 精确整数匹配
exact <- matrix(FALSE, ncol(Mp), ncol(Mg))
for (i in seq_len(ncol(Mp))) for (j in seq_len(ncol(Mg))) {
  exact[i,j] <- all(Mp[,i] == Mg[,j])
}
cat("exact matches:", sum(exact), "\n")
if (sum(exact) > 0) {
  w <- which(exact, arr.ind=TRUE)
  map <- data.frame(paper=colnames(Mp)[w[,1]], geo=colnames(Mg)[w[,2]], stringsAsFactors=FALSE)
  print(map)
  write.csv(map, "D:/claude-study/Spinal Ependymoma/analysis/bulk/paper_geo_map.csv", row.names=FALSE)
}

# 3) 相关度（无论是否精确匹配，看全局最大值以确认 1:1）
cors <- cor(Mp, Mg)
cat("cor max per paper sample: min", round(min(apply(cors,1,max)),5), " max", round(max(cors),5), "\n")
w2 <- which(cors > 0.999, arr.ind=TRUE)
cat("cor>0.999 pairs:", nrow(w2), "\n")
saveRDS(list(Mp=Mp, Mg=Mg, cors=cors, exact=exact), "D:/claude-study/Spinal Ependymoma/analysis/bulk/match_paper_geo.rds")
cat("saved\n")
