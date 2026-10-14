Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6", TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(), "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages(library(readxl))
x <- read_excel("D:/claude-study/Spinal Ependymoma/GEO_data/GSE242994/GSE242994_SP-EPN_annotation.xlsx", sheet=1)
cat("RNA seq col:\n"); print(table(as.character(x[["RNA seq"]]), useNA="ifany"))
cat("DNA seq col:\n"); print(table(as.character(x[["DNA seq"]]), useNA="ifany"))
# 有 RNA seq 的行
rna_rows <- which(!is.na(x[["RNA seq"]]) & x[["RNA seq"]] != "")
cat("RNA-seq flagged rows:", length(rna_rows), "\n")
if (length(rna_rows)) {
  y <- x[rna_rows, ]
  cat("subtype among RNA-flagged:\n")
  print(table(y[["Trasncriptional_subtype"]], useNA="ifany"))
  print(table(y[["Molecular_subtype"]], useNA="ifany"))
  cat("RNA-flagged SP-EPN IDs:", paste(head(y[["SP-EPN ID"]], 10), collapse=", "), "...\n")
  # 与相关匹配的 paper_id 对照
  m <- readRDS("D:/claude-study/Spinal Ependymoma/analysis/bulk/match_paper_geo.rds")
  w <- which(m$cors > 0.999, arr.ind=TRUE)
  map <- data.frame(paper_id=rownames(m$cors)[w[,1]], geo_id=colnames(m$cors)[w[,2]], stringsAsFactors=FALSE)
  map <- map[!duplicated(map$paper_id), ]
  cat("mapped paper ids (n=", nrow(map), "):\n", sep="")
  cat("  overlap with RNA-flagged xlsx rows:", sum(map$paper_id %in% y[["SP-EPN ID"]]), "\n")
  cat("  mapped ids NOT in xlsx:", paste(head(setdiff(map$paper_id, x[["SP-EPN ID"]]), 8), collapse=", "), "\n")
}
