# ML 标签表：GSE242992 72 例 RNA 样本 × xlsx 注释（分子亚型 A/B、NF2 状态、组织学）
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6", TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(), "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages(library(readxl))
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"

# 1) RNA 表达（列名 = SP-EPN ID）
obj <- readRDS(file.path(outdir, "GSE242994_expr.rds"))
rna_cols <- colnames(obj$expr)
cat("RNA columns:", length(rna_cols), "\n")

# 2) xlsx 注释
x <- read_excel("D:/claude-study/Spinal Ependymoma/GEO_data/GSE242994/GSE242994_SP-EPN_annotation.xlsx", sheet=1)
cat("xlsx rows:", nrow(x), "\n")
cat("SP-EPN ID 与 RNA 列匹配数:", sum(x[["SP-EPN ID"]] %in% rna_cols), "\n")
cat("RNA 列中在 xlsx 找到的:", sum(rna_cols %in% x[["SP-EPN ID"]]), "\n")
miss <- setdiff(rna_cols, x[["SP-EPN ID"]])
cat("未匹配 RNA 列（前几个）:", paste(head(miss, 5), collapse=" | "), "\n")

# 3) 标签表
labs <- data.frame(
  sample = rna_cols,
  stringsAsFactors = FALSE
)
xi <- match(rna_cols, x[["SP-EPN ID"]])
labs$Molecular_subtype <- x[["Molecular_subtype"]][xi]
labs$Transcriptional_subtype <- x[["Trasnscriptional_subtype"]][xi]
labs$NF2_status <- x[["NF2 status"]][xi]
labs$Histology <- x[["Histology"]][xi]
labs$Sex <- x[["Sex"]][xi]
labs$Age <- x[["Age (years)"]][xi]
labs$Localization <- x[["Localization"]][xi]
cat("== 72 例 RNA 样本标签分布 ==\n")
print(table(labs$Molecular_subtype, useNA="ifany"))
print(table(labs$NF2_status, useNA="ifany"))
write.csv(labs, file.path(outdir, "GSE242994_rna_labels.csv"), row.names=FALSE)
cat("saved labels\n")
