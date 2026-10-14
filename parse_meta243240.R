# 解析 GSE243240 元数据（表头第 23 行，2,023 行全队列）+ 生存列整理
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6", TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(), "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages(library(readxl))
src <- "D:/claude-study/Spinal Ependymoma/GEO_data/GSE243240/GSE243240_metadata_epn_pohl_130923.xls.gz"
dst <- "D:/claude-study/Spinal Ependymoma/analysis/GSE243240_metadata.csv"

tmp <- "C:/Users/nuo'nuo/AppData/Local/Temp/meta243240.xls"
raw <- readBin(gzfile(src, "rb"), "raw", n=50e6)
writeBin(raw, tmp)
d <- read_excel(tmp, sheet=1, skip=22)
cat("dim:", nrow(d), "x", ncol(d), "\n")
print(names(d))
# 生存列转为数值
for (cn in c("OS [months]", "PFS [months]")) {
  j <- grep(cn, names(d), fixed=TRUE)
  if (length(j)) {
    d[[j]] <- suppressWarnings(as.numeric(gsub(",", ".", as.character(d[[j]]))))
  }
}
# 状态列转 0/1
for (cn in c("OS status", "PFS status")) {
  j <- grep(cn, names(d), fixed=TRUE)
  if (length(j)) {
    cat("==", cn, "==\n"); print(table(d[[j]], useNA="ifany"))
    d[[j]] <- ifelse(grepl("dead|death|yes|^1$|event", as.character(d[[j]]), ignore.case=TRUE), 1,
              ifelse(grepl("alive|no|^0$|censor", as.character(d[[j]]), ignore.case=TRUE), 0, NA))
  }
}
write.csv(d, dst, row.names=FALSE)
cat("saved:", dst, "\n")
j1 <- grep("OS [months]", names(d), fixed=TRUE); j2 <- grep("OS status", names(d), fixed=TRUE)
cat("OS 有效样本:", sum(!is.na(d[[j1]]) & !is.na(d[[j2]])), "\n")
jd <- grep("Patho Diagnosis", names(d), fixed=TRUE)
cat("诊断分布:\n"); print(table(d[[jd]], useNA="ifany"))
cat("分类分布:\n"); print(table(d[[grep("Classification result", names(d), fixed=TRUE)]], useNA="ifany"))
