# 最终标签表 v4：确认 Column1 = SP-EPN 249，补齐 72 例
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6", TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(), "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages({library(readxl); library(data.table)})

# 1) 相关匹配
m <- readRDS("D:/claude-study/Spinal Ependymoma/analysis/bulk/match_paper_geo.rds")
w <- which(m$cors > 0.999, arr.ind=TRUE)
map <- data.frame(paper_id=rownames(m$cors)[w[,1]], geo_id=colnames(m$cors)[w[,2]],
                  cor=round(m$cors[w], 5), stringsAsFactors=FALSE)
map <- map[order(map$cor, decreasing=TRUE), ]
map <- map[!duplicated(map$paper_id) & !duplicated(map$geo_id), ]

# 2) Column1 -> SP-EPN 249（论文表3表头错位：基因列被命名为 SP-EPN 249，真实样本列无表头）
g <- map$geo_id[map$paper_id == "Column1"]
cat("Column1 matched GEO:", g, "\n")
# 验证：该 GEO 列与论文所有样本列的最大相关正是 Column1 列
cc <- sort(m$cors[, g], decreasing=TRUE)
cat("top cor for", g, ":", round(cc[1], 5), "(", names(cc)[1], ")  next:", round(cc[2], 5), "(", names(cc)[2], ")\n")
map$paper_id[map$paper_id == "Column1"] <- "SP-EPN 249"
map <- map[order(map$paper_id), ]

# 3) xlsx 标签
x <- read_excel("D:/claude-study/Spinal Ependymoma/GEO_data/GSE242994/GSE242994_SP-EPN_annotation.xlsx", sheet=1)
xi <- match(map$paper_id, x[["SP-EPN ID"]])
cat("xlsx hit:", sum(!is.na(xi)), "/", nrow(map), "\n")
lab <- data.frame(sample=map$geo_id, paper_id=map$paper_id, cor=map$cor, stringsAsFactors=FALSE)
lab$Transcriptional_subtype <- x[["Trasncriptional_subtype"]][xi]
lab$Molecular_subtype <- x[["Molecular_subtype"]][xi]
lab$NF2_status <- x[["NF2 status"]][xi]
lab$Histology <- x[["Histology"]][xi]
lab$Sex <- x[["Sex"]][xi]
lab$Age <- x[["Age (years)"]][xi]
lab$Localization <- x[["Localization"]][xi]
lab$array_type <- x[["array_type"]][xi]
lab$`22q_loss` <- x[["22q loss"]][xi]
lab$Score_mnp_v12.5 <- x[["Score_mnp_v12.5"]][xi]

cat("== 最终分布 ==\n")
print(table(lab$Transcriptional_subtype, useNA="ifany"))
print(table(lab$Molecular_subtype, useNA="ifany"))
print(table(lab$NF2_status, useNA="ifany"))
print(summary(lab$cor))
# 交叉表：转录亚型 vs NF2
print(table(Transcriptional=lab$Transcriptional_subtype, NF2=lab$NF2_status))
write.csv(lab, "D:/claude-study/Spinal Ependymoma/analysis/bulk/GSE242994_rna_labels.csv", row.names=FALSE)
cat("saved rows:", nrow(lab), "\n")
