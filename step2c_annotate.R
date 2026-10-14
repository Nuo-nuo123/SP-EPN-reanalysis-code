# step2c：细胞亚群注释 + Ro/e 计算
# marker 来源：
#   - 室管膜瘤肿瘤细胞/脑区细胞: Gojo et al. Cancer Cell 2020 (PMID 32888491, 儿童后颅窝室管膜瘤 scRNA)
#   - 胶质/神经元/血管: 通用 CNS 标志物 (Allen Brain Atlas; PanglaoDB)
#   - 免疫细胞: Zhang et al. Genome Biol 2019 (PMID 30674793) 免疫细胞标志物综述
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6", TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(), "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages(library(Seurat))

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/GSE298501"
logf <- file.path(outdir, "step2c_annotate.log")
sink(logf, split=TRUE); cat("== step2c start", format(Sys.time()), "==\n")

obj <- readRDS(file.path(outdir, "GSE298501_integrated.rds"))
cat("cells:", ncol(obj), "\n")

# ---- 经典 marker 面板（每类注释来源）----
markers <- list(
  `Ependymal-tumor` = c("SPARCL1","AQP4","GFAP","CD44","EPCAM","CD24","FABP7","CST3"), # Gojo2020 室管膜瘤细胞 + 室管膜细胞
  Neuron = c("RBFOX3","SNAP25","SYT1","GRIN1"),      # CNS 神经元
  Oligodendrocyte = c("MBP","MOG","PLP1","OLIG2","CNP"),  # 少突
  OPC = c("PDGFRA","CSPG4","OLIG1","BCAN"),          # OPC
  Astrocyte = c("AQP4","GFAP","SLC1A3","ALDH1L1"),   # 星形（与肿瘤重叠，注释时看共表达）
  Microglia = c("PTPRC","CSF1R","C1QA","CX3CR1","P2RY12","TMEM119"),
  Macrophage = c("CD68","CD163","MRC1","LYZ"),
  `T cell` = c("CD3D","CD3E","CD8A","IL7R","NKG7"),
  B_cell = c("MS4A1","CD79A","CD79B","IGHM"),
  Plasma = c("MZB1","XBP1","JCHAIN"),
  Endothelial = c("PECAM1","VWF","CLDN5","FLT1","CDH5"),
  Pericyte = c("PDGFRB","RGS5","ACTA2","DES"),
  Fibroblast = c("COL1A1","DCN","LUM","COL3A1"),
  Schwann = c("MPZ","SOX10","NGFR","MBP"),
  Mast = c("TPSAB1","CPA3","TPSB2","KIT"),
  RBC = c("HBB","HBA1","ALAS2")
)

# ---- 每个簇的 marker 平均表达打分 ----
DefaultAssay(obj) <- "RNA"
mtx <- as.matrix(GetAssayData(obj, layer="data"))
score <- sapply(names(markers), function(ct) {
  g <- intersect(markers[[ct]], rownames(mtx))
  if (length(g) == 0) return(rep(0, length(levels(Idents(obj)))))
  colMeans(as.matrix(mtx[g, , drop=FALSE]))
})
cl.scores <- t(sapply(levels(Idents(obj)), function(cl) {
  cells <- WhichCells(obj, idents=cl)
  colMeans(score[cells, , drop=FALSE])
}))
cl.scores <- data.frame(cluster=rownames(cl.scores), cl.scores, check.names=FALSE)

# 每个簇 top 类型（按得分归一化排序取最高者）
assign_cl <- function(row) {
  s <- row[-1]
  s[is.na(s)] <- 0
  names(s)[which.max(s)]
}
cl.scores$assign <- apply(cl.scores, 1, assign_cl)
cat("== 聚类注释（按平均 marker 得分）==\n")
print(cl.scores[, c("cluster","assign")], row.names=FALSE)
write.csv(cl.scores, file.path(outdir, "cluster_celltype_scores.csv"), row.names=FALSE)

# ---- Ro/e 计算 ----
# Ro/e = 观察到的细胞类型 marker 在该簇 DEG 中的比例 / 该细胞类型 marker 占全部基因的比例
# DEG 来源：cluster_markers.csv（step2b）
mk <- read.csv(file.path(outdir, "cluster_markers.csv"))
if (nrow(mk) > 0) {
  allgenes <- rownames(obj)
  roe <- sapply(levels(Idents(obj)), function(cl) {
    degs <- mk$gene[mk$cluster == cl]
    sapply(names(markers), function(ct) {
      g <- intersect(markers[[ct]], allgenes)
      obs <- mean(degs %in% g)          # 观察比例
      exp_ <- length(g) / length(allgenes)  # 期望比例
      if (exp_ == 0) return(0)
      obs / exp_
    })
  })
  roe <- data.frame(cluster=colnames(roe), t(roe), check.names=FALSE)
  write.csv(roe, file.path(outdir, "cluster_roe.csv"), row.names=FALSE)
  cat("== Ro/e 矩阵（前 8 簇 × 前 6 型）==\n")
  print(head(roe[, 1:7], 8))
} else {
  cat("! cluster_markers.csv 为空，跳过 Ro/e\n")
}

# 保存注释 - 用 cell barcode 作为名称
cl2ct <- setNames(cl.scores$assign, cl.scores$cluster)
cell_barcodes <- colnames(obj)
cell_idents <- as.character(Idents(obj))
celltype_vec <- cl2ct[cell_idents]
names(celltype_vec) <- cell_barcodes
obj$celltype <- celltype_vec
saveRDS(obj, file.path(outdir, "GSE298501_integrated.rds"))
cat("celltype table:\n"); print(table(obj$celltype))
cat("== step2c done", format(Sys.time()), "==\n")
sink()
