# Step2: GSE298501 标准化 + Harmony 去批次 + UMAP + 聚类 + marker
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6", TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(), "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(harmony))

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/GSE298501"
logf <- file.path(outdir, "step2_integrate.log")
sink(logf, split=TRUE); cat("== step2 start", format(Sys.time()), "==\n")

obj <- readRDS(file.path(outdir, "GSE298501_merged_raw.rds"))
cat("cells:", ncol(obj), "\n")

obj <- NormalizeData(obj, normalization.method="LogNormalize", scale.factor=10000)
obj <- FindVariableFeatures(obj, selection.method="vst", nfeatures=2000)
obj <- ScaleData(obj)                     # 只对变量基因
obj <- RunPCA(obj, npcs=50, verbose=FALSE)
obj <- RunHarmony(obj, group.by.vars="lib", max_iter=20)
obj <- RunUMAP(obj, reduction="harmony", dims=1:30, verbose=FALSE)
obj <- FindNeighbors(obj, reduction="harmony", dims=1:30, verbose=FALSE)
obj <- FindClusters(obj, resolution=0.8, verbose=FALSE)
cat("clusters:", paste(levels(Idents(obj)), collapse=","), "\n")
saveRDS(obj, file.path(outdir, "GSE298501_integrated.rds"))

# 聚类 marker（全量 wilcox 较慢，取 min.pct=0.25）
mk <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, verbose=FALSE)
write.csv(mk, file.path(outdir, "cluster_markers.csv"), row.names=FALSE)
cat("markers saved, top:\n")
print(head(mk, 10))
cat("== step2 done", format(Sys.time()), "==\n")
sink()
