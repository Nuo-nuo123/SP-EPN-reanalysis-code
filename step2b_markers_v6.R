# step2b v6: 50K cell subsample to fit in 28GB RAM
# Keep this file pure ASCII.

Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp",
           TEMP="D:/claude-study/Rtmp",
           TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/GSE298501"
logf <- file.path(outdir, "step2b_markers_v6.log")
sink(logf, split=TRUE)
cat("== step2b v6 start", format(Sys.time()), "==\n")

# Load merged_raw
cat("Loading merged_raw...\n")
flush.console()
obj <- readRDS(file.path(outdir, "GSE298501_merged_raw.rds"))
cat("Original cells:", ncol(obj), "\n")
flush.console()

# Subsample to 50K cells (stratified by cluster would be ideal but we don't have clusters yet)
set.seed(42)
sub_cells <- sample(colnames(obj), 50000)
obj <- subset(obj, cells=sub_cells)
cat("Subsampled to:", ncol(obj), "cells\n")
flush.console()
gc(verbose=FALSE)

# Integration
cat("Normalizing...\n")
flush.console()
obj <- NormalizeData(obj, verbose=FALSE)

cat("Finding variable features...\n")
flush.console()
obj <- FindVariableFeatures(obj, nfeatures=2000, verbose=FALSE)

cat("Scaling...\n")
flush.console()
obj <- ScaleData(obj, verbose=FALSE)

cat("PCA...\n")
flush.console()
obj <- RunPCA(obj, npcs=50, verbose=FALSE)

cat("Harmony...\n")
flush.console()
obj <- RunHarmony(obj, group.by.vars="lib", max_iter=20, verbose=FALSE)

cat("UMAP...\n")
flush.console()
obj <- RunUMAP(obj, reduction="harmony", dims=1:30, verbose=FALSE)

cat("Neighbors...\n")
flush.console()
obj <- FindNeighbors(obj, reduction="harmony", dims=1:30, verbose=FALSE)

cat("Clusters...\n")
flush.console()
obj <- FindClusters(obj, resolution=0.8, verbose=FALSE)
cat("Clusters found:", length(levels(obj$seurat_clusters)), "\n")
flush.console()
gc(verbose=FALSE)

# Join layers
cat("Joining layers...\n")
flush.console()
L <- Layers(obj[["RNA"]])
if (length(L) > 3) {
  obj <- JoinLayers(obj)
  cat("JoinLayers done, layers:", paste(Layers(obj[["RNA"]]), collapse=","), "\n")
}
flush.console()
gc(verbose=FALSE)

# Find markers
Idents(obj) <- "seurat_clusters"
clusters <- levels(Idents(obj))
cat("Finding markers for", length(clusters), "clusters...\n")
flush.console()

allmk <- NULL
ckpt <- file.path(outdir, "cluster_markers_ckpt.csv")
for (cl in clusters) {
  cat("-- cluster", cl, format(Sys.time()), "--\n")
  flush.console()
  tryCatch({
    mk <- FindMarkers(obj, ident.1=cl, only.pos=TRUE,
                      min.pct=0.25, logfc.threshold=0.25,
                      test.use="wilcox", verbose=FALSE)
    mk$cluster <- cl
    mk$gene <- rownames(mk)
    allmk <- rbind(allmk, mk)
    write.csv(allmk, ckpt, row.names=FALSE)
    cat("   markers:", nrow(mk), "total:", nrow(allmk), "\n")
    flush.console()
  }, error=function(e) {
    cat("   ERROR:", e$message, "\n")
    flush.console()
  })
  gc(verbose=FALSE)
}

# Save final
write.csv(allmk, file.path(outdir, "cluster_markers.csv"), row.names=FALSE)
cat("Total markers:", nrow(allmk), "\n")
flush.console()

saveRDS(obj, file.path(outdir, "GSE298501_integrated.rds"))
cat("Object saved\n")

cat("== step2b v6 done", format(Sys.time()), "==\n")
sink()
