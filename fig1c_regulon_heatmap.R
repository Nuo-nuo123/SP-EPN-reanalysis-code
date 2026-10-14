# Fig1C: Regulon activity heatmap by cell type (from SCENIC regulon matrix)
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/GSE298501"
figdir <- "D:/claude-study/Spinal Ependymoma/manuscript/figures"
cat("Loading regulon activity matrix (512MB)...\n")
ram <- fread(file.path(outdir, "regulon_activity_matrix.csv"))
cat("Matrix:", nrow(ram), "x", ncol(ram), "\n")

# rows = regulons, first col = regulon name
reg_names <- ram[[1]]
ram[[1]] <- NULL
ram_mat <- as.matrix(ram)
rownames(ram_mat) <- reg_names
rm(ram); gc()

# column names are cell barcodes with library prefix like GSM9015499_MPE_A_1_...
# sample 5000 cells for speed
set.seed(42)
if (ncol(ram_mat) > 5000) {
  idx <- sample(ncol(ram_mat), 5000)
  ram_sub <- ram_mat[, idx]
} else {
  ram_sub <- ram_mat
}
cell_bc <- colnames(ram_sub)
lib_prefix <- sub("_.*", "", cell_bc)  # not informative; group by entity

# Load cell type annotation from Seurat object metadata (subsample aligned?)
# Simpler: derive entity (MPE/SPEPN) from barcode, and use celltype from integrated object
# The regulon matrix columns = 50K cells of integrated object order.
# Load celltype from the Seurat object metadata instead (fast via readRDS of just metadata)
cat("Loading Seurat object for celltype metadata...\n")
obj <- readRDS(file.path(outdir, "GSE298501_integrated.rds"))
# match barcodes
mt <- obj@meta.data
common_cells <- intersect(colnames(ram_sub), rownames(mt))
cat("Matched cells:", length(common_cells), "\n")
ct <- mt[common_cells, "celltype"]
ram_ct <- ram_sub[, common_cells]

# Mean regulon activity per cell type
cts <- sort(unique(ct))
mean_act <- sapply(cts, function(c) rowMeans(ram_ct[, ct == c, drop=FALSE]))
colnames(mean_act) <- cts

# Select top 30 most variable regulons across cell types
rv <- apply(mean_act, 1, sd)
top30 <- names(sort(rv, decreasing=TRUE))[1:30]
mat30 <- mean_act[top30, ]
# z-score rows
mat30_z <- t(apply(mat30, 1, scale))
colnames(mat30_z) <- cts

pdf(file.path(figdir, "Fig1C_regulon_heatmap.pdf"), width=9, height=9)
pheatmap(mat30_z, cluster_rows=TRUE, cluster_cols=TRUE,
         main="Top 30 regulons: mean activity by cell type (z-score)",
         fontsize=9, color=colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()
cat("Fig1C saved\n")

# also save the mean activity matrix
write.csv(mean_act, file.path(outdir, "regulon_mean_by_celltype.csv"))
cat("Mean activity table saved\n")
