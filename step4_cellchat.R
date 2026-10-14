# Step 4: CellChat - Cell-cell communication analysis
# Subsample per celltype to keep memory manageable
# Keep this file pure ASCII.

Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp",
           TEMP="D:/claude-study/Rtmp",
           TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/GSE298501"
sink(file.path(outdir, "step4_cellchat.log"), append=FALSE)
cat("== step4 start", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# Load integrated object
cat("Loading Seurat object...\n")
obj <- readRDS(file.path(outdir, "GSE298501_integrated.rds"))
cat("Cells:", ncol(obj), "\n")
DefaultAssay(obj) <- "RNA"
flush.console()

# Subsample per celltype (max 2000 per type, min 50)
set.seed(42)
cell_meta <- obj@meta.data
ct_table <- table(cell_meta$celltype)
cat("Cell types:\n"); print(ct_table)

sampled_cells <- c()
for (ct in names(ct_table)) {
  cells_ct <- which(cell_meta$celltype == ct)
  n_take <- min(2000, length(cells_ct))
  n_take <- max(n_take, min(50, length(cells_ct)))
  sampled_cells <- c(sampled_cells, sample(cells_ct, n_take))
}
cat("Subsampled cells:", length(sampled_cells), "\n")

obj_cc <- subset(obj, cells=colnames(obj)[sampled_cells])
rm(obj); gc(verbose=FALSE)
cat("CellChat object cells:", ncol(obj_cc), "\n")
flush.console()

# Prepare data
expr_mat <- as.matrix(GetAssayData(obj_cc, layer="data"))
meta <- obj_cc@meta.data
rm(obj_cc); gc(verbose=FALSE)

# Create CellChat object
cat("Creating CellChat object...\n")
cellchat <- createCellChat(object=expr_mat, meta=meta, group.by="celltype")
rm(expr_mat); gc(verbose=FALSE)

CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# Identify overexpressed genes and interactions
cat("Identifying overexpressed genes...\n")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
flush.console()

# Compute communication probability
cat("Computing communication probability...\n")
cellchat <- computeCommunProb(cellchat, population.size=TRUE)

# Filter low-communication (default params only)
cellchat <- filterCommunication(cellchat)

# Compute pathway-level communication
cat("Computing pathway-level communication...\n")
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate network
cellchat <- aggregateNet(cellchat)

# Compute centrality
cat("Computing network centrality...\n")
cellchat <- netAnalysis_computeCentrality(cellchat)
flush.console()

# Save CellChat object
saveRDS(cellchat, file.path(outdir, "cellchat_object.rds"))
cat("CellChat object saved\n")

# Extract ligand-receptor pairs
tryCatch({
  lr_pairs <- extractEnrichedLR(cellchat, signaling=names(cellchat@netP$pathways),
                                 geneLR.returned=FALSE)
  write.csv(lr_pairs, file.path(outdir, "lr_pairs.csv"), row.names=FALSE)
  cat("LR pairs:", nrow(lr_pairs), "\n")
}, error=function(e) {
  cat("LR extraction error:", e$message, "\n")
})

# Visualization
cat("Generating plots...\n")
flush.console()

# Communication heatmap
pdf(file.path(outdir, "communication_heatmap.pdf"), width=10, height=10)
tryCatch({
  netVisual_heatmap(cellchat, title.name="Num. interactions")
}, error=function(e) cat("Heatmap1 error:", e$message, "\n"))
tryCatch({
  netVisual_heatmap(cellchat, measure="weight", title.name="Interaction strength")
}, error=function(e) cat("Heatmap2 error:", e$message, "\n"))
dev.off()

# Centrality scores
tryCatch({
  centrality <- cellchat@netAnalysis$centrality
  write.csv(as.data.frame(centrality), file.path(outdir, "centrality_scores.csv"))
  cat("Centrality scores saved\n")
}, error=function(e) cat("Centrality error:", e$message, "\n"))

# Summary statistics
cat("\nSummary:\n")
cat("  Cell types:", length(unique(meta$celltype)), "\n")
cat("  Pathways:", length(cellchat@netP$pathways), "\n")
cat("  Total interactions:", sum(cellchat@net$count), "\n")

cat("== step4 done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
