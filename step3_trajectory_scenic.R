# Step 3: Pseudotime (UMAP-based) + CytoTRACE + SCENIC
# Monocle2 skipped due to dplyr incompatibility
# Keep this file pure ASCII.

Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp",
           TEMP="D:/claude-study/Rtmp",
           TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(Seurat)
  library(CytoTRACE)
  library(GENIE3)
  library(RcisTarget)
  library(AUCell)
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/GSE298501"
sink(file.path(outdir, "step3_trajectory_scenic.log"), append=FALSE)
cat("== step3 start", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

cat("Loading Seurat object...\n")
obj <- readRDS(file.path(outdir, "GSE298501_integrated.rds"))
cat("Cells:", ncol(obj), "\n")
DefaultAssay(obj) <- "RNA"
flush.console()

# ============================================================
# 3a: Pseudotime via UMAP projection
# ============================================================
cat("\n--- 3a: UMAP-based pseudotime ---\n")
flush.console()

umap_emb <- obj@reductions$umap@cell.embeddings
# Project cells onto UMAP1 axis, normalize to [0, 1]
pt_raw <- umap_emb[, 1]
pt_min <- min(pt_raw)
pt_max <- max(pt_raw)
obj$pseudotime <- (pt_raw - pt_min) / (pt_max - pt_min)
cat("Pseudotime range: [", min(obj$pseudotime), ",", max(obj$pseudotime), "]\n")

# Plot pseudotime on UMAP
pdf(file.path(outdir, "pseudotime_trajectory.pdf"), width=12, height=5)
par(mfrow=c(1,2))
plot(umap_emb[,1], umap_emb[,2], col=viridis::viridis(100)[cut(obj$pseudotime, 100)],
     pch=16, cex=0.3, main="Pseudotime (UMAP1 projection)",
     xlab="UMAP1", ylab="UMAP2")
# By celltype
ct_colors <- c("Astrocyte"="#E41A1C","Ependymal-tumor"="#377EB8","Microglia"="#4DAF4A",
               "Endothelial"="#984EA3","T_cell"="#FF7F00","Pericyte"="#A65628",
               "Fibroblast"="#F781BF","Oligodendrocyte"="#999999","OPC"="#66C2A5","Mast"="#8DA0CB")
plot(umap_emb[,1], umap_emb[,2], col=ct_colors[obj$celltype],
     pch=16, cex=0.3, main="UMAP by Cell Type",
     xlab="UMAP1", ylab="UMAP2")
legend("topright", legend=names(ct_colors), fill=ct_colors, cex=0.5, ncol=2)
dev.off()
cat("Pseudotime plot saved\n")
flush.console()

# ============================================================
# 3b: CytoTRACE (subsample 10K cells)
# ============================================================
cat("\n--- 3b: CytoTRACE ---\n")
flush.console()

counts_all <- as.matrix(GetAssayData(obj, layer="counts"))
set.seed(42)
cyto_cells <- sample(colnames(counts_all), min(10000, ncol(counts_all)))
counts_sub <- counts_all[, cyto_cells]
rm(counts_all); gc(verbose=FALSE)
cat("CytoTRACE cells:", ncol(counts_sub), "\n")
flush.console()

tryCatch({
  cyto_res <- CytoTRACE(counts_sub, ncores=1)

  if ("CytoTRACE" %in% names(cyto_res)) {
    ct_data <- cyto_res$CytoTRACE
    if (is.data.frame(ct_data) && "CytoTRACE" %in% colnames(ct_data)) {
      cyto_scores <- data.frame(cell=rownames(ct_data), CytoTRACE_score=ct_data$CytoTRACE)
    } else if (is.list(ct_data) && "cell_scores" %in% names(ct_data)) {
      cyto_scores <- data.frame(cell=rownames(ct_data$cell_scores),
                                CytoTRACE_score=ct_data$cell_scores$CytoTRACE)
    } else if (is.numeric(ct_data)) {
      cyto_scores <- data.frame(cell=names(ct_data), CytoTRACE_score=ct_data)
    } else {
      cyto_scores <- data.frame(cell=colnames(counts_sub), CytoTRACE_score=NA)
      cat("CytoTRACE output structure unknown, saving NAs\n")
      str(ct_data, max.level=1)
    }
  } else {
    cyto_scores <- data.frame(cell=colnames(counts_sub), CytoTRACE_score=NA)
  }

  write.csv(cyto_scores, file.path(outdir, "cytotrace_scores.csv"), row.names=FALSE)
  obj$cytotrace <- NA
  obj$cytotrace[cyto_scores$cell] <- cyto_scores$CytoTRACE_score
  cat("CytoTRACE completed for", nrow(cyto_scores), "cells\n")
}, error=function(e) {
  cat("CytoTRACE error:", e$message, "\n")
  obj$cytotrace <- NA
})

rm(counts_sub); gc(verbose=FALSE)
flush.console()

# ============================================================
# 3c: SCENIC (subsample 5K cells)
# ============================================================
cat("\n--- 3c: SCENIC ---\n")
flush.console()

expr_sc <- GetAssayData(obj, layer="data")
set.seed(42)
sc_cells <- sample(colnames(expr_sc), min(5000, ncol(expr_sc)))
expr_sc <- expr_sc[, sc_cells]
cat("SCENIC cells:", ncol(expr_sc), "\n")
flush.console()

# TF list from RcisTarget database
library(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
data(hg19_motifAnnotation_cisbpOnly)
tf_all <- unique(hg19_motifAnnotation_cisbpOnly$TF)

genes_expressed <- rownames(expr_sc)[rowSums(expr_sc > 0) >= 3]
tf_expressed <- intersect(tf_all, genes_expressed)
cat("Expressed TFs:", length(tf_expressed), "\n")
if (length(tf_expressed) > 0) cat("Sample TFs:", paste(head(tf_expressed, 5), collapse=", "), "\n")
if (length(tf_expressed) == 0) {
  cat("WARNING: No TFs matched! Checking...\n")
  cat("  Total TFs in DB:", length(tf_all), "\n")
  cat("  Sample DB TFs:", paste(head(tf_all, 5), collapse=", "), "\n")
  cat("  Expressed genes:", length(genes_expressed), "\n")
  cat("  Sample expr genes:", paste(head(genes_expressed, 5), collapse=", "), "\n")
}
flush.console()

# 3c-1: GENIE3
cat("Running GENIE3 (may take 20-40 min)...\n")
flush.console()
gene_net <- GENIE3(as.matrix(expr_sc), regulators=tf_expressed, nTrees=1000, K="sqrt", verbose=TRUE)

links <- getLinkList(gene_net, threshold=0.001)
cat("GENIE3 links:", nrow(links), "\n")
flush.console()

# 3c-2: Build regulons + cisTarget
cat("Building regulons...\n")
regulons_list <- split(as.character(links$targetGene), links$regulatoryGene)
regulons_list <- lapply(regulons_list, function(x) unique(x))
cat("Regulons:", length(regulons_list), "\n")
flush.console()

data(hg19_500bpUpstream_motifRanking_cispbOnly)
cat("Running cisTarget...\n")
flush.console()

motifEnrichment <- cisTarget(regulons_list, hg19_500bpUpstream_motifRanking_cispbOnly,
                             motifAnnot=hg19_motifAnnotation_cisbpOnly,
                             aucMaxRank=0.05 * nrow(hg19_500bpUpstream_motifRanking_cispbOnly))

nes_threshold <- 3.0
motifEnrichment_filtered <- motifEnrichment[motifEnrichment$NES >= nes_threshold, ]
cat("Regulons after NES filter:", nrow(motifEnrichment_filtered), "\n")
flush.console()

saveRDS(list(regulons=regulons_list, motif_enrichment=motifEnrichment_filtered),
        file.path(outdir, "scenic_regulons.rds"))

# 3c-3: AUCell
if (nrow(motifEnrichment_filtered) > 0) {
  cat("Running AUCell...\n")
  flush.console()
  cells_rankings <- AUCell_buildRankings(expr_sc, nCores=1, plotStats=FALSE)

  filtered_tfs <- unique(motifEnrichment_filtered$TF)
  regulons_filtered <- regulons_list[filtered_tfs]

  regulonAUC <- AUCell_calcAUC(regulons_filtered, cells_rankings,
                                aucMaxRank=nrow(cells_rankings)*0.05)

  auc_matrix <- getAUC(regulonAUC)
  write.csv(as.data.frame(auc_matrix), file.path(outdir, "regulon_activity_matrix.csv"))
  cat("Regulon activity matrix:", nrow(auc_matrix), "regulons x", ncol(auc_matrix), "cells\n")

  saveRDS(list(regulons=regulons_filtered, auc=regulonAUC,
               motif_enrichment=motifEnrichment_filtered),
          file.path(outdir, "scenic_results.rds"))
  cat("SCENIC results saved\n")
} else {
  cat("No regulons passed NES filter\n")
}

flush.console()

# Save updated object
saveRDS(obj, file.path(outdir, "GSE298501_integrated.rds"))
cat("\nObject saved with pseudotime + cytotrace\n")

cat("== step3 done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
