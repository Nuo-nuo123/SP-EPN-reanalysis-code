# Step 3c only: SCENIC (optimized for speed)
# Keep this file pure ASCII.

Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp",
           TEMP="D:/claude-study/Rtmp",
           TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(Seurat)
  library(GENIE3)
  library(RcisTarget)
  library(AUCell)
  library(Matrix)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/GSE298501"
sink(file.path(outdir, "step3c_scenic.log"), append=FALSE)
cat("== step3c SCENIC start", format(Sys.time()), "==\n")
flush.console()

# Load object
cat("Loading Seurat object...\n")
obj <- readRDS(file.path(outdir, "GSE298501_integrated.rds"))
cat("Cells:", ncol(obj), "\n")
flush.console()

# Subsample 5K cells
set.seed(42)
sc_cells <- sample(colnames(obj), 5000)
obj_sc <- subset(obj, cells=sc_cells)
expr_sc <- GetAssayData(obj_sc, layer="data")
cat("SCENIC cells:", ncol(expr_sc), "\n")
flush.console()

# TF list
library(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
data(hg19_motifAnnotation_cisbpOnly)
tf_all <- unique(hg19_motifAnnotation_cisbpOnly$TF)

# Filter expressed genes
genes_expressed <- rownames(expr_sc)[rowSums(expr_sc > 0) >= 3]
tf_expressed <- intersect(tf_all, genes_expressed)
cat("Expressed TFs:", length(tf_expressed), "\n")
flush.console()

# Select top variable genes as targets (much faster than all genes)
cat("Selecting top variable genes as targets...\n")
expr_dense <- as.matrix(expr_sc)
gene_var <- apply(expr_dense, 1, var)
top_var <- names(sort(gene_var, decreasing=TRUE))[1:min(2000, sum(gene_var > 0))]
# Make sure TFs are included as potential targets
target_genes <- unique(c(top_var, tf_expressed))
target_genes <- intersect(target_genes, rownames(expr_dense))
cat("Target genes:", length(target_genes), "\n")
flush.console()

# Subset expression to target genes only
expr_target <- expr_dense[target_genes, ]
rm(expr_dense); gc(verbose=FALSE)

# GENIE3 with reduced trees
cat("Running GENIE3 (100 trees,", length(target_genes), "targets)...\n")
flush.console()
gene_net <- GENIE3(expr_target, regulators=tf_expressed, nTrees=100, K="sqrt", verbose=TRUE)

# Extract links
links <- getLinkList(gene_net, threshold=0.001)
cat("GENIE3 links:", nrow(links), "\n")
flush.console()

# Build regulons
regulons_list <- split(as.character(links$targetGene), links$regulatoryGene)
regulons_list <- lapply(regulons_list, function(x) unique(x))
cat("Regulons:", length(regulons_list), "\n")
flush.console()

# cisTarget
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

saveRDS(list(regulons=regulons_list, motif_enrichment=motifEnrichment,
             motif_enrichment_filtered=motifEnrichment_filtered),
        file.path(outdir, "scenic_regulons.rds"))

# AUCell
if (nrow(motifEnrichment_filtered) > 0) {
  cat("Running AUCell...\n")
  flush.console()
  cells_rankings <- AUCell_buildRankings(expr_sc, nCores=1, plotStats=FALSE)

  # Try multiple ways to match TFs from cisTarget to regulons_list
  cat("Matching regulons for AUCell...\n")
  cat("  motifEnrichment cols:", paste(colnames(motifEnrichment_filtered), collapse=", "), "\n")
  cat("  regulons_list names (first 5):", paste(head(names(regulons_list), 5), collapse=", "), "\n")

  # Method 1: Use TF column if it exists
  if ("TF" %in% colnames(motifEnrichment_filtered)) {
    filtered_tfs <- unique(motifEnrichment_filtered$TF)
    cat("  TFs from cisTarget:", length(filtered_tfs), "\n")
    regulons_filtered <- regulons_list[intersect(filtered_tfs, names(regulons_list))]
  } else {
    # Method 2: Use row names
    rn <- rownames(motifEnrichment_filtered)
    regulons_filtered <- regulons_list[intersect(rn, names(regulons_list))]
  }

  # Remove empty/null entries
  regulons_filtered <- regulons_filtered[sapply(regulons_filtered, function(x) length(x) > 0)]
  cat("  Regulons for AUCell:", length(regulons_filtered), "\n")
  flush.console()

  if (length(regulons_filtered) == 0) {
    cat("WARNING: No regulons matched. Using all regulons for AUCell.\n")
    regulons_filtered <- regulons_list
  }

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
cat("== step3c done", format(Sys.time()), "==\n")
sink()
