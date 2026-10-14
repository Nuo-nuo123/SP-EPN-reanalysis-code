# Fig 1A: UMAP colored by cell type
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

figdir <- "D:/claude-study/Spinal Ependymoma/manuscript/figures"
cat("Loading Seurat object...\n")
obj <- readRDS("D:/claude-study/Spinal Ependymoma/analysis/GSE298501/GSE298501_integrated.rds")
cat("Object loaded:", ncol(obj), "cells\n")

# celltype colors
ct_colors <- c(
  Astrocyte="#1f77b4", `Ependymal-tumor`="#d62728", Microglia="#2ca02c",
  Endothelial="#9467bd", T_cell="#ff7f0e", Pericyte="#8c564b",
  Fibroblast="#e377c2", Oligodendrocyte="#17becf", OPC="#bcbd22", Mast="#7f7f7f"
)

pdf(file.path(figdir, "Fig1A_umap.pdf"), width=8, height=7)
p <- DimPlot(obj, reduction="umap", group.by="celltype", cols=ct_colors,
             pt.size=0.3, label=FALSE) +
  ggtitle("50,000 nuclei from 14 spinal ependymomas") +
  theme(legend.position="right") +
  guides(color=guide_legend(override.aes=list(size=3), ncol=1))
print(p)
dev.off()
cat("Fig1A saved\n")
