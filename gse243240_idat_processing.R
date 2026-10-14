# GSE243240: process IDATs -> beta for top DMPs (validation cohort)
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(minfi)
  library(data.table)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
idat_dir <- "D:/claude-study/Spinal Ependymoma/GEO_data/GSE243240/IDAT"
sink(file.path(outdir, "gse243240_idat.log"), append=FALSE)
cat("== GSE243240 IDAT processing start", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# ============================================================
# 1. Load our top 1000 DMPs (defined in 72-sample cohort)
# ============================================================
dmp <- read.csv(file.path(outdir, "dmp_results.csv"), stringsAsFactors=FALSE)
top_dmps <- dmp$cpg[order(dmp$pvalue)][1:1000]
cat("Top DMPs to extract:", length(top_dmps), "\n")

# ============================================================
# 2. Map samples to platforms
# ============================================================
smap <- read.csv(file.path(outdir, "gse243240_287_map.csv"), stringsAsFactors=FALSE)
gsm_450k <- smap$gsm[smap$platform == "GPL13534"]
gsm_epic <- smap$gsm[smap$platform == "GPL21145"]
cat("450K samples:", length(gsm_450k), "EPIC samples:", length(gsm_epic), "\n")

# ============================================================
# 3. Process 450K IDATs
# ============================================================
process_platform <- function(gsms, plat_label) {
  cat("\n--- Processing", plat_label, "---\n")
  targets <- data.frame(
    Basename = file.path(idat_dir, gsms)
  )
  targets$Basename <- sub("\\.idat$", "", targets$Basename)
  # find actual basenames
  all_files <- list.files(idat_dir, pattern="_Grn.idat$", full.names=TRUE)
  base_names <- sub("_Grn\\.idat$", "", all_files)
  gsm_of_base <- sub(".*/(GSM[0-9]+)_.*", "\\1", base_names)
  keep <- gsm_of_base %in% gsms
  bases <- base_names[keep]
  cat("Found", length(bases), "samples for", plat_label, "\n")
  if (length(bases) == 0) return(NULL)

  # read in chunks to limit memory
  beta_list <- list()
  chunk <- 50
  for (start in seq(1, length(bases), by=chunk)) {
    end <- min(start+chunk-1, length(bases))
    tb <- data.frame(Basename = bases[start:end])
    cat("  Reading samples", start, "-", end, "...\n")
    rg <- read.metharray.exp(base=NULL, targets=tb, verbose=FALSE, extended=FALSE, force=TRUE)
    cat("  Normalizing...\n")
    beta <- getBeta(rg)
    rm(rg); gc()
    # subset to top DMPs
    avail <- intersect(top_dmps, rownames(beta))
    beta_list[[length(beta_list)+1]] <- beta[avail, , drop=FALSE]
    cat("  Chunk beta:", nrow(beta_list[[length(beta_list)]]), "probes\n")
  }
  # combine
  beta_all <- do.call(cbind, beta_list)
  cat(plat_label, "final:", nrow(beta_all), "probes x", ncol(beta_all), "samples\n")
  return(beta_all)
}

f450 <- file.path(outdir, "gse243240_beta_450k.rds")
if (file.exists(f450)) {
  cat("450K already processed, loading...\n")
  b450 <- readRDS(f450)
} else {
  b450 <- process_platform(gsm_450k, "450K")
  saveRDS(b450, f450)
  cat("450K saved\n")
}

b_epic <- process_platform(gsm_epic, "EPIC")
saveRDS(b_epic, file.path(outdir, "gse243240_beta_epic.rds"))
cat("EPIC saved\n")

# ============================================================
# 4. Combine both platforms on common probes
# ============================================================
common_probes <- intersect(rownames(b450), rownames(b_epic))
cat("\nCommon probes:", length(common_probes), "\n")
beta_comb <- cbind(b450[common_probes, , drop=FALSE],
                   b_epic[common_probes, , drop=FALSE])
cat("Combined beta:", nrow(beta_comb), "x", ncol(beta_comb), "\n")
saveRDS(beta_comb, file.path(outdir, "gse243240_beta_dmp.rds"))
write.csv(beta_comb, file.path(outdir, "gse243240_beta_dmp.csv"), row.names=TRUE)
cat("Combined saved\n")

cat("\n== IDAT processing done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
