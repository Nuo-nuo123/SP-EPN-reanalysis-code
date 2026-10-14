# GSE243240 noob reprocessing keeping 20,000 random autosomal probes (for the
# random-probe null models). Chunked like gse243240_noob_reprocess.R.
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages({library(minfi); library(data.table)})
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
idat_dir <- "D:/claude-study/Spinal Ependymoma/GEO_data/GSE243240/IDAT"
sink(file.path(outdir, "gse243240_random20k.log"), append=FALSE)
cat("== GSE243240 random-20k noob reprocessing", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# probe pool: the deposited 404,256 universe minus the 1,722 significant DMPs,
# fixed-seed sample of 20,000
dmp <- read.csv(file.path(outdir, "dmp_significant_annotated.csv"), stringsAsFactors=FALSE)
beta72 <- fread(file.path(outdir, "beta_72samples.csv"))
universe <- beta72[[1]]
pool <- setdiff(universe, dmp$cpg)
set.seed(20260919)
keep_probes <- sample(pool, 20000)
cat("pool:", length(pool), "| kept:", length(keep_probes), "\n")

smap <- read.csv(file.path(outdir, "gse243240_287_map.csv"), stringsAsFactors=FALSE)
gsm_450k <- smap$gsm[smap$platform == "GPL13534"]
gsm_epic <- smap$gsm[smap$platform == "GPL21145"]

process_platform <- function(gsms, plat_label, out_rds) {
  cat("\n--- noob:", plat_label, "---\n")
  all_files <- list.files(idat_dir, pattern="_Grn.idat$", full.names=TRUE)
  base_names <- sub("_Grn\\.idat$", "", all_files)
  gsm_of_base <- sub(".*/(GSM[0-9]+)_.*", "\\1", base_names)
  bases <- base_names[gsm_of_base %in% gsms]
  cat("Found", length(bases), "of", length(gsms), "samples\n")
  beta_list <- list()
  chunk <- 40
  for (start in seq(1, length(bases), by=chunk)) {
    end <- min(start+chunk-1, length(bases))
    tb <- data.frame(Basename = bases[start:end])
    cat("  Reading", start, "-", end, "...\n")
    rg <- read.metharray.exp(base=NULL, targets=tb, verbose=FALSE, extended=FALSE, force=TRUE)
    cat("  detectionP + preprocessNoob ...\n")
    detP <- detectionP(rg)
    mset <- preprocessNoob(rg)
    beta <- getBeta(mset)
    rm(rg, mset); gc()
    beta[detP > 0.01] <- NA
    rm(detP); gc()
    avail <- intersect(keep_probes, rownames(beta))
    beta_list[[length(beta_list)+1]] <- beta[avail, , drop=FALSE]
    cat("  Chunk:", nrow(beta_list[[length(beta_list)]]), "probes,",
        round(100*mean(is.na(beta_list[[length(beta_list)]])), 2), "% NA\n")
  }
  beta_all <- do.call(cbind, beta_list)
  cat(plat_label, "final:", nrow(beta_all), "x", ncol(beta_all), "\n")
  saveRDS(beta_all, file.path(outdir, out_rds))
  beta_all
}
b450n <- process_platform(gsm_450k, "450K", "gse243240_beta_random20k_450k.rds")
bepicn <- process_platform(gsm_epic, "EPIC", "gse243240_beta_random20k_epic.rds")
vcom <- intersect(rownames(b450n), rownames(bepicn))
Vall <- cbind(b450n[vcom, , drop=FALSE], bepicn[vcom, , drop=FALSE])
saveRDS(Vall, file.path(outdir, "gse243240_beta_random20k_noob.rds"))
cat("\nCombined random-20k noob:", nrow(Vall), "probes x", ncol(Vall), "samples\n")
cat("== done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
