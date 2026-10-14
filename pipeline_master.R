# Master pipeline: runs steps 4-9 sequentially
# Each step is run via source() with error handling
# Keep this file pure ASCII.

Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp",
           TEMP="D:/claude-study/Rtmp",
           TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

src_dir <- "D:/claude-study/Spinal Ependymoma/R_pkg_logs/src"
log_dir <- "D:/claude-study/Spinal Ependymoma/analysis"
pipeline_log <- file.path(log_dir, "pipeline_master.log")

log_msg <- function(...) {
  msg <- paste0(..., "\n")
  cat(msg)
  cat(msg, file=pipeline_log, append=TRUE)
  flush.console()
}

log_msg("== MASTER PIPELINE start ", Sys.time(), " ==")

steps <- list(
  list(name="step4_cellchat",       script=file.path(src_dir, "step4_cellchat.R")),
  list(name="step5_dmr_integration", script=file.path(src_dir, "step5_dmr_integration.R")),
  list(name="step6_ml_ensemble",    script=file.path(src_dir, "step6_ml_ensemble.R")),
  list(name="step7_validation",     script=file.path(src_dir, "step7_validation.R")),
  list(name="step8_shap",           script=file.path(src_dir, "step8_shap.R")),
  list(name="step9_nomogram_shiny", script=file.path(src_dir, "step9_nomogram_shiny.R"))
)

results <- list()
for (step in steps) {
  log_msg("========================================")
  log_msg("RUNNING: ", step$name, " ", Sys.time())
  log_msg("========================================")

  start_time <- Sys.time()
  tryCatch({
    source(step$script, local=TRUE)
    elapsed <- difftime(Sys.time(), start_time, units="mins")
    results[[step$name]] <- list(status="OK", time=as.numeric(elapsed))
    log_msg(">> ", step$name, " DONE in ", round(as.numeric(elapsed), 1), " min")
  }, error=function(e) {
    elapsed <- difftime(Sys.time(), start_time, units="mins")
    results[[step$name]] <<- list(status="ERROR", error=e$message, time=as.numeric(elapsed))
    log_msg(">> ", step$name, " ERROR: ", e$message)
    log_msg(">> Continuing to next step...")
  })
  gc(verbose=FALSE)
}

log_msg("")
log_msg("========================================")
log_msg("PIPELINE SUMMARY")
log_msg("========================================")
for (nm in names(results)) {
  r <- results[[nm]]
  log_msg(sprintf("%-30s %s (%.1f min)", nm, r$status, r$time))
}

log_msg("")
log_msg("== MASTER PIPELINE done ", Sys.time(), " ==")
