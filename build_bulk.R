# 构建 4 个 bulk 队列的表达矩阵（series matrix -> 基因水平 -> RDS）
# v2：注释表用 readLines 逐行解析（fread 会被 GO 超长字段行搞崩）
Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6", TMP="D:/claude-study/Rtmp", TEMP="D:/claude-study/Rtmp", TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(), "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))
library(data.table)
outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

parse_series <- function(gzpath) {
  txt <- readLines(gzfile(gzpath), n=20000)
  beg <- which(txt == "!series_matrix_table_begin")
  hdr <- strsplit(txt[beg+1], "\t")[[1]]
  hdr <- gsub('"', '', hdr)
  dat <- fread(cmd=paste("gzip -dc", shQuote(gzpath)), skip=beg, header=FALSE)
  colnames(dat) <- hdr
  expr <- as.matrix(dat[, -1, with=FALSE]); rownames(expr) <- dat[[1]]
  suppressWarnings(storage.mode(expr) <- "double")
  meta <- txt[1:(beg-1)]
  g <- grepl("^!Sample_", meta)
  m2 <- strsplit(gsub('"', '', meta[g]), "\t")
  snames <- sapply(m2, function(z) z[1])
  vals <- lapply(m2, function(z) z[-1])
  sm <- data.frame(sample=colnames(expr), stringsAsFactors=FALSE)
  for (i in seq_along(vals)) {
    if (length(vals[[i]]) == ncol(expr)) sm[[gsub("!Sample_", "", snames[i])]] <- vals[[i]]
  }
  list(expr=expr, samples=sm)
}

read_annot <- function(annotgz) {
  # 表头在第 28 行（^Annotation / !Annotation_* 前缀 27 行），只取 ID 与 Gene symbol 两列
  lines <- readLines(gzfile(annotgz))
  hdr <- strsplit(lines[28], "\t")[[1]]
  idi <- which(hdr == "ID"); gi <- which(hdr == "Gene symbol")
  if (!length(idi) || !length(gi)) stop("annot header cols not found")
  parts <- strsplit(lines[29:length(lines)], "\t", fixed=TRUE)
  nmax <- max(idi, gi)
  probe <- vapply(parts, function(p) if (length(p) >= nmax) p[[idi]] else NA_character_, "")
  gene  <- vapply(parts, function(p) if (length(p) >= nmax) p[[gi]] else NA_character_, "")
  df <- data.frame(probe=probe, gene=gene, stringsAsFactors=FALSE)
  df[!is.na(df$gene) & df$gene != "" & df$gene != "---" & !is.na(df$probe) & df$probe != "", ]
}

build <- function(gse, platform) {
  gz <- list.files(paste0("D:/claude-study/Spinal Ependymoma/GEO_data/", gse),
                   pattern="series_matrix[.]txt[.]gz$", full.names=TRUE)[1]
  cat("== ", gse, " from ", gz, "\n", sep="")
  ps <- parse_series(gz)
  cat("  probes:", nrow(ps$expr), "samples:", ncol(ps$expr), "\n")
  agz <- paste0("D:/claude-study/Spinal Ependymoma/GEO_data/platforms/", platform, ".annot.gz")
  if (!file.exists(agz)) { cat("  ! annot missing, saving probe-level only\n"); saveRDS(ps, file.path(outdir, paste0(gse, "_probe.rds"))); return(invisible(NULL)) }
  df <- read_annot(agz)
  df <- df[df$probe %in% rownames(ps$expr), ]
  cat("  annot mapped:", nrow(df), "probes\n")
  e <- ps$expr[df$probe, , drop=FALSE]
  df$m <- rowMeans(e, na.rm=TRUE)
  df <- df[order(df$gene, -df$m), ]
  df <- df[!duplicated(df$gene), ]
  expr <- e[df$probe, , drop=FALSE]
  rownames(expr) <- df$gene
  cat("  gene-level:", nrow(expr), "genes\n")
  saveRDS(list(expr=expr, samples=ps$samples, probes=df), file.path(outdir, paste0(gse, "_expr.rds")))
  write.csv(ps$samples, file.path(outdir, paste0(gse, "_samples.csv")), row.names=FALSE)
}

build("GSE125861", "GPL570")
build("GSE64415", "GPL570")
build("GSE66354", "GPL570")
build("GSE66787", "GPL11532")
cat("== bulk build done ==\n")
