# Step 8: SHAP Explainability Analysis
# Input: stacking_model.rds + training data (from step6)
# Output: shap_values.csv, shap_bar.pdf, shap_beeswarm.pdf, shap_waterfall.pdf
# Keep this file pure ASCII.

Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp",
           TEMP="D:/claude-study/Rtmp",
           TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(fastshap)
  library(randomForest)
  library(data.table)
  library(ggplot2)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "step8_shap.log"), append=TRUE)
cat("== step8 start", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# Load model and data
model_info <- readRDS(file.path(outdir, "stacking_model.rds"))
feature_genes <- model_info$feature_genes

expr_raw <- readRDS(file.path(outdir, "GSE242994_expr.rds"))
if (is.list(expr_raw) && "expr" %in% names(expr_raw)) { expr <- expr_raw$expr } else { expr <- expr_raw }
labels <- fread(file.path(outdir, "GSE242994_rna_labels.csv"))
common <- intersect(colnames(expr), labels$sample)

X <- t(as.matrix(expr[feature_genes, common]))
y <- as.factor(labels$Transcriptional_subtype[match(common, labels$sample)])

# Train a Random Forest for SHAP (tree-based models work best with SHAP)
cat("Training RF for SHAP...\n")
rf_model <- randomForest(X, y, ntree=500, importance=TRUE)

# ============================================================
# Compute SHAP values
# ============================================================
cat("Computing SHAP values...\n")

# Prediction function for fastshap
pred_wrap <- function(object, newdata) {
  predict(object, newdata, type="prob")[,2]
}

shap_vals <- explain(rf_model,
                     X=X,
                     pred_wrapper=pred_wrap,
                     nsim=50,
                     shap_only=FALSE)

# Save SHAP values
shap_raw <- unclass(shap_vals)
shap_df <- as.data.frame(shap_raw)
shap_df$sample <- common
write.csv(shap_df, file.path(outdir, "shap_values.csv"), row.names=FALSE)
cat("SHAP values saved:", nrow(shap_df), "samples x", ncol(shap_df)-1, "genes\n")

# ============================================================
# Visualizations
# ============================================================

# 1. Bar plot (global importance)
cat("Generating bar plot...\n")
pdf(file.path(outdir, "shap_bar.pdf"), width=10, height=8)
shap_mat <- as.matrix(shap_raw)
shap_long <- reshape2::melt(shap_mat)
mean_abs_shap <- aggregate(value ~ Var2, data=shap_long, FUN=function(x) mean(abs(x)))
colnames(mean_abs_shap) <- c("gene", "mean_abs_shap")
mean_abs_shap <- mean_abs_shap[order(mean_abs_shap$mean_abs_shap, decreasing=TRUE), ]
top_genes <- head(mean_abs_shap, 20)

p <- ggplot(top_genes, aes(x=reorder(gene, mean_abs_shap), y=mean_abs_shap)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  labs(title="Top 20 Features by Mean |SHAP value|",
       x="Gene", y="Mean |SHAP value|") +
  theme_minimal()
print(p)
dev.off()

# 2. Beeswarm plot
cat("Generating beeswarm plot...\n")
pdf(file.path(outdir, "shap_beeswarm.pdf"), width=10, height=10)

# Prepare data for beeswarm
shap_long$gene <- factor(shap_long$Var2,
                          levels=rev(top_genes$gene))
shap_long <- shap_long[shap_long$Var2 %in% top_genes$gene, ]

# Add feature values
X_df <- as.data.frame(X[, top_genes$gene])
X_long <- reshape2::melt(X_df)
X_long$gene <- factor(X_long$Var2, levels=rev(top_genes$gene))
shap_long$feature_value <- X_long$value

p <- ggplot(shap_long, aes(x=gene, y=value, color=feature_value)) +
  geom_jitter(width=0.2, size=0.8, alpha=0.6) +
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=median(shap_long$feature_value)) +
  coord_flip() +
  labs(title="SHAP Beeswarm Plot",
       x="Gene", y="SHAP value", color="Feature\nvalue") +
  theme_minimal() +
  geom_hline(yintercept=0, linetype="dashed", color="gray50")
print(p)
dev.off()

# 3. Waterfall plot (single sample)
cat("Generating waterfall plot...\n")
pdf(file.path(outdir, "shap_waterfall.pdf"), width=12, height=8)

# Pick first sample
sample_idx <- 1
shap_sample <- shap_mat[sample_idx, ]
base_value <- mean(predict(rf_model, X, type="prob")[,2])

# Sort by absolute SHAP
shap_sorted <- sort(shap_sample, decreasing=TRUE)
top10 <- head(names(shap_sorted), 10)

waterfall_data <- data.frame(
  gene=top10,
  shap=shap_sorted[top10],
  cumsum=cumsum(shap_sorted[top10])
)

p <- ggplot(waterfall_data, aes(x=reorder(gene, cumsum), y=shap)) +
  geom_bar(stat="identity", fill=ifelse(waterfall_data$shap > 0, "red", "blue")) +
  coord_flip() +
  labs(title=paste("Waterfall Plot - Sample:", common[sample_idx]),
       x="Gene", y="SHAP value") +
  theme_minimal() +
  geom_hline(yintercept=0, linetype="dashed")
print(p)
dev.off()

cat("\n== step8 done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
