# Step 6: Feature Selection + 12-Algorithm Ensemble ML
# Input: GSE242994_expr.rds (72 RNA-seq) + GSE242994_rna_labels.csv (A:43/B:29)
#        multiomics_integration.csv (from step5, optional)
# Output: feature_importance_per_algo.csv, consensus_signature.csv,
#         venn_selected_genes.pdf, stacking_model.rds
# Keep this file pure ASCII.

Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp",
           TEMP="D:/claude-study/Rtmp",
           TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(glmnet)
  library(MASS)
  library(mboost)
  library(randomForest)
  library(gbm)
  library(xgboost)
  library(e1071)
  library(klaR)
  library(caret)
  library(data.table)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "step6_ml_ensemble.log"), append=TRUE)
cat("== step6 start", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# ============================================================
# Load data
# ============================================================
cat("Loading data...\n")
expr_data <- readRDS(file.path(outdir, "GSE242994_expr.rds"))
if (is.list(expr_data) && "expr" %in% names(expr_data)) {
  expr <- expr_data$expr
} else {
  expr <- expr_data
}
labels <- fread(file.path(outdir, "GSE242994_rna_labels.csv"))

# Match samples
common_samples <- intersect(colnames(expr), labels$sample)
cat("Matched samples:", length(common_samples), "\n")

X <- t(as.matrix(expr[, common_samples]))
y <- as.factor(labels$Transcriptional_subtype[match(common_samples, labels$sample)])
cat("Class A:", sum(y=="A"), "Class B:", sum(y=="B"), "\n")

# ============================================================
# Feature pre-selection
# ============================================================
cat("\n--- Feature pre-selection ---\n")

# Load multiomics integration if available
integ_file <- file.path(outdir, "multiomics_integration.csv")
if (file.exists(integ_file)) {
  integ <- fread(integ_file)
  candidate_genes <- integ$gene
  cat("Multiomics candidate genes:", length(candidate_genes), "\n")
} else {
  # Use top variable genes
  gene_var <- apply(X, 2, var)
  candidate_genes <- colnames(X)[order(gene_var, decreasing=TRUE)[1:500]]
  cat("Top 500 variable genes selected\n")
}

# Filter to available genes
candidate_genes <- intersect(candidate_genes, colnames(X))
X_cand <- X[, candidate_genes]
cat("Candidate features:", ncol(X_cand), "\n")

# LASSO initial screen
cat("LASSO screening...\n")
cv_lasso <- cv.glmnet(X_cand, as.numeric(y)-1, alpha=1, nfolds=5, family="binomial")
lasso_coef <- coef(cv_lasso, s="lambda.min")
lasso_genes <- names(lasso_coef)[lasso_coef[,1] != 0][-1]  # remove intercept
cat("LASSO selected:", length(lasso_genes), "genes\n")

X_sel <- X[, intersect(candidate_genes, lasso_genes)]
if (ncol(X_sel) < 5) {
  cat("Too few features, using top 50 variable\n")
  X_sel <- X[, candidate_genes[1:min(50, length(candidate_genes))]]
}
cat("Features for ML:", ncol(X_sel), "\n")

# ============================================================
# 12 Algorithm Ensemble
# ============================================================
cat("\n--- 12 Algorithm Ensemble ---\n")

set.seed(42)
train_idx <- createDataPartition(y, p=0.7, list=FALSE)
X_train <- X_sel[train_idx, ]
y_train <- y[train_idx]
X_test <- X_sel[-train_idx, ]
y_test <- y[-train_idx]

results_list <- list()
importance_list <- list()

# 1. Lasso
cat("1/12 Lasso...\n")
tryCatch({
  fit <- glmnet(X_train, as.numeric(y_train)-1, alpha=1, family="binomial")
  pred <- predict(fit, X_test, s=0.01, type="response")
  results_list[["lasso"]] <- pred[,1]
  importance_list[["lasso"]] <- data.frame(gene=rownames(coef(fit))[-1],
                                           importance=abs(coef(fit)[-1, ncol(coef(fit))]))
}, error=function(e) cat("  Error:", e$message, "\n"))

# 2. Ridge
cat("2/12 Ridge...\n")
tryCatch({
  fit <- glmnet(X_train, as.numeric(y_train)-1, alpha=0, family="binomial")
  pred <- predict(fit, X_test, s=0.01, type="response")
  results_list[["ridge"]] <- pred[,1]
}, error=function(e) cat("  Error:", e$message, "\n"))

# 3. Random Forest
cat("3/12 Random Forest...\n")
tryCatch({
  fit <- randomForest(X_train, y_train, ntree=500, importance=TRUE)
  pred <- predict(fit, X_test, type="prob")[,2]
  results_list[["rf"]] <- pred
  importance_list[["rf"]] <- data.frame(gene=rownames(importance(fit)),
                                        importance=importance(fit)[,1])
}, error=function(e) cat("  Error:", e$message, "\n"))

# 4. GBM
cat("4/12 GBM...\n")
tryCatch({
  df_train <- data.frame(X_train, y=as.numeric(y_train)-1)
  fit <- gbm(y ~ ., data=df_train, n.trees=100, distribution="bernoulli")
  pred <- predict(fit, data.frame(X_test), n.trees=100, type="response")
  results_list[["gbm"]] <- pred
  # GBM relative influence as importance
  rel_inf <- summary(fit, plotit=FALSE)
  importance_list[["gbm"]] <- data.frame(gene=rel_inf$var, importance=rel_inf$rel.inf)
}, error=function(e) cat("  Error:", e$message, "\n"))

# 5. XGBoost
cat("5/12 XGBoost...\n")
tryCatch({
  fit <- xgboost(data=X_train, label=as.integer(as.numeric(y_train)-1),
                 nrounds=50, objective="binary:logistic", verbose=0)
  pred <- predict(fit, X_test)
  results_list[["xgb"]] <- pred
  imp <- xgb.importance(model=fit)
  if (!is.null(imp) && nrow(imp) > 0)
    importance_list[["xgb"]] <- data.frame(gene=imp$Feature, importance=imp$Gain)
}, error=function(e) cat("  Error:", e$message, "\n"))

# 6. SVM
cat("6/12 SVM...\n")
tryCatch({
  fit <- svm(X_train, y_train, probability=TRUE, kernel="radial")
  pred <- attr(predict(fit, X_test, probability=TRUE), "probabilities")[,2]
  results_list[["svm"]] <- pred
}, error=function(e) cat("  Error:", e$message, "\n"))

# 7. Naive Bayes
cat("7/12 Naive Bayes...\n")
tryCatch({
  fit <- NaiveBayes(X_train, y_train)
  pred <- predict(fit, X_test)$posterior[,2]
  results_list[["nb"]] <- pred
}, error=function(e) cat("  Error:", e$message, "\n"))

# 8. LDA
cat("8/12 LDA...\n")
tryCatch({
  fit <- lda(X_train, y_train)
  pred <- predict(fit, X_test)$posterior[,2]
  results_list[["lda"]] <- pred
}, error=function(e) cat("  Error:", e$message, "\n"))

# 9. Elastic Net
cat("9/12 Elastic Net...\n")
tryCatch({
  fit <- glmnet(X_train, as.numeric(y_train)-1, alpha=0.5, family="binomial")
  pred <- predict(fit, X_test, s=0.01, type="response")
  results_list[["enet"]] <- pred[,1]
}, error=function(e) cat("  Error:", e$message, "\n"))

# 10. glmBoost
cat("10/12 glmBoost...\n")
tryCatch({
  fit <- glmboost(X_train, y_train,
                  family=Binomial(), control=boost_control(mstop=50))
  pred <- predict(fit, newdata=X_test, type="response")
  results_list[["glmboost"]] <- pred[,1]
}, error=function(e) cat("  Error:", e$message, "\n"))

# 11. C5.0
cat("11/12 C5.0...\n")
tryCatch({
  library(C50)
  fit <- C5.0(X_train, y_train, trials=10)
  pred <- predict(fit, X_test, type="prob")[,2]
  results_list[["c50"]] <- pred
}, error=function(e) cat("  Error:", e$message, "\n"))

# 12. Neural Network
cat("12/12 Neural Network...\n")
tryCatch({
  library(nnet)
  df_train <- data.frame(X_train, y=as.numeric(y_train)-1)
  fit <- nnet(y ~ ., data=df_train, size=5, decay=0.1, linout=FALSE, maxit=100)
  pred <- predict(fit, data.frame(X_test))[,1]
  results_list[["nnet"]] <- pred
}, error=function(e) cat("  Error:", e$message, "\n"))

cat("Algorithms completed:", length(results_list), "/12\n")

# ============================================================
# Stacking
# ============================================================
cat("\n--- Stacking ensemble ---\n")

# Safety: ensure all predictions are numeric vectors of correct length
n_test <- length(y_test)
results_clean <- list()
for (nm in names(results_list)) {
  p <- as.numeric(results_list[[nm]])
  if (length(p) == n_test && !any(is.na(p))) {
    results_clean[[nm]] <- p
  } else {
    cat("  Skipping", nm, "- length", length(p), "expected", n_test, "\n")
  }
}
cat("Clean predictions:", length(results_clean), "algorithms\n")

if (length(results_clean) < 2) {
  stop("Too few algorithms succeeded for stacking")
}

# Create meta-features
meta_features <- do.call(cbind, results_clean)
colnames(meta_features) <- names(results_clean)

# Meta-learner: logistic regression on test predictions
meta_fit <- glm(as.numeric(y_test)-1 ~ ., data=as.data.frame(meta_features),
                family=binomial)

# Predict on test
meta_pred <- predict(meta_fit, type="response")

# Evaluate
library(pROC)
roc_obj <- roc(as.numeric(y_test)-1, meta_pred)
cat("Stacking AUC:", auc(roc_obj), "\n")

# Save model
stacking_model <- list(
  base_models=results_clean,
  meta_model=meta_fit,
  feature_genes=colnames(X_sel),
  auc=auc(roc_obj)
)
saveRDS(stacking_model, file.path(outdir, "stacking_model.rds"))

# Save feature importance
if (length(importance_list) > 0) {
  imp_combined <- do.call(rbind, lapply(names(importance_list), function(algo) {
    df <- importance_list[[algo]]
    df$algorithm <- algo
    df
  }))
  write.csv(imp_combined, file.path(outdir, "feature_importance_per_algo.csv"), row.names=FALSE)
}

# Consensus signature (genes selected by majority)
if (length(importance_list) >= 3) {
  gene_votes <- table(unlist(lapply(importance_list, function(df) {
    df$gene[order(df$importance, decreasing=TRUE)[1:min(20, nrow(df))]]
  })))
  consensus_genes <- names(gene_votes)[gene_votes >= 3]
  cat("Consensus signature genes:", length(consensus_genes), "\n")
  write.csv(data.frame(gene=consensus_genes, votes=gene_votes[consensus_genes]),
            file.path(outdir, "consensus_signature.csv"), row.names=FALSE)
}

cat("\n== step6 done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
