# Step 9: Nomogram + Decision Curve Analysis + Shiny App
# Input: stacking_model.rds + clinical data
# Output: nomogram.pdf, dca_plot.pdf, shiny_app/app.R
# Keep this file pure ASCII.

Sys.setenv(R_LIBS_USER="D:/claude-study/Rlib/4.6",
           TMP="D:/claude-study/Rtmp",
           TEMP="D:/claude-study/Rtmp",
           TMPDIR="D:/claude-study/Rtmp")
.libPaths(unique(c("D:/claude-study/Rlib/4.6", .libPaths(),
                   "C:/Users/nuo'nuo/AppData/Local/R/win-library/4.6")))

suppressPackageStartupMessages({
  library(rms)
  library(ggplot2)
})

outdir <- "D:/claude-study/Spinal Ependymoma/analysis/bulk"
sink(file.path(outdir, "step9_nomogram.log"), append=TRUE)
cat("== step9 start", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")

# ============================================================
# 9a: Nomogram
# ============================================================
cat("\n--- 9a: Nomogram ---\n")

# Load data
model_info <- readRDS(file.path(outdir, "stacking_model.rds"))
consensus <- read.csv(file.path(outdir, "consensus_signature.csv"))
signature_genes <- consensus$gene

expr_raw <- readRDS(file.path(outdir, "GSE242994_expr.rds"))
if (is.list(expr_raw) && "expr" %in% names(expr_raw)) { expr <- expr_raw$expr } else { expr <- expr_raw }
labels <- read.csv(file.path(outdir, "GSE242994_rna_labels.csv"))
common <- intersect(colnames(expr), labels$sample)

# Select top signature genes (limit to 5-8 for nomogram readability)
top_genes <- head(signature_genes, min(6, length(signature_genes)))
X <- t(as.matrix(expr[top_genes, common]))
y <- as.numeric(labels$Transcriptional_subtype[match(common, labels$sample)] == "B")

# Fit logistic regression
df <- as.data.frame(X)
df$subtype <- y

# Use rms for nomogram
dd <- datadist(df)
options(datadist="dd")

fit <- lrm(subtype ~ ., data=df[, c(top_genes, "subtype")], x=TRUE, y=TRUE)
cat("LR model R2:", fit$stats["R2"], "\n")

# Create nomogram
nom <- nomogram(fit, fun=function(x) 1/(1+exp(-x)),
                funlabel="P(Subtype B)",
                lp=TRUE)

pdf(file.path(outdir, "nomogram.pdf"), width=12, height=8)
plot(nom, cex.var=0.8, cex.axis=0.7)
dev.off()
cat("Nomogram saved\n")

# ============================================================
# 9b: Decision Curve Analysis (DCA)
# ============================================================
cat("\n--- 9b: DCA ---\n")

tryCatch({
  library(rmda)

  # Build explicit formula
  dca_formula <- as.formula(paste("subtype ~", paste(top_genes, collapse=" + ")))
  dca_obj <- decision_curve(dca_formula,
                            data=df[, c(top_genes, "subtype")],
                            family=binomial(link="logit"),
                            thresholds=seq(0, 1, by=0.01))

  pdf(file.path(outdir, "dca_plot.pdf"), width=8, height=6)
  plot_decision_curve(dca_obj,
                      cost.benefit.axis=TRUE,
                      legend.position="topright")
  dev.off()
  cat("DCA plot saved\n")
}, error=function(e) {
  cat("DCA error:", e$message, "\n")
})

# ============================================================
# 9c: Shiny App
# ============================================================
cat("\n--- 9c: Shiny app ---\n")

shiny_dir <- file.path(outdir, "shiny_app")
dir.create(shiny_dir, showWarnings=FALSE, recursive=TRUE)

# Write app.R
app_code <- sprintf('
library(shiny)
library(rms)

# Load model
load("model.rds")

ui <- fluidPage(
  titlePanel("Spinal Ependymoma Subtype Predictor"),

  sidebarLayout(
    sidebarPanel(
      h4("Enter Gene Expression Values"),
%s,
      actionButton("predict", "Predict", class="btn-primary"),
      hr(),
      h4("Prediction Result"),
      verbatimTextOutput("result")
    ),

    mainPanel(
      h3("Model Information"),
      p("This nomogram predicts the probability of transcriptional subtype B vs A."),
      h4("Signature Genes:"),
      verbatimTextOutput("genes"),
      h4("Instructions:"),
      p("1. Enter normalized expression values for each gene"),
      p("2. Click Predict to get the probability"),
      p("3. Probability > 0.5 suggests Subtype B")
    )
  )
)

server <- function(input, output, session) {
  output$genes <- renderPrint({
    cat(paste(genes, collapse="\\n"))
  })

  output$result <- renderPrint({
    input$predict
    isolate({
      newdata <- data.frame(
%s
      )
      pred <- predict(fit, newdata=newdata, type="fitted")
      cat(sprintf("Probability of Subtype B: %%.3f\\n", pred))
      if (pred > 0.5) {
        cat("Predicted: Subtype B\\n")
      } else {
        cat("Predicted: Subtype A\\n")
      }
    })
  })
}

shinyApp(ui, server)
',
paste(sapply(top_genes, function(g) sprintf('      numericInput("%s", "%s", value=0, min=-10, max=10, step=0.1)', g, g)), collapse=",\n"),
paste(sapply(top_genes, function(g) sprintf('        %s = input$%s', g, g)), collapse=",\n")
)

writeLines(app_code, file.path(shiny_dir, "app.R"))

# Save model for Shiny
save(list=c("fit", "top_genes"), file=file.path(shiny_dir, "model.rds"))
cat("Shiny app created at:", shiny_dir, "\n")

cat("\n== step9 done", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==\n")
sink()
