# Analysis code and processed tables for:
# "The deposited methylome separates SP-EPN subtypes A and B, whereas the deposited transcriptome does not: a leakage-free re-analysis"

This archive contains the reproducible analysis code and processed tables for a leakage-free re-analysis of public multi-omics data on spinal ependymoma (SP-EPN).

## Contents

```
code/               Analysis scripts (R 4.6.0 and one Python 3 helper)
processed_tables/   Supply material tables 1-8 (CSV) + table descriptions
```

## Data availability

All raw data are publicly available at NCBI GEO:

| Accession | Contents |
|---|---|
| GSE242991 | Discovery methylation matrix (450K/EPIC), 226 sample columns |
| GSE242992 | Discovery RNA-seq counts, 72 subtyped SP-EPN tumors (A/B labels) |
| GSE298501 | Single-nucleus RNA-seq (10x), SP-EPN and MPE tumors |
| GSE243240 | Validation methylation cohort (287 samples; 77 450K + 210 EPIC IDATs) |
| GSE66787, GSE64415, GSE125861, GSE66354 | Context cohorts referenced in the manuscript |

Raw data are not redistributed here; download from GEO. The 72-sample A/B labels are from the source study (Neyazi et al., *Nature Communications* 2023); the mapping between GEO samples and published sample IDs is reproduced in `processed_tables` and in `code/build_labels_v4.R` / `code/match_by_counts3.R`.

## Environment

- Windows 11 x64; R 4.6.0 (2026-04-24 ucrt); Rtools45; Bioconductor 3.23
- Python 3 (for `chr22_ci_mb.py` and `md_to_docx.py`)
- Key R packages: Seurat v5, SingleCellExperiment, DESeq2, edgeR, limma, minfi, methylKit, DMRcate, data.table, randomForest, pROC, caret, glmnet, xgboost, mclust, pheatmap, ComplexHeatmap, SCENIC/RcisTarget, CytoTRACE, Monocle 3, CellChat, DoubletFinder, scDblFinder, harmony, sctransform, clusterProfiler, missMethyl, survival, rms, shiny, NMF, fastshap
- Full installed-package inventory: see `R_pkg_logs/phaseE_verify.csv` in the project repository (GitHub link, see below)

**Note:** scripts contain absolute Windows paths (`D:/claude-study/...`); adjust the path variables at the top of each script before re-running. R is invoked via PowerShell (`powershell.exe -Command "& 'C:/Program Files/R/R-4.6.0/bin/Rscript.exe' 'script.R'"`).

## Script map (analysis modules)

| Script(s) | Result |
|---|---|
| `step1_qc.R`, `step2_integrate.R`, `step2b_markers_v6.R`, `step2c_annotate.R` | Single-nucleus QC, Seurat v5 integration, marker-based annotation of 10 cell types (Fig. 1A) |
| `scdblfinder_run.R` | Doublet detection (sensitivity analysis) |
| `cytotrace_by_group.R`, `step3_trajectory_scenic.R` | CytoTRACE / Monocle developmental-potential scores (Fig. 1B, Fig. 3A) |
| `step3c_scenic_only.R`, `fig1c_regulon_heatmap.R` | SCENIC regulon activity (Fig. 1C) |
| `step4_cellchat.R` | Cell-cell communication (not retained in final manuscript) |
| `build_242994.R`, `build_bulk.R`, `build_labels_v4.R`, `match_by_counts3.R`, `crosscheck_subtype.R` | Discovery bulk RNA/methylation matrices, 72-sample A/B label mapping |
| `deg_reanalysis.R`, `fig2a_final.R`, `fig2c_review_fix.R` | Bulk DEG re-analysis: 18,556 genes, 0 at FDR < 0.05 (Fig. 2A/2C) |
| `methylation_nested_cv.R`, `prep_ml_labels.R`, `fig2b_nested_roc.R` | Methylation nested CV (top-100 DMPs, 5 reps; pooled OOF AUC 0.956, CI 0.936-0.976) (Fig. 2B/2D) |
| `verify_metrics.R`, `verify_transcriptome_pooled.R` | Pooled OOF AUCs and threshold metrics for both modalities |
| `step5b_methylation_dmp.R`, `dmp_annotation_overlap.R`, `fig_chromosome.R`, `chr22_ci_mb.py`, `snp_fisher.R` | 1,722 DMPs, gene-region/island enrichment, chromosome 22 depletion (Fig. S1) |
| `gse243240_idat_processing.R`, `gse243240_noob_reprocess.R`, `gse243240_validation.R`, `gse243240_spine_subtypes.R`, `gse243240_consistent.R`, `parse_meta243240.R`, `noob_sens91.R`, `noob_sensitivity.R` | Validation cohort: label-free directional transfer of the A/B dichotomy; 14 overlap samples resolved as identical arrays (Fig. 3, Fig. S2) |
| `nullmodel_correct.R`, `noob_nullmodel.R`, `gse243240_random20k_noob.R` | Random-probe null models (20,000-probe noob pool; age 0/100, ARI median 0.002) |
| `pseudobulk_platform.R`, `tumor_frac_median.R`, `composition_check.R`, `singlecell_AB_comparison.R` | Pseudo-bulk testing, tumor-purity quantification, composition checks, single-nucleus A/B comparisons |
| `external_validation_v2.R`, `nf2_replication.R` | Cross-cohort structural replication, NF2 differential-expression replication |
| `m10_22q_regroup.R` | 22q/NF2-status regrouping sensitivity (201 DMPs vs 0 DEGs) |
| `review1010_analyses.R` | Detection rates of DEG panel genes; doublet-excluded CytoTRACE sensitivity |
| `fig3A3B_review_fix.R`, `fig3c_scores.R`, `manuscript_umap.R` | Figure 3 panels |
| `md_to_docx.py` | Manuscript Word generation from the Markdown source |

## License

- Code: MIT
- Processed tables: CC-BY 4.0

## Contact

[Corresponding author email - placeholder]
