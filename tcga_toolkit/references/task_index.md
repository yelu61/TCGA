# Task Index

Detailed schema, inputs, outputs and gotchas for every task understood by
[`run_task.R`](../scripts/run_task.R). Run `Rscript tcga_toolkit/scripts/list_tasks.R --json` to get
the same information as machine-readable JSON.

Conventions:

- `required` — must be in the config or the task aborts.
- `optional` — has a sensible default.
- Outputs always land under `tcga_runs/<task_id>/{results,plots,objects}/`.

---

## audit_data

Inventory local TCGA / GTEx files and flag projects with expression `*.rda`
older than the GDC release date you cite.

- required: none
- optional: `gdc_release_version` (default `45.0`), `gdc_release_date`
  (default `2025-12-04`)
- outputs: `project_audit.csv`, `gtex_audit.csv`, `summary.json`

Use first whenever you suspect a project may be stale before chaining
DEG / survival downstream.

## prepare_bulk_rna

Materialise a per-project bundle (counts + tpm + clinical) from
`0-Data/<project>_mrna.rda` and cache it in `.tcga_cache/`.

- required: `project`
- outputs: cached object in `.tcga_cache/<project>_prepared.rds` (subsequent
  tasks reuse it transparently).

## run_deg

DESeq2 differential expression between two clinical groups.

- required: `project`, `numerator`, `denominator`
- optional: `group_column` (default `Group`), `padj_cutoff` (0.05),
  `log2fc_cutoff` (1), `min_count` (10), `min_fraction` (0.25)
- outputs: `<project>_<num>_vs_<den>_deg.csv`, `_volcano.png`/`.pdf`
- gotcha: if `Group` is missing, the toolkit derives it from
  `definition`/`sample_type`. Inspect with `inspect_project.R` first.

## run_enrichment

ORA or GSEA against local GMT files in `1-Input/`.

- required: `input_deg` (path to a DEG CSV produced by `run_deg`)
- optional: `mode` (`ora` | `gsea`, default `ora`), `gmt_paths`,
  `padj_cutoff`, `log2fc_cutoff`
- outputs: one `<gmt_stem>_<mode>_enrichment.csv` per GMT collection plus
  dotplots.

## run_clinical_assoc

Associate a gene's expression or a signature score with clinical variables
(non-parametric tests) and overall survival (Cox + median-split KM).

- required: `project` and one of `gene` / `signature_file`
- optional: `clinical_columns`, `signature_method` (`mean_z` default)
- outputs: `<project>_clinical_associations.csv`, `<project>_survival_cox.csv`,
  `<project>_clinical_boxplot.{png,pdf}`, `<project>_clinical_km.{png,pdf}`.

## run_tme

IOBR deconvolution of bulk RNA into tumour microenvironment fractions.

- required: `project`
- optional: `methods` (default `c("estimate")`; can include `cibersort`,
  `quantiseq`, etc), `perm` (only used by CIBERSORT, default 100)
- outputs: one `<project>_<method>_tme.csv` per method.
- requires the GitHub package `IOBR/IOBR`; run `check_deps.R --task run_tme`.

## run_gtex_compare

Compare TCGA tumour/normal expression with matched GTEx tissue from
`GTEX/gtex_RSEM_gene_tpm.gz`.

- required: `project`, `gene`
- optional: `gtex_tissues` (defaults to `GTEX/samplepair.txt` lookup)
- outputs: long-format and summary CSVs plus a grouped boxplot.

## validate_sc_signatures

Score a user-supplied single-cell signature on a TCGA cohort and test
clinical / survival associations of the score.

- required: `project`, `signature_file` (gene column required, optional
  `weight`, `direction`)
- optional: `signature_method`, `clinical_columns`
- outputs: signature score CSV, association CSV, survival CSV, boxplot,
  KM curve.

## pan_cancer_expression

Pan-cancer comparison of tumour vs normal for a gene or signature.

- required: `gene` or `signature_file`
- optional: `projects` (defaults to all TCGA-* projects with data)
- outputs: per-sample values, project summary (with tumour-vs-normal
  Wilcoxon), pan-cancer boxplot.

## stage_analysis

Test gene/signature against AJCC pathologic stage in a single project.

- required: `project` and one of `gene` / `signature_file`
- optional: `stage_column` (default `ajcc_pathologic_stage`)
- outputs: stage values CSV, stage summary CSV, stage-stratified boxplot.

## survival_map

Pan-cancer overall-survival map (median split + Cox per project).

- required: `gene` or `signature_file`
- optional: `projects`
- outputs: per-project Cox table, KM heatmap-style figure.

## subtype_analysis

Compare gene/signature across pan-cancer molecular subtypes
(`0-Data/PanCancer_subtypes.rda`).

- required: `project` and one of `gene` / `signature_file`
- optional: `subtype_column` (default `Subtype_Selected`)
- outputs: values CSV, subtype summary, Kruskal stat CSV, boxplot.

## gene_correlation_heatmap

Pearson/Spearman correlation between one target gene and a gene list
(e.g. an Excel sheet) across one or more projects.

- required: `target_gene`, `gene_list_file` (CSV/TSV/XLSX)
- optional: `sheet`, `gene_column`, `projects`, `method` (`pearson` |
  `spearman`), `sample_filter` (e.g. `c("Tumor")`)
- outputs: per-pair correlation CSV, project x gene heatmap.

## pipeline

Run several tasks sequentially with output passing.

- required: `steps` (array of task configs)
- variable substitution available in any string field of later steps:
  - `{{stepN.run_dir}}`, `{{stepN.results_dir}}`,
    `{{stepN.plots_dir}}`, `{{stepN.objects_dir}}`,
    `{{stepN.task_id}}`
  - `{{stepN.glob:<glob>}}` resolves to the first file matching the glob
    inside step N's `results_dir` (then `objects_dir`).

## render_report

Bundle several `tcga_runs/<task_id>/` directories into one
self-contained `report.html` (markdown + tables + base64-embedded PNGs).

- required: `run_dirs`
- optional:
  - `format` — `html` (default) or `quarto`.
  - When `format: "quarto"`:
    - `output_formats` — array of Quarto output formats; defaults to
      `["html", "pdf"]`. Other valid values: `docx`.
    - `title`, `author` — set the Quarto front-matter.
    - Falls back to bundled HTML when the `quarto` binary is not on PATH.
- outputs:
  - HTML mode: `report.html` at the run root.
  - Quarto mode: `<task_id>_report.qmd` plus rendered `.html` / `.pdf`
    files in the same directory.

## maf_summary

maftools-based summary of mutation burden, oncoplot and top genes.

- required: `project`
- optional: `top_n` (default 20)
- outputs: per-gene mutation frequency CSV, oncoplot PDF/PNG.

## mutation_survival

Mutant vs wildtype overall-survival comparison for a single gene.

- required: `project`, `gene`
- outputs: per-sample mutation status CSV, KM curve, Cox table.

## prognostic_model

Multivariable prognostic model from a panel of features.

- required: `project` and one of `feature_genes`, `signature_file`,
  `gene_list_file`
- optional: `method` (`lasso_cox` default | `rsf`), `test_fraction` (0.3),
  `seed` (20260512), `time_points` (`[365,1095,1825]`), `covariates`
  (clinical columns to keep as covariates), `ntree`, `nodesize`
- outputs:
  - `<project>_prognostic_coefficients.csv` — non-zero coefficients (lasso)
    or variable importance (RSF).
  - `<project>_prognostic_risk_scores.csv` — per-sample risk, with
    `split = train|test`.
  - `<project>_prognostic_performance.csv` — c-index per split.
  - `<project>_prognostic_timeROC_auc.csv` — time-dependent AUC (if
    `timeROC` is installed and method is `lasso_cox`).
  - `objects/<project>_prognostic_model.rds` — model object.
- gotcha: requires ≥30 samples with complete survival + expression. The
  task filters to `Group == "Tumor"` when at least 30 tumour samples are
  available; otherwise it keeps the full cohort.

## ssgsea_score

Sample-by-pathway scoring with GSVA (ssGSEA / GSVA / z-score / PLAGE).

- required: `project`, `gmt_paths`
- optional: `method` (`ssgsea` default), `min_size` (10), `max_size` (500)
- outputs: one `<project>_<gmt_stem>_<method>.csv` per collection. Rows
  are samples; columns are gene sets.
- works with both GSVA ≥1.50 (param API) and older versions.

## tmb_analysis

Per-sample tumour mutation burden and TMB-stratified survival.

- required: `project`
- optional:
  - `tmb_source`: `auto` (default; PanCancer_TMB → MAF fallback),
    `pancancer`, or `maf`.
  - `tmb_threshold`: explicit cutoff in mutations/Mb. Default = 75th
    percentile of the cohort.
  - `exome_size_mb`: used when computing from MAF (default 38 Mb).
- outputs: per-sample TMB CSV, cohort summary, TMB-stratified Cox + KM.

## cohort_qc

Per-sample QC for a project.

- required: `project`
- optional: `outlier_sd` (3), `min_library_size` (1e6), `batch_column`,
  `run_combat_seq` (boolean; needs `sva`)
- outputs: per-sample QC CSV (library size, detection rate, PCA flag),
  PCA scatter plot, library size histogram, optional ComBat-seq adjusted
  counts under `objects/`.

## cnv_summary

Gene-level copy-number burden, top amplifications and deletions, and an
optional CNV-burden Cox regression.

- required: `project`
- optional:
  - `cnv_file` — gene x sample matrix (`.rda`, `.tsv`, `.csv`). Discrete
    GISTIC-style values (`-2..+2`) or continuous log2 ratios both work.
    Auto-detected in `0-Data/<project>_cnv.rda` /
    `<project>_gistic.rda` / `.tsv` / `.csv`.
  - `top_n` (default 25).
  - `amp_threshold` / `del_threshold` — used when values are continuous
    (default ±0.3 log2).
  - `focal_threshold` — used when values are discrete (default 1).
- outputs:
  - `<project>_cnv_gene_frequency.csv` — amp/del frequency per gene.
  - `<project>_cnv_top_amplifications.csv`, `<project>_cnv_top_deletions.csv`.
  - `<project>_cnv_per_sample.csv` — per-sample amp/del count + burden.
  - `<project>_cnv_burden_cox.csv` — Cox regression on CNV burden if survival
    metadata is available.
  - Bar plots for top amp / top del + a CNV-burden histogram.
- gotcha: the task aborts cleanly when no CNV file is found and prints the
  TCGAbiolinks query to fetch one.

## methylation_diff

Differential methylation between two clinical groups on a user-provided
beta-value matrix.

- required: `project`
- optional:
  - `meth_file` — CpG x sample matrix (`.rda`, `.tsv`, `.csv`). 0..1 betas
    or 0..100 percentages are both accepted (auto-rescaled).
  - `group_column` (default `Group`), `numerator` (`Tumor`),
    `denominator` (`Normal`).
  - `delta_beta_cutoff` (0.1), `padj_cutoff` (0.05), `top_n` (50).
- outputs:
  - `<project>_methylation_diff.csv` — per-CpG mean β, Δβ, p-value, padj.
  - `<project>_methylation_top_hyper.csv`,
    `<project>_methylation_top_hypo.csv`.
  - `<project>_methylation_volcano.png/pdf`.
- gotcha: per-CpG Wilcoxon on a full 450k panel is slow (5–15 min). Use a
  smaller probe subset (e.g. promoter-only) for faster iteration.

## external_validate

Validate a signature or risk score on an external cohort (CGGA, GEO,
METABRIC, ICGC, ...). Either supply a `signature_file` (gene list, used
with z-score mean) or a `weight_file` (gene + coefficient, used as a
weighted linear score). The latter is typically the
`<project>_prognostic_coefficients.csv` written by `prognostic_model`.

- required: `cohort_name`, `expression_file`, `clinical_file`
- optional: `feature_genes`, `signature_file`, `weight_file`
- expression_file: gene x sample matrix; gene IDs must be symbols (or
  match the supplied signature/weight rows).
- clinical_file: must contain a sample-id column (one of `sample_id`,
  `Sample`, `barcode`, ...) plus a survival time column
  (e.g. `OS.time`) and an event column (e.g. `OS`).
- outputs: per-sample risk + survival CSV, Cox HR + 95% CI + concordance,
  median-split KM curve.

## immune_phenotype

Immune phenotype scoring. Primary path uses IOBR TME deconvolution
(ESTIMATE + CIBERSORT by default), joined with the pan-cancer immune
subtype label from `0-Data/PanCancer_subtypes.rda`. When IOBR is not
installed the task falls back to ssGSEA on the immune-related sets in the
Hallmark GMT.

- required: `project`
- optional: `tme_methods`, `perm`, `immune_subtype_file`, `immune_gmt`
- outputs: per-sample immune metric matrix, summary statistics, optional
  immune-subtype scatter (`StromalScore` vs `ImmuneScore`).

## wgcna_modules

WGCNA co-expression module detection + module-trait correlation.

- required: `project`
- optional:
  - `n_top_genes` (5000), `min_module_size` (30), `merge_cut_height` (0.25).
  - `soft_power` — pass an integer to skip auto-selection.
  - `sample_filter` — e.g. `["Tumor"]` to restrict to tumour samples.
  - `trait_columns` — clinical columns to correlate against module
    eigengenes (default `Group`, `age_at_index`, `ajcc_pathologic_stage`).
- outputs:
  - `<project>_wgcna_module_assignment.csv` — gene → module label/color.
  - `<project>_wgcna_module_sizes.csv` — module size summary.
  - `<project>_wgcna_module_eigengenes.csv` — sample × module matrix.
  - `<project>_wgcna_module_trait_correlation.csv` — when trait_columns
    yield enough non-missing rows.
  - `objects/<project>_wgcna_network.rds` — the full network object.
- gotcha: WGCNA is memory-hungry on large gene sets; `n_top_genes = 5000`
  is the recommended ceiling on a laptop.

## drug_response

Predict per-sample drug IC50 using oncoPredict against a GDSC or CTRP
training set, then test for difference by `Group`.

- required: `project`
- optional:
  - `training_dir` — points to the unpacked oncoPredict training data.
    Auto-detected in `1-Input/oncoPredict_training_data/` and
    `1-Input/DataFiles/Training Data/`.
  - `drug_set` — `GDSC2` (default), `GDSC1`, or `CTRP2`.
- outputs:
  - `<project>_<drug_set>_drug_predictions.csv` — sample × drug matrix.
  - `<project>_<drug_set>_drug_summary.csv` — per-drug median/mean/SD.
  - `<project>_<drug_set>_drug_group_diff.csv` — Kruskal test by `Group`
    (when ≥2 groups present).
  - Top differential drugs boxplot.
- gotcha: oncoPredict's `calcPhenotype` can take 10–60 minutes depending
  on the cohort and drug set; consider running with `--overwrite` only
  after a smaller test cohort succeeds.
