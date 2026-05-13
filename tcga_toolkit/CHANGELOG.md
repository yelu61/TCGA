# Changelog

All notable changes to the TCGA toolkit are tracked here.

## 0.3.0 — 2026-05-13

### Added — analysis coverage

- `cnv_summary`: gene-level copy-number burden, top amplifications /
  deletions, and CNV-burden survival.
- `methylation_diff`: differential methylation between two clinical
  groups on a user-provided beta matrix (450k or EPIC).
- `external_validate`: validate a signature or risk score on an external
  cohort (CGGA, GEO, METABRIC, ICGC, ...). Accepts either a feature list
  or a weight CSV produced by `prognostic_model`.
- `immune_phenotype`: IOBR TME deconvolution joined with the PanCancer
  immune subtype; falls back to ssGSEA on immune Hallmark sets when
  IOBR is missing.
- `wgcna_modules`: WGCNA co-expression module detection with optional
  module-trait correlation.
- `drug_response`: oncoPredict drug IC50 prediction against GDSC2 / GDSC1
  / CTRP2, with group-difference testing.

### Added — infrastructure

- `task_render_report.R` gained a `format: "quarto"` mode that emits a
  Quarto `.qmd` and renders to HTML + PDF (configurable). Falls back to
  the bundled HTML renderer when the Quarto binary is unavailable.
- Python adapter under `tcga_toolkit/adapters/python/`:
  - `runs.py` reads any `tcga_runs/<task_id>/` directory.
  - `survival.py` wraps scikit-survival (Cox, KM, time-dependent AUC).
  - `shap_explainer.py` provides SHAP explainability via a LightGBM
    surrogate model fitted on exported risk scores.
- Smoke test suite under `tcga_toolkit/tests/`:
  - `test_introspection.R`, `test_templates.R`, `test_dispatch.R`,
    `test_python_adapter.py`, plus `run_all.sh` runner.

### Added — templates

- `brca_cnv_summary.json`, `brca_methylation_diff.json`,
  `cgga_external_validate.json`, `brca_immune_phenotype.json`,
  `brca_wgcna.json`, `brca_drug_response.json`, `fa_quarto_report.json`.

### Changed

- `specs/README.md`, `references/task_index.md`,
  `references/pipeline_recipes.md`, `references/data_inventory.md`,
  `references/troubleshooting.md`, `references/README.md`, `SKILL.md`
  expanded to cover the new tasks, Quarto report mode, optional CNV /
  methylation / external data layout, and the smoke test runner.
- `validate_config.R` accepts the new task names. `list_tasks.R` and
  `check_deps.R` enumerate them with the right dependency lists.

### Notes

- New tasks degrade cleanly when their auxiliary data is absent: CNV,
  methylation, external cohort and oncoPredict training data are all
  optional, and the task aborts with a clear message + the
  TCGAbiolinks / OSF URL to fetch them.

## 0.2.0 — 2026-05-12

### Added
- `tcga_toolkit/VERSION` file; `toolkit_version` is now written into every `run_metadata.json`.
- New tasks:
  - `prognostic_model`: Lasso-Cox / Random Survival Forest with train-test split, time-dependent ROC, risk score export.
  - `ssgsea_score`: GSVA / ssGSEA / single-sample pathway scoring against any GMT collection.
  - `tmb_analysis`: tumour mutation burden summary and survival stratification (uses local MAF or `0-Data/PanCancer_TMB.Rdata`).
  - `cohort_qc`: per-project QC (sample counts, library size, gene detection, PCA outlier flag, optional ComBat-seq batch correction).
- Introspection scripts under `tcga_toolkit/scripts/`:
  - `list_tasks.R` — print every task plus required/optional fields.
  - `inspect_project.R` — summarise available clinical columns, subtypes, sample counts.
  - `validate_config.R` — dry-run config schema validation (no analysis side effects).
  - `list_runs.R` — list `tcga_runs/` contents with status filter.
  - `check_deps.R` — detect missing R/Bioconductor packages per task and emit install commands.
- Expanded `tcga_toolkit/references/` into per-topic files: `task_index.md`, `pipeline_recipes.md`, `data_inventory.md`, `troubleshooting.md`, `extension_guide.md`.
- New config templates for the four new tasks under `tcga_toolkit/templates/`.

### Changed
- `SKILL.md` rewritten with a decision tree, cookbook recipes, an explicit list of when *not* to use the skill, and an expanded keyword surface for trigger matching.
- `specs/README.md` lists the new tasks and their required fields.
- `run_task.R` registers the new task dispatchers and propagates `toolkit_version` into run metadata.

### Notes
- Existing run directories under `tcga_runs/` remain compatible. The new `toolkit_version` field is additive.

## 0.1.0 — 2026-04-16

Initial public toolkit: `audit_data`, `prepare_bulk_rna`, `run_deg`, `run_enrichment`, `run_clinical_assoc`, `run_tme`, `run_gtex_compare`, `validate_sc_signatures`, `pan_cancer_expression`, `stage_analysis`, `survival_map`, `subtype_analysis`, `gene_correlation_heatmap`, `pipeline`, `render_report`, `maf_summary`, `mutation_survival`.
