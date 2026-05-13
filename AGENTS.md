# Repository Guidelines

## Project Structure & Module Organization
This repository uses `tcga_toolkit/` as the canonical entrypoint for all TCGA and GTEx analysis workflows. The legacy R notebooks (`0-TCGA_DATAdownload.ipynb`, `1-TCGA_DEGanalysis.ipynb`, `2-TCGA_GENEanalysis.ipynb`, `3-TCGA_SNVanalysis.ipynb`) and `deg_analysis.R` are retained for historical reference but should not be used for new analyses.

Primary data directories are path-sensitive:
- `0-Data/`: serialized `.rda` datasets (expression, clinical, MAF, pan-cancer subtypes)
- `1-Input/`: curated GMT collections and pre-computed DEG inputs
- `GDCdata/`: raw GDC downloads by project
- `GTEX/`: GTEx matrix (`gtex_RSEM_gene_tpm.gz`), tissue annotation, and TCGA-GTEx tissue mapping
- `tcga_runs/`: all toolkit outputs, organized by `<timestamp>-<task_id>/`

## Build, Test, and Development Commands
Run commands from the repository root so relative paths resolve correctly.

- `Rscript tcga_toolkit/scripts/run_task.R --config tcga_toolkit/templates/<config>.json`: run a toolkit task
- `R -q -e "source('tcga_toolkit/scripts/common.R')"`: validate common helpers before running tasks
- `jupyter lab`: open legacy notebooks (read-only reference)

## Coding Style & Naming Conventions
Follow the existing R style in `tcga_toolkit/scripts/`: 2-space indentation, `<-` for assignment, and snake_case for functions. Keep task scripts self-contained: accept `config` + `ctx`, write outputs to `ctx$results_dir`, `ctx$plots_dir`, and `ctx$objects_dir`, and finish with `write_report()`.

Preserve naming patterns already used in the repo: cancer-specific files such as `TCGA-BRCA_mrna.rda`, task configs such as `deg_target_os.json`, and output names such as `TARGET-OS_CD24_gsea_KEGG.csv`.

## Testing Guidelines
There is no automated test suite. Validate changes by running the affected task config and checking that expected artifacts appear in `tcga_runs/<task_id>/`. For differential analysis changes, verify sample alignment, factor levels, and contrast labels before saving results.

## Commit & Pull Request Guidelines
Use short, imperative commit subjects such as `Add subtype_analysis task` or `Fix GTEx tissue lookup fallback`.

For pull requests, identify the affected task(s), list regenerated outputs, and include figure previews when plots change. Avoid committing notebook output churn or OS metadata such as `.DS_Store`.
