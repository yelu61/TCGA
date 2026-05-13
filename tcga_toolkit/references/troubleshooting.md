# Troubleshooting

Lookup table for the most common failure modes when running tasks.

---

## Missing R packages

Error pattern: `<task> requires the <pkg> package.` or `there is no package called '...'`.

Resolution:

```bash
Rscript tcga_toolkit/scripts/check_deps.R --task <task>
```

Then run the install command the script emits (CRAN / Bioconductor / GitHub).

Notes:

- `run_tme` depends on the GitHub package `IOBR/IOBR` — install via
  `remotes::install_github("IOBR/IOBR")`.
- `ssgsea_score` works against both old (`min.sz=`) and new (`*Param`)
  GSVA APIs; if you upgrade GSVA the task auto-adapts.

## "Gene X not found in project Y"

Causes:

- The project's `gene_map` is keyed on Ensembl IDs and the user gave a
  symbol that's not in the map.
- The gene was renamed; see `gene_alias_map` in
  [`scripts/common.R`](../scripts/common.R).

Resolution:

- Run `inspect_project.R --project Y` to confirm the project loads.
- Try the Ensembl ID, the official HGNC symbol, or add the alias to
  `gene_alias_map`.

## "No genes passed the count filter"

Defaults: `min_count = 10` reads in at least `min_fraction = 0.25` of samples.

Resolution:

- For small cohorts (e.g. TARGET-OS, TCGA-CHOL, TCGA-DLBC), lower
  `min_fraction` to 0.1 in the `run_deg` config.
- Verify the project actually has a meaningful Tumor/Normal split with
  `inspect_project.R` — many TCGA cohorts have ≤ 5 normals.

## "No subtype-matched samples available for X"

Causes:

- `0-Data/PanCancer_subtypes.rda` may not cover the project, or
- The `cancer.type` field uses a different code than expected
  (e.g. `KIRC` vs `TCGA-KIRC`).

Resolution:

- Open `PanCancer_subtypes.rda` and inspect the `cancer.type` levels.
- Override `subtype_column` if your cohort uses a custom column name.

## GTEx tissue lookup returns empty

Resolution:

- Confirm `GTEX/samplepair.txt` contains your project. If not, pass
  `gtex_tissues` explicitly in the config.

## Pipeline variable did not resolve

Common error: `Pipeline glob pattern '*.csv' matched no files in step N
results/objects`.

Causes:

- Earlier step failed silently before writing outputs.
- Glob pattern is too restrictive (e.g. you typed
  `{{step0.glob:*_deg_filtered.csv}}` but `run_deg` writes
  `*_deg.csv`).

Resolution:

- Read `tcga_runs/<task_id>/run_metadata.json` for each step to confirm
  `status: "success"`.
- Use `inspect_project.R` and `list_runs.R` to confirm.

## Out-of-memory when loading a large project

Symptoms: R aborts while loading `TCGA-BRCA_mrna.rda` or similar.

Resolution:

- Cached prepared bundle in `.tcga_cache/` is read on subsequent runs;
  make sure you don't keep multiple R sessions live at once.
- For pure exploration, prefer `inspect_project.R` (which only reports
  metadata) over manually `load()`-ing the rda.

## "No PanCancer TMB record and no MAF file"

`tmb_analysis` first tries `PanCancer_TMB.Rdata`, then falls back to the
project MAF. If both are missing the task aborts.

Resolution:

- Confirm `0-Data/<project>_maf.rda` exists, and that `maftools` is
  installed.

## Stale `*_mrna.rda` after GDC refresh

`audit_data` flags `mrna_stale_vs_gdc45 = TRUE` when the file pre-dates
the GDC release you specify. The toolkit will still run, but you should
refresh the project's rda from `GDCdata/` to avoid downstream mismatches.

## Cache invalidation

The cache lives at `.tcga_cache/<project>_prepared.rds`. It is auto-
invalidated when the underlying `*_mrna.rda` is newer. Force a rebuild by
deleting the cache file.

## "cnv_summary requires a CNV matrix"

`cnv_summary` and `methylation_diff` deliberately fail clearly when no
data file is found rather than silently fabricating one.

Resolution:

- Place the matrix at `0-Data/<project>_cnv.rda` (CNV) or
  `0-Data/<project>_methylation.rda` (methylation), **or**
- Pass `cnv_file` / `meth_file` in the config pointing at any location.
- See [data_inventory.md](data_inventory.md) for the TCGAbiolinks query
  to fetch them.

## "external_validate matched only X of Y features"

The signature / weight file refers to symbols (or IDs) that don't exist
in the external cohort's expression matrix.

Resolution:

- Make sure both matrices use the same identifier type (HGNC symbol vs
  Ensembl). Bridge via the `gencode.v23.annotation.gene.probemap` in
  `GTEX/` if needed.
- Lowering the matching threshold rarely helps — if fewer than half of
  the features are present, the validation will be biased.

## "drug_response needs GDSC training data"

`oncoPredict` is not bundled. Download once and the task auto-detects it
in `1-Input/oncoPredict_training_data/`.

## WGCNA out-of-memory

WGCNA stores a dense similarity matrix in RAM. With 16 GB of memory
`n_top_genes ≥ 8000` typically OOMs.

Resolution:

- Lower `n_top_genes` to 5000 (default) or 3000.
- Restrict to tumour samples via `sample_filter: ["Tumor"]`.
- Run on a machine with more RAM.

## Quarto render fails

`render_report` with `format: "quarto"` requires the standalone Quarto
binary on PATH and (for PDF) a working LaTeX installation.

Resolution:

- Verify with `quarto --version`. If missing, install from
  <https://quarto.org/>.
- For PDF output also need TinyTeX (`quarto install tinytex`).
- The task falls back to the bundled HTML when Quarto is unavailable —
  check the run report to confirm which mode was used.

## Re-running with `--overwrite`

By default each run gets a fresh timestamped directory. If you want to
overwrite the previous run for the same `task_id`:

```bash
Rscript tcga_toolkit/scripts/run_task.R --config X.json --overwrite
```

This deletes and recreates `tcga_runs/<task_id>/`.
