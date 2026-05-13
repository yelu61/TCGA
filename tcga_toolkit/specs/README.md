# TCGA Toolkit Specs

All tasks are config-driven. Each run uses:

```bash
Rscript tcga_toolkit/scripts/run_task.R --config <config.json>
```

Required top-level field:

- `task`: one of `audit_data`, `prepare_bulk_rna`, `run_deg`, `run_enrichment`,
  `run_clinical_assoc`, `run_tme`, `run_gtex_compare`, `validate_sc_signatures`,
  `pan_cancer_expression`, `stage_analysis`, `survival_map`, `subtype_analysis`,
  `gene_correlation_heatmap`, `pipeline`, `render_report`, `maf_summary`,
  `mutation_survival`, `prognostic_model`, `ssgsea_score`, `tmb_analysis`,
  `cohort_qc`, `cnv_summary`, `methylation_diff`, `external_validate`,
  `immune_phenotype`, `wgcna_modules`, `drug_response`.

Common optional fields:

- `task_id`: custom suffix for the run folder
- `project`: TCGA or TARGET project id, e.g. `TCGA-BRCA` or `TARGET-OS`

For a fuller description of each task and its config schema, see
[../references/task_index.md](../references/task_index.md).

## Task-specific notes

- `audit_data`: optional `gdc_release_version`, `gdc_release_date`
- `run_deg`: requires `project`, `numerator`, `denominator`; optional
  `group_column`, `padj_cutoff`, `log2fc_cutoff`, `min_count`, `min_fraction`
- `run_enrichment`: requires `input_deg`; optional `mode` (`ora` or `gsea`)
  and `gmt_paths`
- `run_clinical_assoc`: requires `project` plus either `gene` or
  `signature_file`
- `run_tme`: requires `project`; optional `methods`, `perm`
- `run_gtex_compare`: requires `project`, `gene`; optional `gtex_tissues`
- `validate_sc_signatures`: requires `project`, `signature_file`; optional
  `signature_method`, `clinical_columns`
- `pan_cancer_expression`: requires `gene` or `signature_file`; optional
  `projects`
- `stage_analysis`: requires `project` plus `gene` or `signature_file`;
  optional `stage_column`
- `survival_map`: requires `gene` or `signature_file`; optional `projects`
- `subtype_analysis`: requires `project` plus `gene` or `signature_file`;
  optional `subtype_column`
- `gene_correlation_heatmap`: requires `target_gene`, `gene_list_file`;
  optional `sheet`, `gene_column`, `projects`, `method`, `sample_filter`
- `pipeline`: requires `steps` (array of task configs); each step can
  reference previous step outputs via `{{stepN.run_dir}}`,
  `{{stepN.results_dir}}`, `{{stepN.plots_dir}}`, `{{stepN.objects_dir}}`,
  `{{stepN.task_id}}`, or `{{stepN.glob:*_deg.csv}}`
- `render_report`: requires `run_dirs` (array of run directory paths);
  generates a unified `report.html`
- `maf_summary`: requires `project`; optional `top_n` (default 20)
- `mutation_survival`: requires `project` and `gene`
- `prognostic_model`: requires `project` and one of `feature_genes`,
  `signature_file`, or `gene_list_file`; optional `method` (`lasso_cox` |
  `rsf`), `test_fraction`, `seed`, `time_points`, `covariates`, `ntree`,
  `nodesize`
- `ssgsea_score`: requires `project`, `gmt_paths`; optional `method`
  (`ssgsea` | `gsva` | `zscore` | `plage`), `min_size`, `max_size`
- `tmb_analysis`: requires `project`; optional `tmb_source`
  (`auto` | `maf` | `pancancer`), `tmb_threshold`, `exome_size_mb`
- `cohort_qc`: requires `project`; optional `batch_column`, `outlier_sd`,
  `min_library_size`, `run_combat_seq`
- `cnv_summary`: requires `project`; optional `cnv_file` (path to gene-level
  CNV matrix), `top_n`, `amp_threshold`, `del_threshold`, `focal_threshold`
- `methylation_diff`: requires `project`; optional `meth_file` (beta-value
  matrix CpG x sample), `group_column`, `numerator`, `denominator`,
  `delta_beta_cutoff`, `padj_cutoff`, `top_n`
- `external_validate`: requires `cohort_name`, `expression_file`,
  `clinical_file`; optional `feature_genes`, `signature_file`,
  `weight_file` (CSV of gene / coefficient pairs)
- `immune_phenotype`: requires `project`; optional `tme_methods`, `perm`,
  `immune_subtype_file`, `immune_gmt`
- `wgcna_modules`: requires `project`; optional `n_top_genes`,
  `min_module_size`, `merge_cut_height`, `soft_power`, `sample_filter`,
  `trait_columns`
- `drug_response`: requires `project`; optional `training_dir`,
  `drug_set` (`GDSC2` default; also `GDSC1` / `CTRP2`)

`render_report` accepts an optional `format` field. With `format: "quarto"`
the task generates a Quarto (`.qmd`) document and renders it to both HTML
and PDF (override with `output_formats: ["html", "pdf", "docx"]`). When
Quarto is not installed the task falls back to the self-contained HTML
bundle.

## Outputs

Every run writes to:

- `tcga_runs/<task_id>/results/`
- `tcga_runs/<task_id>/plots/`
- `tcga_runs/<task_id>/objects/`
- `tcga_runs/<task_id>/report.md`
- `tcga_runs/<task_id>/run_metadata.json` (includes `toolkit_version`)

## Companion scripts (introspection, no analysis side effects)

- `list_tasks.R` — print every task and its required/optional fields
- `inspect_project.R --project TCGA-XXX` — summarize sample counts, group
  levels, available clinical columns and subtypes for a project
- `validate_config.R --config X.json` — dry-run schema validation against
  this spec
- `list_runs.R [--status failed] [--task run_deg] [--since 2026-05-01]` —
  summarize existing `tcga_runs/` entries
- `check_deps.R [--task run_deg]` — check R / Bioconductor / GitHub
  dependencies and emit install commands
