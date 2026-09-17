# TCGA Toolkit Instructions

This repository uses a shared, agent-neutral execution model for TCGA analysis. Codex should treat `tcga_toolkit/` as the canonical entrypoint for reusable workflow logic.

## Default Workflow

1. Inspect `tcga_toolkit/specs/README.md`.
2. Reuse or create a config file under `tcga_toolkit/templates/` or a nearby working path.
3. Run:

```bash
Rscript tcga_toolkit/scripts/run_task.R --config <config.json> [--output-root <dir>]
```

   - `--output-root`: write `tcga_runs/` to another project instead of the TCGA repo.
   - Native run bundles are immutable. `--overwrite` is rejected; directory collisions fail before execution.
   - Config snapshots, declared inputs, output hashes and the R session are recorded for each run, including pipeline steps. See `tcga_toolkit/references/result_contract.md` for scope and limitations.

4. Read results from `tcga_runs/<task_id>/`.

5. Preview archive candidates in `tcga_runs/` (filters combine with AND):

```bash
Rscript tcga_toolkit/scripts/cleanup_runs.R --keep 20
Rscript tcga_toolkit/scripts/cleanup_runs.R --older-than 30 --failed-only --dry-run
```

   Cleanup defaults to preview. Applying an archive requires current/released run protection and fresh project lifecycle evidence; the archive is reversible and never deletes native run files. See `tcga_toolkit/references/result_contract.md`.

## Expectations

- Do not duplicate stable analysis logic inside notebooks when a toolkit script can own it.
- Keep public-data preparation and cohort/multi-omics tasks here; generic matrix analyses can be handed to RNAseq-Templates with explicit scale, gene IDs, sample metadata and processing history.
- `bulk-rnaseq-analysis` routes scientific work and adapts native results. `manage-research-project` owns generic artifact registration, current selections, freely named collections and frozen deliveries. Panel IDs are optional manuscript layout information, not a required output structure.
- Prefer updating `tcga_toolkit/scripts/` over creating one-off analysis fragments.
- Reuse the repository's local datasets in `0-Data/`, `1-Input/`, `GDCdata/`, and `GTEX/`.
- Assume `*_mrna.rda` inputs are `RangedSummarizedExperiment` objects unless inspection shows otherwise.
- For GTEx comparisons, prefer targeted extraction from `GTEX/gtex_RSEM_gene_tpm.gz` and tissue matching from `GTEX/samplepair.txt`.

## Supported Tasks

- `audit_data`
- `prepare_bulk_rna`
- `run_deg`
- `run_enrichment`
- `run_clinical_assoc`
- `run_tme`
- `run_gtex_compare`
- `validate_sc_signatures`
- `pan_cancer_expression`
- `stage_analysis`
- `survival_map`
- `subtype_analysis`
- `gene_correlation_heatmap`
- `pipeline` (multi-step orchestration with variable substitution)
- `render_report` (generate a unified HTML report from one or more `run_dirs`)
- `maf_summary` (mutation burden, oncoplot, gene frequency from MAF)
- `mutation_survival` (mutant vs wildtype survival for a given gene)
