# Data Inventory

A map of the static datasets this toolkit assumes are present in the
repository. Paths are relative to the project root (the directory that
contains `tcga_toolkit/`).

Use `audit_data` (or `inspect_project.R --list`) to programmatically
verify presence and freshness.

---

## `0-Data/` — Per-project bundles

For each TCGA / TARGET project the toolkit expects up to four files:

| File | Type | Used by |
| --- | --- | --- |
| `<PROJECT>_mrna.rda` | `RangedSummarizedExperiment` (counts + tpm + clinical) | every task that needs expression |
| `<PROJECT>_clinical.rda` | data.frame | sometimes loaded separately for richer clinical fields |
| `<PROJECT>_maf.rda` | maftools-compatible MAF | `maf_summary`, `mutation_survival`, `tmb_analysis` |
| `1-Input/<PROJECT>_DEGinput.rda` | optional pre-prepared bundle (`Counts`, `TPM`, `clinical_clean`, `id2symbol`) | shortcut for `prepare_bulk_rna` / `run_deg` |

Pan-cancer artefacts in `0-Data/`:

- `PanCancer_subtypes.rda` — used by `subtype_analysis`
  (column `Subtype_Selected`, key `pan.samplesID`)
- `PanCancer_TMB.Rdata` — preferred source for `tmb_analysis`
  (`tmb_source: auto` checks here first)

Use `inspect_project.R --project TCGA-BRCA` to see which of these files
exist for a given cohort and what clinical columns / subtype levels are
available before authoring a config.

## `1-Input/` — Curated inputs

| Path | Purpose |
| --- | --- |
| `Homo-HALLMARK-gmts.txt` | MSigDB Hallmark, used by `run_enrichment` and `ssgsea_score` by default |
| `Homo-GOCC-gmts.txt` | Gene Ontology Cellular Component |
| `Homo-GOMF-gmts.txt` | Gene Ontology Molecular Function |
| `<PROJECT>_DEGinput.rda` | optional precomputed bundle |
| user-provided signatures / GMTs | safe to drop here |

## `GDCdata/` — Raw GDC downloads

One subdirectory per project (`GDCdata/TCGA-BRCA/`). These are the
unprocessed files retrieved with TCGAbiolinks; the toolkit does not read
them directly but `audit_data` reports their presence and last-modified
times.

## `GTEX/` — GTEx reference panel

| File | Used for |
| --- | --- |
| `gtex_RSEM_gene_tpm.gz` | streamed per-gene extraction in `run_gtex_compare` (the toolkit `zgrep`s this file rather than loading it) |
| `tissue_gtex.csv` | sample → tissue mapping |
| `samplepair.txt` | TCGA-project → GTEx-tissue lookup |
| `gencode.v23.annotation.gene.probemap` | Ensembl ID ↔ symbol bridge |
| `exp_gtex.tpm.csv` | the full 8.5 GB matrix; **do not load** unless absolutely needed |

## Repository root extras

- `Human DNA Repair Genes.xlsx` — referenced from the
  `gene_correlation_heatmap` template and recipe R7.
- `tcga_runs/<task_id>/` — output of every run.
- `.tcga_cache/<project>_prepared.rds` — auto-rebuilt cache; safe to
  delete to force a fresh load.

## Optional / user-provided datasets

These are not shipped with the repository but several tasks accept them
when present:

| Location | Used by | How to get it |
| --- | --- | --- |
| `0-Data/<PROJECT>_cnv.rda` (or `_gistic.rda` / `.tsv` / `.csv`) | `cnv_summary` | `TCGAbiolinks::GDCquery(project='TCGA-XXX', data.category='Copy Number Variation', data.type='Gene Level Copy Number')` then save as a gene × sample matrix |
| `0-Data/<PROJECT>_methylation.rda` (or `.tsv` / `.csv`) | `methylation_diff` | `TCGAbiolinks::GDCquery(project='TCGA-XXX', data.category='DNA Methylation', data.type='Methylation Beta Value', platform='Illumina Human Methylation 450')` |
| `1-Input/external/<COHORT>_expression.tsv` + `<COHORT>_clinical.tsv` | `external_validate` | Manually download CGGA / METABRIC / GEO series; clinical needs survival time + event columns |
| `1-Input/oncoPredict_training_data/<set>_Expr.rds` + `<set>_Res.rds` | `drug_response` | Download GDSC1 / GDSC2 / CTRP2 training data from <https://osf.io/c6tfx/> |

## Pre-flight checklist

Before composing a new config:

1. `Rscript tcga_toolkit/scripts/inspect_project.R --project <PROJECT>` to
   confirm which group / subtype / clinical columns are available.
2. `Rscript tcga_toolkit/scripts/check_deps.R --task <TASK>` to confirm
   the dependencies for the task you want.
3. `Rscript tcga_toolkit/scripts/list_runs.R --task <TASK>` to see whether
   the same analysis has already been run.
