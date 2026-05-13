# tcga-toolkit

A config-driven analysis toolkit for **TCGA / TARGET / GTEx** data in R, with a
companion Python adapter for downstream survival modelling and SHAP
explainability.

The toolkit exposes **27 reusable tasks** (DEG, enrichment, ssGSEA, TME
deconvolution, GTEx comparison, pan-cancer survival maps, Lasso-Cox prognostic
models, WGCNA, CNV / methylation / immune phenotype, drug response,
external-cohort validation, etc.) behind a single JSON-driven dispatcher.

> **Heads-up**: This repository ships only the toolkit. It does **not** ship
> TCGA / GTEx data — you must obtain those yourself (see [Data setup](#data-setup)).

---

## Quick start

```bash
# 1. Clone
git clone https://github.com/yelu61/TCGA.git
cd TCGA

# 2. List available tasks
Rscript tcga_toolkit/scripts/list_tasks.R

# 3. Check R / Bioconductor / GitHub dependencies
Rscript tcga_toolkit/scripts/check_deps.R

# 4. Inspect a project (after you have data in place)
Rscript tcga_toolkit/scripts/inspect_project.R --project TCGA-BRCA

# 5. Validate a config (dry-run, no analysis)
Rscript tcga_toolkit/scripts/validate_config.R \
  --config tcga_toolkit/templates/brca_cohort_qc.json

# 6. Run a task
Rscript tcga_toolkit/scripts/run_task.R \
  --config tcga_toolkit/templates/brca_cohort_qc.json
```

Outputs land in `tcga_runs/<timestamp-task_id>/{results,plots,objects}/` plus a
`report.md` and `run_metadata.json` (containing `toolkit_version`).

---

## The 27 tasks

Run `Rscript tcga_toolkit/scripts/list_tasks.R` for the full schema. Highlights:

| Category | Tasks |
| --- | --- |
| **Data audit / prep** | `audit_data`, `prepare_bulk_rna`, `cohort_qc` |
| **Differential expression** | `run_deg`, `run_enrichment`, `ssgsea_score` |
| **Clinical / survival** | `run_clinical_assoc`, `survival_map`, `mutation_survival`, `tmb_analysis` |
| **Cross-cohort comparisons** | `run_gtex_compare`, `pan_cancer_expression`, `external_validate` |
| **Subtype / stage** | `subtype_analysis`, `stage_analysis` |
| **Mutation / CNV / methylation** | `maf_summary`, `cnv_summary`, `methylation_diff` |
| **Modelling** | `prognostic_model` (Lasso-Cox / RSF), `wgcna_modules`, `validate_sc_signatures` |
| **TME / immune** | `run_tme`, `immune_phenotype` |
| **Drug response** | `drug_response` (oncoPredict / GDSC / CTRP2) |
| **Visualisation / reports** | `gene_correlation_heatmap`, `render_report` (HTML / Quarto-PDF) |
| **Orchestration** | `pipeline` (multi-step with variable substitution) |

See [`tcga_toolkit/references/task_index.md`](tcga_toolkit/references/task_index.md)
for the long-form documentation of each task.

---

## Architecture

```
tcga_toolkit/
├── scripts/
│   ├── run_task.R                     # main dispatcher
│   ├── common.R                       # shared helpers
│   ├── task_*.R                       # 27 task implementations
│   ├── list_tasks.R                   # introspection: tasks
│   ├── inspect_project.R              # introspection: project metadata
│   ├── validate_config.R              # introspection: config dry-run
│   ├── list_runs.R                    # introspection: tcga_runs/
│   ├── check_deps.R                   # introspection: R / Bioc / GH deps
│   ├── cleanup_runs.R                 # housekeeping
│   └── enhance_*.R                    # post-processing utilities
├── templates/                         # ready-to-copy JSON configs
├── specs/README.md                    # one-line schema per task
├── references/                        # long-form docs
│   ├── task_index.md
│   ├── pipeline_recipes.md
│   ├── data_inventory.md
│   ├── troubleshooting.md
│   └── extension_guide.md
├── adapters/python/                   # Python downstream adapter
├── tests/                             # smoke tests
├── VERSION                            # 0.3.0
└── CHANGELOG.md
```

---

## Data setup

The toolkit expects data at fixed locations relative to the project root:

```
0-Data/          # one *_mrna.rda, *_clinical.rda, *_maf.rda per project
1-Input/         # MSigDB GMT files for enrichment / ssGSEA
GDCdata/         # optional raw GDC downloads
GTEX/            # GTEx matrix + tissue mapping (for run_gtex_compare)
```

Minimal download recipe (adapt to your needs):

```r
# TCGA-BRCA bulk RNA + clinical
library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")
GDCdownload(query)
data <- GDCprepare(query)
saveRDS(data, "0-Data/TCGA-BRCA_mrna.rda")

# MSigDB Hallmark GMT (download from MSigDB; CC-BY 4.0)
# https://www.gsea-msigdb.org/gsea/msigdb/
```

See [`tcga_toolkit/references/data_inventory.md`](tcga_toolkit/references/data_inventory.md)
for the full data layout the tasks expect, including optional CNV /
methylation matrices and the OSF link for oncoPredict training data.

---

## Smoke tests

```bash
bash tcga_toolkit/tests/run_all.sh
```

Verifies that all task scripts source cleanly, every advertised task has a
function defined, every template passes schema validation, and the Python
adapter imports.

---

## Python adapter

After running an analysis you can pick up the outputs from Python:

```python
import sys; sys.path.insert(0, "tcga_toolkit/adapters")
from python import load_run, load_risk_scores, fit_external_cox, time_dependent_auc

run = load_run("tcga_runs/20260513-XXXX-brca-demo-lasso-cox")
risk = load_risk_scores(run)
print(fit_external_cox(risk[risk["split"] == "test"]))
```

See [`tcga_toolkit/adapters/python/README.md`](tcga_toolkit/adapters/python/README.md)
for scikit-survival + SHAP integration.

---

## Adding a new task

See [`tcga_toolkit/references/extension_guide.md`](tcga_toolkit/references/extension_guide.md).
The short version:

1. Create `tcga_toolkit/scripts/task_<name>.R`.
2. Source it in `run_task.R` and add to the dispatch list.
3. Add a schema entry to `validate_config()` in `common.R`.
4. Document in `list_tasks.R`, `check_deps.R`, `task_index.md`,
   `specs/README.md`, drop a template in `templates/`, append a CHANGELOG entry.
5. Run the smoke tests.

---

## Citation

If you use this toolkit in published work, please cite:

```
<TBD — fill in once a preprint / paper / Zenodo DOI is available>
```

The underlying data sources should also be cited:

- **TCGA**: Weinstein, J. N. et al. *The Cancer Genome Atlas Pan-Cancer
  analysis project.* Nat Genet 45, 1113–1120 (2013).
- **GTEx**: GTEx Consortium. *The GTEx Consortium atlas of genetic regulatory
  effects across human tissues.* Science 369, 1318–1330 (2020).
- **MSigDB**: Liberzon, A. et al. *Molecular signatures database (MSigDB)
  3.0.* Bioinformatics 27, 1739–1740 (2011).
- **maftools** (Mayakonda et al. 2018), **DESeq2** (Love et al. 2014),
  **clusterProfiler** (Wu et al. 2021), **WGCNA** (Langfelder & Horvath
  2008), **IOBR** (Zeng et al. 2021), **oncoPredict** (Maeser et al. 2021),
  **GSVA** (Hänzelmann et al. 2013), **survival** (Therneau), **glmnet**
  (Friedman et al. 2010), **randomForestSRC** (Ishwaran et al.), **timeROC**
  (Blanche et al.).

---

## License

Apache License 2.0 — see [`LICENSE`](LICENSE).

This license covers the **code** in this repository. It does **not** cover
any TCGA / GTEx / MSigDB / GDC / oncoPredict data that users download
separately; those data are subject to their own terms of use.

---

## Contributing

Issues and pull requests welcome. See [`CONTRIBUTING.md`](CONTRIBUTING.md).
