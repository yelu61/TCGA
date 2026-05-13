# TCGA Toolkit References

The toolkit assumes the repository keeps these source directories stable:

- `0-Data/`: serialized project-level expression, clinical, and MAF objects
- `1-Input/`: local GMT collections and curated intermediate inputs
- `GDCdata/`: raw GDC downloads by project
- `GTEX/`: GTEx matrix, tissue annotation, and TCGA-GTEx tissue mapping

Current local assumptions encoded in the scripts:

- `*_mrna.rda` is a `RangedSummarizedExperiment`
- preferred assays are `unstranded` for counts and `tpm_unstrand` for TPM
- GTEx expression comes from `GTEX/gtex_RSEM_gene_tpm.gz`
- GTEx tissue matching defaults to `GTEX/samplepair.txt`

If a future refresh changes these formats, update
[`tcga_toolkit/scripts/common.R`](../scripts/common.R) first.

---

## Topic-specific docs

Read these by topic when planning or debugging a task — start with the
short SKILL.md, then drill into whichever topic is closest to your
problem.

| File | When to read |
| --- | --- |
| [task_index.md](task_index.md) | Full schema, inputs, outputs and gotchas for every task. |
| [pipeline_recipes.md](pipeline_recipes.md) | Cookbook of common multi-step pipelines (DEG → enrichment, gene deep-dive, mutation landscape, etc). |
| [data_inventory.md](data_inventory.md) | Where the input rda / GMT / GTEx files live and what they contain. |
| [troubleshooting.md](troubleshooting.md) | Common failure modes and fixes. |
| [extension_guide.md](extension_guide.md) | How to add a new task without breaking the rest. |

## Introspection scripts

Lightweight scripts that print machine-readable answers about the toolkit
state without running an analysis:

- `Rscript tcga_toolkit/scripts/list_tasks.R [--task X] [--json]`
- `Rscript tcga_toolkit/scripts/inspect_project.R --project TCGA-XXX [--json]`
- `Rscript tcga_toolkit/scripts/validate_config.R --config X.json [--json]`
- `Rscript tcga_toolkit/scripts/list_runs.R [--status failed] [--task X] [--since YYYY-MM-DD] [--json]`
- `Rscript tcga_toolkit/scripts/check_deps.R [--task X] [--json]`

## Smoke tests

```bash
bash tcga_toolkit/tests/run_all.sh
```

Verifies that all task scripts source, every advertised task has a
function defined, every template passes schema validation, and the
Python adapter imports cleanly. See [`../tests/README.md`](../tests/README.md).

## Python adapter

Read the R toolkit's `tcga_runs/` outputs from Python and run downstream
scikit-survival / SHAP analyses:

```bash
python tcga_toolkit/adapters/python/...  # see adapters/python/README.md
```
