# Private research repo

This is your **personal** working copy of the TCGA project. It is intended
to live in a **private** GitHub repository alongside your public
[`tcga-toolkit`](https://github.com/<YOUR-HANDLE>/tcga-toolkit). Keep it
private — it contains research-direction breadcrumbs (specific gene
targets, pathway hypotheses, in-progress run outputs) that should not be
public until any associated paper is published.

## What this repo contains

| Path | Purpose | Committed? |
| --- | --- | --- |
| `tcga_toolkit/` | Git submodule pointing to your public toolkit | yes (as submodule) |
| `templates/` *(local copies)* | Your research-specific configs (mybl2, FA, hsa03460, etc) | yes |
| `tcga_runs/` | Result directories from `run_task.R` | yes (CSVs / PNGs / JSON only via `.gitignore` exceptions) |
| `notebooks/` | Legacy ipynb references | yes |
| `2-Output/` | Pre-toolkit analysis outputs | no (gitignored — uncomment in `.gitignore` to include) |
| `AGENTS.md`, `CLAUDE.md`, `.claude/` | Per-machine assistant configuration | partial (settings.json yes, settings.local.json no) |
| `0-Data/`, `GDCdata/`, `GTEX/`, `1-Input/Homo-*.txt` | Bulk data | **no** — backed up separately |
| `Human DNA Repair Genes.xlsx` | Third-party gene list | **no** — proprietary |

## Data backup (handled outside git)

The 21 GB of TCGA / GTEx / GDC data **never** goes into git, even in a
private repo. Recommended backup options:

- **Time Machine** — automatic, incremental, local.
- **External SSD / HDD** — quarterly snapshot.
- **Cloud archive** (Backblaze B2, S3 Glacier, university NAS) — yearly.

To re-create the working tree on a new machine after cloning this repo:

1. Clone with submodules:

   ```bash
   git clone --recurse-submodules https://github.com/<YOUR-HANDLE>/tcga-research.git
   ```

2. Restore the data directories from your backup (`0-Data/`, `GDCdata/`,
   `GTEX/`, `1-Input/`).

3. Verify:

   ```bash
   Rscript tcga_toolkit/scripts/audit_data.R   # not exposed — use:
   Rscript tcga_toolkit/scripts/run_task.R \
     --config tcga_toolkit/templates/audit_tcga.json
   bash tcga_toolkit/tests/run_all.sh
   ```

## Day-to-day workflow

- **Toolkit changes** (anything inside `tcga_toolkit/`): switch into the
  submodule directory, commit there, push to the public repo, then back
  in this repo `git add tcga_toolkit && git commit -m "bump toolkit"`.
- **Config / template changes** (research-specific): commit here.
- **Run outputs**: tracked CSV / PNG / JSON inside `tcga_runs/` are
  committed (see `.gitignore` exceptions). Large `*.rda` / `*.rds`
  objects are blocked by the global rules.

## What never to commit, even private

- TCGA / GTEx `*_mrna.rda`, `*_clinical.rda`, `*_maf.rda`. (Re-distributing
  TCGA controlled-access data is a TOS violation; the open-access subset
  is technically OK but still bulky and unnecessary in git.)
- `1-Input/Homo-*-gmts.txt` if they include KEGG / Reactome content —
  those carry their own redistribution restrictions.
- `Human DNA Repair Genes.xlsx` — proprietary curation.
