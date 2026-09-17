# Native runs and project results

TCGA Toolkit owns public-data preparation, cohort-specific clinical/multi-omics
analysis and native execution evidence. RNAseq-Templates owns generic analyses
starting from expression matrices. The bulk RNA-seq skill selects and connects
the backends; the project-management skill owns project artifact registration,
current result selections, flexible collections and milestone releases.

A handoff must state matrix scale (counts, TPM, transformed values), gene-ID
namespace, sample-to-patient mapping, covariates and processing history. A file
path by itself does not authorize a statistical model. This result contract
does not replace task-specific scientific validation.

## Native schema v1

The supported CLI creates a new `tcga_runs/<task_id>/` directory. A collision or
`--overwrite` request fails before execution. The same applies to pipeline
children, whose config snapshots contain resolved step substitutions. Calling
individual historical task functions directly does not create this contract;
use `run_task.R` for registered executions.

`run_metadata.json` preserves the existing `task`, `task_id`, `toolkit_version`,
`config_path`, `run_dir`, `status`, `error` and `completed_at` fields and adds:

| Field | Meaning |
| --- | --- |
| `manifest_schema_version` | Integer `1`. |
| `backend` | `tcga-toolkit`. |
| `backend_revision`, `backend_dirty` | Git revision and worktree state; unavailable identity is explicitly recorded. |
| `config_file`, `config_md5` | `config_used.json` and its digest; effective config saved before execution. |
| `input_manifest_file`, `input_manifest_md5` | `run_inputs.tsv` and its digest. |
| `artifact_manifest_file`, `artifact_manifest_md5` | `run_artifacts.tsv` and its digest. |
| `environment_file` | `sessionInfo.txt`; R/package session at task completion. |
| `provenance_status` | `pending`, `complete`, or `failed`. |
| `provenance_scope` | Explicitly limited discovery scope described below. |
| `parent_run_id` | Pipeline child's parent, when applicable. |
| `dependency_run_ids` | Completed pipeline parent's child runs, when applicable. |

Execution status is `running`, `success` or `failed`. A provenance preparation or
finalization error forces `failed`, even if a handler returned normally. A
failed handler can still have complete provenance and partial output artifacts;
neither state makes those outputs eligible as accepted scientific results.

Both TSV files have columns `path`, `role`, `md5`, `bytes`, `status`. MD5 binds
native files; the generic project registry can independently use SHA-256.

- Input paths are absolute for local files. Config-relative paths are resolved
  against the invocation working directory, consistent with current task code.
  External dependencies are recorded by reference, without copying cohort data.
- Input `role` is `declared:<config-key>` or `cohort_candidate:<file-type>`.
  Declared path fields and `provenance_inputs` are traversed; explicit input
  directories, including report `run_dirs`, expand to individual files. Known
  existing expression/clinical/MAF/prepared files for explicitly named projects
  are inventoried as candidates, without claiming every candidate was consumed.
- Input status is `recorded`, `missing`, `unresolved` (URI or pipeline token),
  or `directory_empty`. Non-file rows have empty hashes and byte counts.
  Declared recorded inputs are checked again after execution for changes.
- Artifact paths are run-relative; every regular output file is hashed. Roles
  are `figure`, `table`, `report`, `object`, `provenance` or `other`. All artifact
  rows have status `recorded`. The inventory includes config/session files,
  and excludes metadata and both manifests to avoid self-reference. Symlinked
  output files/directories fail provenance finalization.

This inventory is **not a complete runtime I/O audit**. Implicit caches,
third-party reference packages, remotely fetched material and projects discovered
inside tasks may need additional explicit input declarations and scientific
review. Session information describes the environment; it does not lock or
restore it. Git dirty state identifies uncommitted execution without preserving
the patch itself. Missing evidence must remain explicit during registration.

## Flexible project results

Native runs retain their backend layout. Select individual artifacts into the
project's current results and freely named collections through the shared
project lifecycle interface. A collection may follow a research question,
comparison, method, audience or manuscript topic. Figure and table names remain
editable; panel IDs are optional layout metadata. Registration alone does not
change the current selection or create a release.

Sharing materializes the selected files and their manifest. A milestone release
freezes exact selections; later plot/layout updates can change current results
without modifying that release. Generated files and manually maintained files
must be distinguished so updates preserve manual work. Native run timestamps
do not determine which result is scientifically accepted.

## Archive preview and recovery

`cleanup_runs.R` defaults to preview, and supplied age, recency and failure
filters combine with AND. Unknown/running directories and symbolic links are
excluded. Applying requires `--archive-plan FILE` from the project lifecycle
helper. Only that plan's explicit `candidate_native_directories` authorize moves;
`unreferenced_registered_directories` remain review suggestions. The script
verifies the project root, SHA-256 hashes of registry/result files, and result
scope file listings to reject stale plans, including newly added manual files.
Unregistered native directories are never candidates for application.

Protected native directories come from that plan. Optional `--protected-runs`
adds one native run ID per line; it cannot authorize moves. Known report-input
and pipeline dependencies extend protection transitively. Incomplete legacy or
manual dependency evidence must be resolved by the lifecycle helper before it
emits an actionable candidate list. No automatic cleanup follows analysis.

```bash
# Preview with fresh project-generated lifecycle evidence:
Rscript tcga_toolkit/scripts/cleanup_runs.R --output-root /path/to/study \
  --older-than 30 --failed-only --archive-plan /path/to/archive-plan.json
# Explicitly apply that same verified selection:
Rscript tcga_toolkit/scripts/cleanup_runs.R --output-root /path/to/study \
  --older-than 30 --failed-only --archive-plan /path/to/archive-plan.json --apply
```

Applying archive moves eligible bundles under
`archive/tcga_runs/<batch>/<task_id>/`, retaining all bytes. Each batch records
`ARCHIVE_MANIFEST.tsv`, protected IDs and `RESTORE.txt`. A failed move leaves
earlier moves explicitly recorded and recoverable. Restore only into an absent
original path, with no overwrite. Archive does not certify scientific validity
or successful reproduction of historical results.
