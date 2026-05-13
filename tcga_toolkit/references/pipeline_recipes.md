# Pipeline Recipes

Reusable cookbook for multi-step analyses. Each recipe is a fragment that
can be lifted into a JSON config file passed to `run_task.R`. All recipes
use `pipeline` so later steps consume earlier outputs via
`{{stepN.glob:*.csv}}` or `{{stepN.run_dir}}`.

For the underlying field schema, see [task_index.md](task_index.md) and
[`tcga_toolkit/specs/README.md`](../specs/README.md).

---

## R1. DEG → Enrichment → Report

End-to-end discovery pipeline for one project.

```json
{
  "task": "pipeline",
  "task_id": "brca-deg-enrich-report",
  "steps": [
    {
      "task": "run_deg",
      "project": "TCGA-BRCA",
      "numerator": "Tumor",
      "denominator": "Normal"
    },
    {
      "task": "run_enrichment",
      "input_deg": "{{step0.glob:*_deg.csv}}",
      "mode": "ora"
    },
    {
      "task": "render_report",
      "run_dirs": ["{{step0.run_dir}}", "{{step1.run_dir}}"]
    }
  ]
}
```

## R2. Single-cell signature → bulk validation → pan-cancer survival

Validate a signature on one cohort, then sweep it across all projects.

```json
{
  "task": "pipeline",
  "task_id": "sc-validate-then-pancancer",
  "steps": [
    {
      "task": "validate_sc_signatures",
      "project": "TCGA-BRCA",
      "signature_file": "1-Input/my_sc_signature.tsv"
    },
    {
      "task": "survival_map",
      "signature_file": "1-Input/my_sc_signature.tsv"
    },
    {
      "task": "render_report",
      "run_dirs": ["{{step0.run_dir}}", "{{step1.run_dir}}"]
    }
  ]
}
```

## R3. Cohort QC → DEG → Prognostic model

Always QC first when you plan to build a model.

```json
{
  "task": "pipeline",
  "task_id": "brca-qc-deg-prognostic",
  "steps": [
    {
      "task": "cohort_qc",
      "project": "TCGA-BRCA"
    },
    {
      "task": "run_deg",
      "project": "TCGA-BRCA",
      "numerator": "Tumor",
      "denominator": "Normal",
      "padj_cutoff": 0.01,
      "log2fc_cutoff": 1.5
    },
    {
      "task": "prognostic_model",
      "project": "TCGA-BRCA",
      "method": "lasso_cox",
      "signature_file": "{{step1.glob:*_deg.csv}}"
    }
  ]
}
```

Note: `signature_file` for `prognostic_model` accepts a marker table; the
DEG CSV already has a `gene` column, so it works directly.

## R4. Gene-of-interest deep dive

Comprehensive characterisation of one gene (e.g. an oncogene candidate).

```json
{
  "task": "pipeline",
  "task_id": "tp53-deep-dive",
  "steps": [
    { "task": "pan_cancer_expression", "gene": "TP53" },
    { "task": "survival_map", "gene": "TP53" },
    { "task": "run_clinical_assoc", "project": "TCGA-BRCA", "gene": "TP53" },
    { "task": "subtype_analysis", "project": "TCGA-BRCA", "gene": "TP53" },
    { "task": "stage_analysis", "project": "TCGA-BRCA", "gene": "TP53" },
    { "task": "run_gtex_compare", "project": "TCGA-BRCA", "gene": "TP53" },
    {
      "task": "render_report",
      "run_dirs": [
        "{{step0.run_dir}}", "{{step1.run_dir}}", "{{step2.run_dir}}",
        "{{step3.run_dir}}", "{{step4.run_dir}}", "{{step5.run_dir}}"
      ]
    }
  ]
}
```

## R5. Mutation landscape + survival impact

For driver-gene candidates.

```json
{
  "task": "pipeline",
  "task_id": "brca-mutation-landscape",
  "steps": [
    { "task": "maf_summary", "project": "TCGA-BRCA", "top_n": 25 },
    { "task": "tmb_analysis", "project": "TCGA-BRCA" },
    { "task": "mutation_survival", "project": "TCGA-BRCA", "gene": "TP53" },
    {
      "task": "render_report",
      "run_dirs": [
        "{{step0.run_dir}}", "{{step1.run_dir}}", "{{step2.run_dir}}"
      ]
    }
  ]
}
```

## R6. Pathway-level pan-cancer survey

ssGSEA per project, then export.

```json
{
  "task": "pipeline",
  "task_id": "brca-pathway-survey",
  "steps": [
    {
      "task": "ssgsea_score",
      "project": "TCGA-BRCA",
      "gmt_paths": ["1-Input/Homo-HALLMARK-gmts.txt"]
    },
    {
      "task": "render_report",
      "run_dirs": ["{{step0.run_dir}}"]
    }
  ]
}
```

You can then use the resulting per-sample pathway matrix as input for
`prognostic_model` by treating one pathway column as a covariate.

## R7. Gene-set correlation map across cohorts

```json
{
  "task": "gene_correlation_heatmap",
  "task_id": "mybl2-dna-repair-correlation",
  "target_gene": "MYBL2",
  "gene_list_file": "Human DNA Repair Genes.xlsx",
  "sheet": "Sheet1",
  "method": "spearman",
  "projects": ["TCGA-BRCA", "TCGA-LUAD", "TCGA-COAD"],
  "sample_filter": ["Tumor"]
}
```

## R8. Lasso-Cox in TCGA → external validation in CGGA → Quarto report

End-to-end "model in TCGA, validate elsewhere, hand the PDF to a
collaborator" pipeline.

```json
{
  "task": "pipeline",
  "task_id": "brca-prognostic-with-external",
  "steps": [
    {
      "task": "cohort_qc",
      "project": "TCGA-BRCA"
    },
    {
      "task": "prognostic_model",
      "project": "TCGA-BRCA",
      "method": "lasso_cox",
      "signature_file": "tcga_toolkit/templates/fa_kegg_hsa03460_genes.txt"
    },
    {
      "task": "external_validate",
      "cohort_name": "CGGA",
      "expression_file": "1-Input/external/CGGA_expression.tsv",
      "clinical_file": "1-Input/external/CGGA_clinical.tsv",
      "weight_file": "{{step1.glob:*_prognostic_coefficients.csv}}"
    },
    {
      "task": "render_report",
      "format": "quarto",
      "output_formats": ["html", "pdf"],
      "title": "Prognostic model: TCGA discovery + CGGA validation",
      "run_dirs": [
        "{{step0.run_dir}}", "{{step1.run_dir}}", "{{step2.run_dir}}"
      ]
    }
  ]
}
```

## R9. CNV + methylation + RNA multi-omic snapshot

When the user supplies CNV and methylation matrices alongside RNA, this
gives a one-shot multi-omic landscape.

```json
{
  "task": "pipeline",
  "task_id": "brca-multiomic-snapshot",
  "steps": [
    { "task": "run_deg", "project": "TCGA-BRCA",
      "numerator": "Tumor", "denominator": "Normal" },
    { "task": "cnv_summary", "project": "TCGA-BRCA" },
    { "task": "methylation_diff", "project": "TCGA-BRCA",
      "numerator": "Tumor", "denominator": "Normal" },
    { "task": "render_report",
      "format": "quarto",
      "run_dirs": [
        "{{step0.run_dir}}", "{{step1.run_dir}}", "{{step2.run_dir}}"
      ]
    }
  ]
}
```

## R10. WGCNA → immune phenotype → drug response

For "find a co-expression module → see if it's immune-related → look for
therapeutic angles" workflows.

```json
{
  "task": "pipeline",
  "task_id": "brca-wgcna-immune-drugs",
  "steps": [
    { "task": "wgcna_modules", "project": "TCGA-BRCA",
      "sample_filter": ["Tumor"] },
    { "task": "immune_phenotype", "project": "TCGA-BRCA" },
    { "task": "drug_response", "project": "TCGA-BRCA",
      "drug_set": "GDSC2" }
  ]
}
```
