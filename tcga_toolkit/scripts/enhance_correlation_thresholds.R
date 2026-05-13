#!/usr/bin/env Rscript

# Multi-threshold filtering for gene correlation results
# Usage: Rscript enhance_correlation_thresholds.R <run_dir>

args <- commandArgs(trailingOnly = TRUE)
run_dir <- args[1]

if (length(args) < 1) {
  stop("Usage: Rscript enhance_correlation_thresholds.R <run_dir>", call. = FALSE)
}

results_dir <- file.path(run_dir, "results")
plots_dir <- file.path(run_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

long_df <- utils::read.csv(file.path(results_dir, "MYBL2_gene_correlation_long.csv"), stringsAsFactors = FALSE)

# Define threshold combinations
cor_thresholds <- c(0.2, 0.3, 0.4, 0.5)
p_thresholds <- c(0.05, 0.01, 0.001)
min_project_counts <- c(1, 5, 10, 15, 20, 25, 30)

threshold_rows <- list()
for (cor_thr in cor_thresholds) {
  for (p_thr in p_thresholds) {
    for (min_proj in min_project_counts) {
      # Count genes meeting criteria
      gene_counts <- by(long_df, long_df$gene, function(df) {
        sum(abs(df$correlation) >= cor_thr & df$p_value < p_thr, na.rm = TRUE)
      })
      n_genes <- sum(gene_counts >= min_proj)
      total_genes <- length(unique(long_df$gene))

      threshold_rows[[length(threshold_rows) + 1L]] <- data.frame(
        cor_threshold = cor_thr,
        p_threshold = p_thr,
        min_projects = min_proj,
        genes_passed = as.integer(n_genes),
        total_genes = as.integer(total_genes),
        stringsAsFactors = FALSE
      )
    }
  }
}

threshold_df <- do.call(rbind, threshold_rows)
utils::write.csv(threshold_df, file.path(results_dir, "MYBL2_threshold_comparison.csv"), row.names = FALSE)

# Generate detailed gene lists for key thresholds
key_thresholds <- list(
  list(cor = 0.3, p = 0.05, min_proj = 5, label = "r0.3_p0.05_min5"),
  list(cor = 0.3, p = 0.05, min_proj = 10, label = "r0.3_p0.05_min10"),
  list(cor = 0.3, p = 0.05, min_proj = 15, label = "r0.3_p0.05_min15"),
  list(cor = 0.3, p = 0.05, min_proj = 20, label = "r0.3_p0.05_min20"),
  list(cor = 0.3, p = 0.01, min_proj = 10, label = "r0.3_p0.01_min10"),
  list(cor = 0.4, p = 0.05, min_proj = 10, label = "r0.4_p0.05_min10"),
  list(cor = 0.5, p = 0.05, min_proj = 10, label = "r0.5_p0.05_min10")
)

gene_list_rows <- list()
for (kt in key_thresholds) {
  gene_counts <- by(long_df, long_df$gene, function(df) {
    sum(abs(df$correlation) >= kt$cor & df$p_value < kt$p, na.rm = TRUE)
  })
  passed_genes <- names(gene_counts)[gene_counts >= kt$min_proj]
  for (g in passed_genes) {
    df_gene <- long_df[long_df$gene == g, ]
    sig_count <- sum(abs(df_gene$correlation) >= kt$cor & df_gene$p_value < kt$p, na.rm = TRUE)
    mean_r <- mean(df_gene$correlation, na.rm = TRUE)
    gene_list_rows[[length(gene_list_rows) + 1L]] <- data.frame(
      threshold_label = kt$label,
      gene = g,
      significant_projects = as.integer(sig_count),
      mean_correlation = mean_r,
      stringsAsFactors = FALSE
    )
  }
}

gene_list_df <- do.call(rbind, gene_list_rows)
utils::write.csv(gene_list_df, file.path(results_dir, "MYBL2_threshold_gene_lists.csv"), row.names = FALSE)

# Plot: threshold comparison bar chart
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_df <- threshold_df[threshold_df$p_threshold == 0.05, ]
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = factor(min_projects), y = genes_passed, fill = factor(cor_threshold))) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::labs(
      title = "MYBL2 vs DNA Repair Genes: Multi-threshold Screening",
      subtitle = "p < 0.05",
      x = "Minimum significant cancer types",
      y = "Number of genes passed",
      fill = "|r| threshold"
    ) +
    ggplot2::scale_fill_brewer(palette = "Blues") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )

  pdf_path <- file.path(plots_dir, "MYBL2_threshold_comparison.pdf")
  png_path <- file.path(plots_dir, "MYBL2_threshold_comparison.png")
  grDevices::pdf(pdf_path, width = 10, height = 6)
  print(p)
  grDevices::dev.off()
  grDevices::png(png_path, width = 1600, height = 960, res = 150)
  print(p)
  grDevices::dev.off()
}

# Summary table for key thresholds
summary_rows <- list()
for (kt in key_thresholds) {
  n <- sum(threshold_df$cor_threshold == kt$cor & threshold_df$p_threshold == kt$p & threshold_df$min_projects == kt$min_proj)
  if (n > 0) {
    row <- threshold_df[threshold_df$cor_threshold == kt$cor & threshold_df$p_threshold == kt$p & threshold_df$min_projects == kt$min_proj, ]
    summary_rows[[length(summary_rows) + 1L]] <- row
  }
}
summary_key <- do.call(rbind, summary_rows)
utils::write.csv(summary_key, file.path(results_dir, "MYBL2_key_threshold_summary.csv"), row.names = FALSE)

message("Threshold comparison complete.")
message(sprintf("Total genes in pool: %s", length(unique(long_df$gene))))
message("Key results:")
for (i in seq_len(nrow(summary_key))) {
  message(sprintf("  |r| >= %.1f, p < %.3f, >= %s projects -> %s genes",
                  summary_key$cor_threshold[i], summary_key$p_threshold[i], summary_key$min_projects[i], summary_key$genes_passed[i]))
}
