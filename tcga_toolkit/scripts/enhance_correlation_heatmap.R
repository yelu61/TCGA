#!/usr/bin/env Rscript

# Post-processing script for gene correlation heatmap results.
# Reads an existing run directory produced by task_gene_correlation_heatmap
# and generates an enhanced heatmap with significance stars + summary tables.
# Optional: <hr_file> adds a forest plot as a ComplexHeatmap rowAnnotation
# (perfectly aligned with the heatmap rows).

args <- commandArgs(trailingOnly = TRUE)
run_dir <- args[1]
min_abs_cor <- as.numeric(args[2] %||% 0.3)
max_p <- as.numeric(args[3] %||% 0.05)
min_significant_projects <- as.integer(args[4] %||% 5)
show_dend <- as.logical(args[5] %||% "TRUE")
color_scheme <- args[6] %||% "blue_white_red"
sort_by <- args[7] %||% "mean_correlation"
highlight_gene <- args[8] %||% NULL
hr_file <- args[9] %||% NULL

if (length(args) < 1) {
  stop("Usage: Rscript enhance_correlation_heatmap.R <run_dir> [min_abs_cor] [max_p] [min_significant_projects] [show_dend] [color_scheme] [sort_by] [highlight_gene] [hr_file]", call. = FALSE)
}

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
})

results_dir <- file.path(run_dir, "results")
plots_dir <- file.path(run_dir, "plots")

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

long_df <- utils::read.csv(file.path(results_dir, "MYBL2_gene_correlation_long.csv"), stringsAsFactors = FALSE)
long_df$p_signif <- long_df$p_value < max_p & abs(long_df$correlation) >= min_abs_cor

# Summarize per gene
gene_summary <- do.call(rbind, by(long_df, long_df$gene, function(df) {
  data.frame(
    gene = unique(df$gene),
    mean_correlation = mean(df$correlation, na.rm = TRUE),
    median_correlation = stats::median(df$correlation, na.rm = TRUE),
    max_correlation = max(df$correlation, na.rm = TRUE),
    min_correlation = min(df$correlation, na.rm = TRUE),
    significant_projects = sum(df$p_value < max_p & abs(df$correlation) >= min_abs_cor, na.rm = TRUE),
    total_projects = nrow(df),
    stringsAsFactors = FALSE
  )
}))
gene_summary$composite_score <- 0.6 * gene_summary$mean_correlation + 0.4 * gene_summary$max_correlation
gene_summary <- switch(
  sort_by,
  composite_score = gene_summary[order(-gene_summary$composite_score, -gene_summary$max_correlation, -gene_summary$mean_correlation), ],
  max_correlation = gene_summary[order(-gene_summary$max_correlation, -gene_summary$mean_correlation), ],
  significant_projects = gene_summary[order(-gene_summary$significant_projects, -gene_summary$mean_correlation), ],
  min_correlation = gene_summary[order(-gene_summary$min_correlation, -gene_summary$mean_correlation), ],
  gene_summary[order(-gene_summary$mean_correlation), ]
)
gene_summary <- gene_summary[!is.na(gene_summary$gene), ]

selected_genes <- gene_summary$gene[gene_summary$significant_projects >= min_significant_projects]
if (!length(selected_genes)) {
  selected_genes <- head(gene_summary$gene, 30)
}

if (!is.null(highlight_gene) && !is.na(highlight_gene) && nzchar(highlight_gene) && highlight_gene != "NULL" && highlight_gene != "null") {
  selected_genes <- c(highlight_gene, setdiff(selected_genes, highlight_gene))
}

cluster_rows_flag <- isTRUE(sort_by == "cluster" || sort_by == "hclust")

# Prepare matrix for selected genes
df_sel <- long_df[long_df$gene %in% selected_genes, ]
df_sel$project_code <- factor(df_sel$project_code, levels = unique(df_sel$project_code[order(df_sel$project)]))
projects <- levels(df_sel$project_code)
genes <- selected_genes

cor_mat <- matrix(NA_real_, nrow = length(genes), ncol = length(projects), dimnames = list(genes, projects))
p_mat <- cor_mat
for (i in seq_len(nrow(df_sel))) {
  cor_mat[df_sel$gene[i], df_sel$project_code[i]] <- df_sel$correlation[i]
  p_mat[df_sel$gene[i], df_sel$project_code[i]] <- df_sel$p_value[i]
}

# Significance annotation layer
stars <- matrix("", nrow = nrow(cor_mat), ncol = ncol(cor_mat), dimnames = dimnames(cor_mat))
stars[!is.na(p_mat) & p_mat < 0.001 & abs(cor_mat) >= min_abs_cor] <- "***"
stars[!is.na(p_mat) & p_mat >= 0.001 & p_mat < 0.01 & abs(cor_mat) >= min_abs_cor] <- "**"
stars[!is.na(p_mat) & p_mat >= 0.01 & p_mat < max_p & abs(cor_mat) >= min_abs_cor] <- "*"

# Shared row annotation pieces
mean_cor_anno <- ComplexHeatmap::anno_barplot(
  gene_summary$mean_correlation[match(genes, gene_summary$gene)],
  axis = TRUE,
  gp = grid::gpar(fill = "#2c3e50"),
  width = grid::unit(1.6, "cm"),
  which = "row"
)
sig_projects_anno <- ComplexHeatmap::anno_barplot(
  gene_summary$significant_projects[match(genes, gene_summary$gene)],
  axis = TRUE,
  gp = grid::gpar(fill = "#e74c3c"),
  width = grid::unit(1.6, "cm"),
  which = "row"
)

colors <- switch(
  color_scheme,
  green_black_red = circlize::colorRamp2(c(-0.6, 0, 0.6), c("red", "black", "green")),
  purple_orange = circlize::colorRamp2(c(-0.6, 0, 0.6), c("#54278f", "#f7f7f7", "#f16913")),
  viridis_like = circlize::colorRamp2(c(-0.6, 0, 0.6), c("#440154", "#fde725", "#21918c")),
  rdylbu = circlize::colorRamp2(c(-0.6, -0.3, 0, 0.3, 0.6), c("#d73027", "#fc8d59", "#ffffbf", "#91bfdb", "#4575b4")),
  spectral = circlize::colorRamp2(c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6), c("#5e4fa2", "#3288bd", "#abdda4", "#ffffbf", "#fdae61", "#f46d43", "#d53e4f")),
  rich_diverging = circlize::colorRamp2(c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6), c("#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac")),
  circlize::colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red"))
)

has_hr <- !is.null(hr_file) && file.exists(hr_file)

# --- Build right_annotation (with or without forest) ---
if (has_hr) {
  hr_df <- utils::read.csv(hr_file, stringsAsFactors = FALSE)
  hr_df <- hr_df[match(genes, hr_df$gene), ]
  hr_df$gene <- genes
  hr_df$sig_text <- ""
  hr_df$sig_text[!is.na(hr_df$p) & hr_df$p < 0.05] <- "*"
  hr_df$sig_text[!is.na(hr_df$p) & hr_df$p < 0.01] <- "**"
  hr_df$sig_text[!is.na(hr_df$p) & hr_df$p < 0.001] <- "***"

  # Auto-adapt x-axis to data: symmetric around HR=1 in log10 space
  finite_vals <- c(hr_df$hr, hr_df$ci_low, hr_df$ci_high)
  finite_vals <- finite_vals[is.finite(finite_vals) & finite_vals > 0]
  if (length(finite_vals) >= 2) {
    log_vals <- log10(finite_vals)
    log_span <- max(abs(log_vals), na.rm = TRUE)
    log_span <- max(log_span, 0.15)
    x_min <- 10^(-log_span * 1.15)
    x_max <- 10^(log_span * 1.15)
  } else {
    x_min <- 0.5
    x_max <- 2.0
  }
  candidate_breaks <- c(0.1, 0.25, 0.33, 0.5, 0.67, 0.8, 1, 1.25, 1.5, 2, 3, 4, 5, 7, 10)
  axis_breaks <- candidate_breaks[candidate_breaks >= x_min * 0.95 & candidate_breaks <= x_max * 1.05]
  if (length(axis_breaks) < 3) axis_breaks <- c(x_min, 1, x_max)

  forest_anno <- ComplexHeatmap::AnnotationFunction(
    fun = function(index, k, n) {
      sub <- hr_df_local[index, , drop = FALSE]
      n_rows <- length(index)

      grid::pushViewport(grid::viewport(
        xscale = c(log10(x_min_local), log10(x_max_local)),
        yscale = c(0.5, n_rows + 0.5)
      ))

      # Faint background gridlines at axis breaks
      for (br in axis_breaks_local) {
        col_br <- if (abs(br - 1) < 1e-6) NA else "gray92"
        if (!is.na(col_br)) {
          grid::grid.lines(
            x = grid::unit(log10(br), "native"),
            y = grid::unit(c(0.5, n_rows + 0.5), "native"),
            gp = grid::gpar(col = col_br, lwd = 0.4)
          )
        }
      }
      # Reference line at HR = 1
      grid::grid.lines(
        x = grid::unit(0, "native"),
        y = grid::unit(c(0.5, n_rows + 0.5), "native"),
        gp = grid::gpar(col = "gray45", lty = "dashed", lwd = 0.7)
      )

      for (i in seq_len(nrow(sub))) {
        hr <- sub$hr[i]
        if (is.na(hr) || hr <= 0) next

        y_pos <- n_rows - i + 1
        is_sig <- isTRUE(sub$p[i] < 0.05)
        color <- if (is_sig) "#c0392b" else "#7f8c8d"

        ci_l_clip <- max(sub$ci_low[i], x_min_local, na.rm = TRUE)
        ci_h_clip <- min(sub$ci_high[i], x_max_local, na.rm = TRUE)

        grid::grid.segments(
          x0 = grid::unit(log10(ci_l_clip), "native"),
          x1 = grid::unit(log10(ci_h_clip), "native"),
          y0 = grid::unit(y_pos, "native"),
          y1 = grid::unit(y_pos, "native"),
          gp = grid::gpar(col = color, lwd = 1.5)
        )
        grid::grid.points(
          x = grid::unit(log10(hr), "native"),
          y = grid::unit(y_pos, "native"),
          pch = 16,
          size = grid::unit(1.7, "mm"),
          gp = grid::gpar(col = color)
        )

        sig_text <- sub$sig_text[i]
        if (nzchar(sig_text)) {
          grid::grid.text(
            sig_text,
            x = grid::unit(log10(ci_h_clip), "native") + grid::unit(0.4, "mm"),
            y = grid::unit(y_pos, "native"),
            just = c("left", "centre"),
            gp = grid::gpar(col = "#c0392b", fontface = "bold", fontsize = 8)
          )
        }
      }

      grid::popViewport()
    },
    which = "row",
    width = grid::unit(3.5, "cm"),
    var_import = list(
      hr_df_local = hr_df,
      x_min_local = x_min,
      x_max_local = x_max,
      axis_breaks_local = axis_breaks
    )
  )

  hr_label <- vapply(hr_df$hr, function(h) if (is.na(h)) "NA" else sprintf("%.2f", h), character(1))
  hr_text_anno <- ComplexHeatmap::anno_text(
    hr_label,
    just = "right",
    location = grid::unit(0.95, "npc"),
    gp = grid::gpar(fontsize = 8, col = "gray25"),
    width = grid::unit(0.95, "cm"),
    which = "row"
  )

  row_ha <- ComplexHeatmap::rowAnnotation(
    mean_cor = mean_cor_anno,
    sig_projects = sig_projects_anno,
    forest = forest_anno,
    HR = hr_text_anno,
    annotation_label = c("mean r", "sig N", "OS HR (95% CI)", "HR"),
    annotation_name_side = "top",
    annotation_name_rot = 0,
    annotation_name_gp = grid::gpar(fontsize = 8)
  )
} else {
  row_ha <- ComplexHeatmap::rowAnnotation(
    mean_cor = mean_cor_anno,
    sig_projects = sig_projects_anno,
    annotation_label = c("mean r", "sig N"),
    annotation_name_side = "top",
    annotation_name_rot = 0,
    annotation_name_gp = grid::gpar(fontsize = 8)
  )
}

ht <- ComplexHeatmap::Heatmap(
  cor_mat,
  name = "R",
  col = colors,
  row_title = "DNA Repair Genes",
  column_title = "TCGA Projects",
  cluster_rows = cluster_rows_flag,
  cluster_columns = TRUE,
  show_row_dend = show_dend,
  show_column_dend = show_dend,
  na_col = "#F0F0F0",
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 9),
  column_names_rot = 90,
  column_names_gp = grid::gpar(fontsize = 9),
  heatmap_legend_param = list(at = c(-0.5, 0, 0.5)),
  right_annotation = row_ha,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(p_mat[i, j]) & p_mat[i, j] < max_p & abs(cor_mat[i, j]) >= min_abs_cor) {
      grid::grid.text(stars[i, j], x, y, gp = grid::gpar(fontsize = 7, col = "black"))
    }
  }
)

# Output dimensions
n_genes <- length(genes)
n_projects <- ncol(cor_mat)
right_w_inches <- if (has_hr) 5.5 else 2.5  # mean_cor + sig_projects + forest + HR text
total_w <- max(12, n_projects * 0.25 + 4) + right_w_inches
total_h <- max(8, n_genes * 0.30) + 1.5  # +1.5 for column title/labels

label_suffix <- sprintf("_min%s_%s_%s", min_significant_projects, color_scheme, sort_by)
if (has_hr) {
  label_suffix <- sprintf("%s_fa", label_suffix)
}
if (!is.null(highlight_gene) && !is.na(highlight_gene) && nzchar(highlight_gene) && highlight_gene != "NULL" && highlight_gene != "null") {
  label_suffix <- sprintf("%s_top%s", label_suffix, highlight_gene)
}

if (has_hr) {
  pdf_name <- sprintf("MYBL2_gene_correlation_fa_heatmap%s.pdf", label_suffix)
  png_name <- sprintf("MYBL2_gene_correlation_fa_heatmap%s.png", label_suffix)
} else {
  pdf_name <- sprintf("MYBL2_gene_correlation_enhanced_heatmap%s.pdf", label_suffix)
  png_name <- sprintf("MYBL2_gene_correlation_enhanced_heatmap%s.png", label_suffix)
}
pdf_path <- file.path(plots_dir, pdf_name)
png_path <- file.path(plots_dir, png_name)

# --- PDF ---
grDevices::pdf(pdf_path, width = total_w, height = total_h)
ComplexHeatmap::draw(ht, merge_legend = TRUE, padding = grid::unit(c(2, 2, 2, 2), "mm"))
grDevices::dev.off()

# --- PNG ---
grDevices::png(png_path, width = total_w * 120, height = total_h * 120, res = 120)
ComplexHeatmap::draw(ht, merge_legend = TRUE, padding = grid::unit(c(2, 2, 2, 2), "mm"))
grDevices::dev.off()

if (has_hr) {
  message(sprintf("FA heatmap + forest plot saved to: %s", png_path))
} else {
  message(sprintf("Enhanced heatmap saved to: %s", png_path))
}

utils::write.csv(gene_summary, file.path(results_dir, "MYBL2_gene_correlation_summary.csv"), row.names = FALSE)

# Top positive and negative associations report
top_positive <- head(gene_summary[order(-gene_summary$mean_correlation), c("gene", "mean_correlation", "significant_projects")], 15)
top_negative <- head(gene_summary[order(gene_summary$mean_correlation), c("gene", "mean_correlation", "significant_projects")], 15)

message(sprintf("Selected genes for enhanced heatmap: %s", length(selected_genes)))
message(sprintf("Top positive: %s", paste(head(top_positive$gene, 5), collapse = ", ")))
message(sprintf("Top negative: %s", paste(head(top_negative$gene, 5), collapse = ", ")))
