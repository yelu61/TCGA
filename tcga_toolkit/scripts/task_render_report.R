task_render_report <- function(config, ctx) {
  run_dirs <- unlist(config$run_dirs %||% character())
  if (!length(run_dirs)) {
    fail("render_report requires run_dirs array.")
  }
  format <- tolower(config$format %||% "html")
  if (format == "quarto") {
    return(render_report_quarto(config, ctx, run_dirs))
  }

  sections <- c()
  for (run_dir in run_dirs) {
    run_dir_abs <- normalizePath(run_dir, mustWork = TRUE)
    meta_path <- file.path(run_dir_abs, "run_metadata.json")
    report_md_path <- file.path(run_dir_abs, "report.md")
    results_dir <- file.path(run_dir_abs, "results")
    plots_dir <- file.path(run_dir_abs, "plots")

    meta <- if (file.exists(meta_path)) jsonlite::fromJSON(meta_path, simplifyVector = FALSE) else list(task = "unknown", task_id = basename(run_dir_abs))

    header <- sprintf("<h2>%s (%s)</h2><p><strong>Run dir:</strong> %s</p>", meta$task, meta$task_id, run_dir_abs)

    md_html <- ""
    if (file.exists(report_md_path)) {
      lines <- read_text_lines(report_md_path)
      md_html <- markdown_to_simple_html(lines)
    }

    tables_html <- ""
    if (dir.exists(results_dir)) {
      csv_files <- sort(list.files(results_dir, pattern = "\\.csv$", full.names = TRUE))
      if (length(csv_files)) {
        tables_html <- "<h3>Tables</h3>"
        for (csv in csv_files) {
          df <- tryCatch(utils::read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE, nrows = 11), error = function(e) NULL)
          if (!is.null(df) && nrow(df) > 0) {
            caption <- sprintf("<h4>%s</h4>", basename(csv))
            if (nrow(df) > 10) {
              note <- "<p><em>Showing first 10 rows.</em></p>"
              df <- head(df, 10)
            } else {
              note <- ""
            }
            tables_html <- paste0(tables_html, caption, note, dataframe_to_html(df))
          }
        }
      }
    }

    plots_html <- ""
    if (dir.exists(plots_dir)) {
      img_files <- sort(list.files(plots_dir, pattern = "\\.(png|jpg|jpeg)$", full.names = TRUE, ignore.case = TRUE))
      if (length(img_files)) {
        plots_html <- "<h3>Plots</h3><div style='display:flex;flex-wrap:wrap;gap:20px;'>"
        for (img in img_files) {
          img_uri <- image_to_data_uri(img)
          plots_html <- paste0(
            plots_html,
            sprintf("<div style='max-width:600px;'><p><strong>%s</strong></p><img src='%s' style='max-width:100%%;border:1px solid #ddd;' /></div>", basename(img), img_uri)
          )
        }
        plots_html <- paste0(plots_html, "</div>")
      }
    }

    sections <- c(sections, "<section style='margin-bottom:40px;border-bottom:2px solid #eee;padding-bottom:20px;'>", header, md_html, tables_html, plots_html, "</section>")
  }

  html <- paste0(
    "<!DOCTYPE html><html><head><meta charset='UTF-8'><title>TCGA Run Report</title>",
    "<style>",
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Helvetica,Arial,sans-serif;line-height:1.6;max-width:1200px;margin:40px auto;padding:0 20px;color:#333;}",
    "h1,h2,h3{color:#2c3e50;}table{border-collapse:collapse;width:100%%;margin-bottom:20px;}th,td{border:1px solid #ddd;padding:8px;text-align:left;font-size:14px;}th{background:#f5f5f5;}tr:nth-child(even){background:#fafafa;}",
    "code{background:#f4f4f4;padding:2px 6px;border-radius:3px;}pre{background:#f4f4f4;padding:12px;border-radius:4px;overflow-x:auto;}",
    "</style></head><body>",
    sprintf("<h1>TCGA Analysis Report</h1><p>Generated: %s</p>", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste(sections, collapse = "\n"),
    "</body></html>"
  )

  out_path <- file.path(ctx$run_dir, "report.html")
  writeLines(html, con = out_path)

  write_report(
    ctx,
    sprintf("Rendered Report for %s run(s)", length(run_dirs)),
    c(
      sprintf("- Runs summarized: `%s`", length(run_dirs)),
      "",
      "## Outputs",
      "- `report.html`"
    )
  )

  out_path
}

dataframe_to_html <- function(df) {
  if (!nrow(df)) {
    return("<p>No data.</p>")
  }
  cols <- colnames(df)
  header <- paste0("<th>", html_escape(cols), "</th>", collapse = "")
  rows <- vapply(seq_len(nrow(df)), function(i) {
    cells <- paste0("<td>", html_escape(as.character(df[i, ])), "</td>", collapse = "")
    paste0("<tr>", cells, "</tr>")
  }, character(1))
  paste0("<table><thead><tr>", header, "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table>")
}

html_escape <- function(x) {
  x <- gsub("\u0026", "\u0026amp;", x)
  x <- gsub("\u003c", "\u0026lt;", x)
  x <- gsub("\u003e", "\u0026gt;", x)
  x
}

markdown_to_simple_html <- function(lines) {
  if (!length(lines)) {
    return("")
  }
  out <- c()
  in_pre <- FALSE
  for (line in lines) {
    if (grepl("^# ", line)) {
      out <- c(out, sprintf("<h1>%s</h1>", html_escape(sub("^# ", "", line))))
    } else if (grepl("^## ", line)) {
      out <- c(out, sprintf("<h2>%s</h2>", html_escape(sub("^## ", "", line))))
    } else if (grepl("^### ", line)) {
      out <- c(out, sprintf("<h3>%s</h3>", html_escape(sub("^### ", "", line))))
    } else if (grepl("^- ", line)) {
      out <- c(out, sprintf("<li>%s</li>", html_escape(sub("^- ", "", line))))
    } else if (grepl("^```", line)) {
      in_pre <- !in_pre
      if (in_pre) {
        out <- c(out, "<pre><code>")
      } else {
        out <- c(out, "</code></pre>")
      }
    } else if (nzchar(line)) {
      if (in_pre) {
        out <- c(out, html_escape(line))
      } else {
        out <- c(out, sprintf("<p>%s</p>", html_escape(line)))
      }
    }
  }
  paste(out, collapse = "\n")
}

image_to_data_uri <- function(path) {
  ext <- tolower(tools::file_ext(path))
  mime <- switch(ext, png = "image/png", jpg = "image/jpeg", jpeg = "image/jpeg", "image/png")
  raw <- readBin(path, what = "raw", n = file.info(path)$size)
  b64 <- jsonlite::base64_enc(raw)
  sprintf("data:%s;base64,%s", mime, b64)
}

quarto_available <- function() {
  bin <- Sys.which("quarto")
  nzchar(bin) && file.exists(bin)
}

render_report_quarto <- function(config, ctx, run_dirs) {
  output_formats <- unlist(config$output_formats %||% c("html", "pdf"))
  title <- config$title %||% sprintf("TCGA Analysis Report (%s)", ctx$task_id)
  author <- config$author %||% Sys.getenv("USER", unset = "TCGA Toolkit")

  if (!quarto_available()) {
    info("Quarto binary not found on PATH; falling back to bundled HTML report.")
    legacy_config <- config
    legacy_config$format <- "html"
    return(task_render_report(legacy_config, ctx))
  }

  qmd_path <- file.path(ctx$run_dir, sprintf("%s_report.qmd", ctx$task_id))

  yaml_header <- c(
    "---",
    sprintf("title: \"%s\"", title),
    sprintf("author: \"%s\"", author),
    sprintf("date: \"%s\"", format(Sys.Date())),
    "format:",
    "  html:",
    "    embed-resources: true",
    "    toc: true",
    "    toc-depth: 3",
    "    code-fold: true",
    "    fig-cap-location: top",
    "  pdf:",
    "    documentclass: scrreprt",
    "    toc: true",
    "    toc-depth: 3",
    "    fig-cap-location: top",
    "    geometry: margin=1in",
    "execute:",
    "  echo: false",
    "  warning: false",
    "---",
    ""
  )

  body <- c()
  for (run_dir in run_dirs) {
    run_dir_abs <- normalizePath(run_dir, mustWork = TRUE)
    meta_path <- file.path(run_dir_abs, "run_metadata.json")
    report_md_path <- file.path(run_dir_abs, "report.md")
    results_dir <- file.path(run_dir_abs, "results")
    plots_dir <- file.path(run_dir_abs, "plots")

    meta <- if (file.exists(meta_path)) {
      jsonlite::fromJSON(meta_path, simplifyVector = FALSE)
    } else {
      list(task = "unknown", task_id = basename(run_dir_abs))
    }
    body <- c(body,
              sprintf("# %s (`%s`) {.unnumbered}", meta$task %||% "task", meta$task_id %||% ""),
              sprintf("*Run dir:* `%s`", run_dir_abs),
              sprintf("*Toolkit version:* `%s` *Status:* `%s` *Completed:* `%s`",
                      meta$toolkit_version %||% "?",
                      meta$status %||% "?",
                      meta$completed_at %||% "?"),
              "")
    if (file.exists(report_md_path)) {
      report_lines <- read_text_lines(report_md_path)
      report_lines <- report_lines[!grepl("^# ", report_lines)]
      body <- c(body, "## Summary", report_lines, "")
    }

    if (dir.exists(results_dir)) {
      csv_files <- sort(list.files(results_dir, pattern = "\\.csv$", full.names = TRUE))
      if (length(csv_files)) {
        body <- c(body, "## Tables")
        for (csv in csv_files) {
          body <- c(body,
                    sprintf("### `%s`", basename(csv)),
                    "```{r}",
                    sprintf("df <- utils::read.csv(\"%s\", stringsAsFactors = FALSE, check.names = FALSE, nrows = 25)", csv),
                    "knitr::kable(df)",
                    "```",
                    "")
        }
      }
    }

    if (dir.exists(plots_dir)) {
      img_files <- sort(list.files(plots_dir, pattern = "\\.(png|jpg|jpeg|pdf)$",
                                    full.names = TRUE, ignore.case = TRUE))
      png_files <- img_files[grepl("\\.png$", img_files, ignore.case = TRUE)]
      if (length(png_files)) {
        body <- c(body, "## Plots")
        for (img in png_files) {
          body <- c(body,
                    sprintf("![%s](%s){fig-alt=\"%s\"}", basename(img), img, basename(img)),
                    "")
        }
      }
    }
    body <- c(body, "\n---\n", "")
  }

  writeLines(c(yaml_header, body), qmd_path)

  rendered <- c()
  for (fmt in output_formats) {
    cmd <- sprintf("quarto render %s --to %s", shQuote(qmd_path), fmt)
    info("Running: %s", cmd)
    status <- suppressWarnings(system(cmd, intern = FALSE))
    if (status != 0) {
      warning(sprintf("quarto render to %s failed with status %s; check qmd at %s", fmt, status, qmd_path), call. = FALSE)
    } else {
      ext <- switch(fmt, html = "html", pdf = "pdf", docx = "docx", fmt)
      out <- sub("\\.qmd$", sprintf(".%s", ext), qmd_path)
      if (file.exists(out)) rendered <- c(rendered, out)
    }
  }

  write_report(
    ctx,
    sprintf("Quarto report for %s run(s)", length(run_dirs)),
    c(
      sprintf("- Runs summarized: `%s`", length(run_dirs)),
      sprintf("- Quarto: `%s`", Sys.which("quarto")),
      sprintf("- Formats rendered: `%s`", paste(output_formats, collapse = ", ")),
      sprintf("- Files: %s", paste(sprintf("`%s`", basename(rendered)), collapse = ", ")),
      "",
      "## Outputs",
      sprintf("- `%s`", basename(qmd_path)),
      paste0("- `", basename(rendered), "`", collapse = "\n")
    )
  )

  list(qmd = qmd_path, rendered = rendered)
}
