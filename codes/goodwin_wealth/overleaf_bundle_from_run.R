# ============================================================
# overleaf_bundle_from_run.R
# Build an Overleaf-ready bundle from a run folder:
# - copies figures, tex tables, and key CSVs
# - writes a starter main.tex that includes everything
# - zips the bundle for upload to Overleaf
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(purrr)
})

# ---- helpers ----
dir_create <- function(x) dir.create(x, showWarnings = FALSE, recursive = TRUE)

copy_if_exists <- function(files, to_dir) {
  files <- files[file.exists(files)]
  if (length(files) == 0) return(invisible(0L))
  file.copy(files, to_dir, overwrite = TRUE)
}

safe_list <- function(path, pattern = NULL, recursive = TRUE) {
  if (!dir.exists(path)) return(character())
  list.files(path, pattern = pattern, recursive = recursive, full.names = TRUE)
}

# ---- main bundler ----
make_overleaf_bundle <- function(
    run_base,
    out_root = "outputs/wealth_goodwin/overleaf_bundles",
    bundle_name = NULL,
    zip_bundle = TRUE
) {
  stopifnot(dir.exists(run_base))
  
  run_id <- basename(run_base)
  if (is.null(bundle_name)) bundle_name <- paste0("overleaf_", run_id)
  
  bundle_dir <- file.path(out_root, bundle_name)
  figs_dir   <- file.path(bundle_dir, "figures")
  tabs_dir   <- file.path(bundle_dir, "tables")
  data_dir   <- file.path(bundle_dir, "data")
  meta_dir   <- file.path(bundle_dir, "meta")
  
  dir_create(figs_dir)
  dir_create(tabs_dir)
  dir_create(data_dir)
  dir_create(meta_dir)
  
  # Collect files from the run folder
  all_files <- list.files(run_base, recursive = TRUE, full.names = TRUE)
  
  figs <- all_files[str_detect(tolower(all_files), "\\.(png|pdf)$")]
  tex  <- all_files[str_detect(tolower(all_files), "\\.tex$")]
  csv  <- all_files[str_detect(tolower(all_files), "\\.csv$")]
  txt  <- all_files[str_detect(tolower(all_files), "\\.txt$|\\.md$|\\.log$")]
  
  # Copy
  n_figs <- copy_if_exists(figs, figs_dir)
  n_tex  <- copy_if_exists(tex,  tabs_dir)
  n_csv  <- copy_if_exists(csv,  data_dir)
  n_txt  <- copy_if_exists(txt,  meta_dir)
  
  # Also snapshot the run manifest/cfg if present
  copy_if_exists(file.path(run_base, "manifest.csv"), meta_dir)
  copy_if_exists(file.path(run_base, "cfg_snapshot.rds"), meta_dir)
  
  # Write an Overleaf-friendly manifest of what we bundled
  bundled <- tibble(
    src_path = c(figs, tex, csv, txt),
    type = c(rep("figure", length(figs)),
             rep("table_tex", length(tex)),
             rep("data_csv", length(csv)),
             rep("meta", length(txt)))
  ) %>%
    mutate(
      file = basename(src_path),
      dest_path = case_when(
        type == "figure"    ~ file.path("figures", file),
        type == "table_tex" ~ file.path("tables", file),
        type == "data_csv"  ~ file.path("data", file),
        TRUE                ~ file.path("meta", file)
      )
    ) %>%
    arrange(type, file)
  
  write_csv(bundled, file.path(bundle_dir, "bundle_manifest.csv"))
  
  # Build a starter main.tex that auto-includes everything we copied
  fig_files <- basename(figs)
  fig_files <- fig_files[order(fig_files)]
  tex_files <- basename(tex)
  tex_files <- tex_files[order(tex_files)]
  
  # Only include figures that are image formats Overleaf handles well (PNG/PDF)
  # For PDFs: includegraphics works fine.
  fig_blocks <- if (length(fig_files) == 0) {
    c("% (No figures found in this run bundle)")
  } else {
    unlist(lapply(fig_files, function(f) {
      cap <- str_replace_all(tools::file_path_sans_ext(f), "_", " ")
      c(
        "\\begin{figure}[H]",
        "\\centering",
        paste0("\\includegraphics[width=0.95\\linewidth]{figures/", f, "}"),
        paste0("\\caption{", cap, "}"),
        paste0("\\label{fig:", tools::file_path_sans_ext(f), "}"),
        "\\end{figure}",
        ""
      )
    }))
  }
  
  tab_blocks <- if (length(tex_files) == 0) {
    c("% (No .tex tables found in this run bundle)")
  } else {
    unlist(lapply(tex_files, function(f) {
      c(
        paste0("% ---- input: ", f),
        paste0("\\input{tables/", f, "}"),
        ""
      )
    }))
  }
  
  main_tex <- c(
    "% Auto-generated Overleaf bundle",
    paste0("% Source run: ", run_base),
    "",
    "\\documentclass[11pt]{article}",
    "\\usepackage[margin=1in]{geometry}",
    "\\usepackage{graphicx}",
    "\\usepackage{booktabs}",
    "\\usepackage{float}",        # for [H]
    "\\usepackage{hyperref}",
    "",
    "\\title{Goodwin--Minsky SFC Hopf Pipeline: Run Appendix}",
    paste0("\\author{Run: ", run_id, "}"),
    "\\date{}",
    "",
    "\\begin{document}",
    "\\maketitle",
    "",
    "\\section*{What this is}",
    "This appendix compiles tables and figures exported from a single run of the Stage~5 system analysis pipeline.",
    "",
    "\\section{Tables}",
    tab_blocks,
    "",
    "\\section{Figures}",
    fig_blocks,
    "",
    "\\section*{Data}",
    "CSV outputs are included in \\texttt{data/}. If you want to show any of them as LaTeX tables, convert them in R or use Overleaf-friendly scripts.",
    "",
    "\\end{document}"
  )
  
  writeLines(main_tex, file.path(bundle_dir, "main.tex"))
  
  # Zip for Overleaf upload
  zip_path <- file.path(out_root, paste0(bundle_name, ".zip"))
  if (zip_bundle) {
    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)
    setwd(out_root)
    
    # zip the folder (relative paths)
    if (file.exists(zip_path)) file.remove(zip_path)
    utils::zip(zipfile = paste0(bundle_name, ".zip"), files = bundle_name)
  }
  
  list(
    run_base = run_base,
    bundle_dir = bundle_dir,
    zip_path = zip_path,
    counts = list(figures = length(figs), tex_tables = length(tex), csv = length(csv), meta = length(txt))
  )
}

# ---- Example usage ----
# Point this at ONE run folder you want to push to Overleaf
# e.g. run_base <- "outputs/wealth_goodwin/system_analysis/runs/run_20251230_211015"
# res <- make_overleaf_bundle(run_base)
# print(res)
# cat("Upload this zip to Overleaf:\n", res$zip_path, "\n")
