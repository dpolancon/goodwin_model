# ============================================================
# Stage 5 Outputs Manifest (folder + file-level index)
# - Scans base_dir recursively
# - Tags Stage 5 artifacts (pool, metrics, robustness, ranking, trajectories)
# - Builds wide summary for quick manipulation
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
})

make_stage5_manifest <- function(base_dir = "outputs/wealth_goodwin/grid_search",
                                 write_files = TRUE,
                                 out_name = "manifest_stage5.csv") {
  stopifnot(dir.exists(base_dir))
  
  # List all files under base_dir
  files <- list.files(base_dir, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  
  manifest <- tibble(path = files) %>%
    mutate(
      rel_path  = stringr::str_replace(path, paste0("^", stringr::fixed(base_dir), "/?"), ""),
      dir       = dirname(rel_path),
      file      = basename(rel_path),
      ext       = tools::file_ext(file),
      bytes     = suppressWarnings(file.size(path)),
      mtime     = suppressWarnings(file.info(path)$mtime),
      is_stage5 = stringr::str_detect(rel_path, "(?i)stage5") |
        stringr::str_detect(file, "(?i)traj_.*_(below|above)\\.csv$") |
        file %in% c("stage5_candidate_pool.csv",
                    "stage5A_hopf_certified.csv",
                    "stage5B_cycle_metrics.csv",
                    "robustness_summary.csv",
                    "stage5_final_ranked.csv",
                    "stage5_top20.csv")
    ) %>%
    filter(is_stage5) %>%
    mutate(
      # Classify artifacts
      artifact = case_when(
        file == "stage5_candidate_pool.csv"     ~ "5A_pool",
        file == "stage5A_hopf_certified.csv"    ~ "5A_hopf_certified",
        file == "stage5B_cycle_metrics.csv"     ~ "5B_cycle_metrics",
        file == "robustness_summary.csv"        ~ "5C_robustness_summary",
        file == "stage5_final_ranked.csv"       ~ "5D_final_ranked",
        file == "stage5_top20.csv"              ~ "5D_top20",
        str_detect(file, "(?i)^traj_\\d+_(below|above)\\.csv$") ~ "5B_trajectory",
        str_detect(file, "(?i)\\.png$")         ~ "plot_png",
        str_detect(file, "(?i)\\.pdf$")         ~ "plot_pdf",
        str_detect(file, "(?i)\\.tex$")         ~ "latex_table",
        TRUE                                    ~ "other"
      ),
      candidate_id = suppressWarnings(as.integer(str_match(file, "(?i)^traj_(\\d+)_")[,2])),
      side = case_when(
        str_detect(file, "(?i)_below\\.csv$") ~ "below",
        str_detect(file, "(?i)_above\\.csv$") ~ "above",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(artifact, rel_path)
  
  # Optional: read key Stage 5 tables if present, build a "summary handle"
  safe_read <- function(p) if (file.exists(p)) readr::read_csv(p, show_col_types = FALSE) else tibble()
  
  # Try common locations first, fall back to manifest path search
  pick_path <- function(filename) {
    hit <- manifest %>% filter(file == filename) %>% slice_head(n = 1) %>% pull(path)
    if (length(hit) == 0) return(NA_character_) else hit
  }
  
  p_metrics <- pick_path("stage5B_cycle_metrics.csv")
  p_robust  <- pick_path("robustness_summary.csv")
  p_rank    <- pick_path("stage5_final_ranked.csv")
  p_pool    <- pick_path("stage5_candidate_pool.csv")
  
  metrics <- if (!is.na(p_metrics)) safe_read(p_metrics) else tibble()
  robust  <- if (!is.na(p_robust))  safe_read(p_robust)  else tibble()
  ranked  <- if (!is.na(p_rank))    safe_read(p_rank)    else tibble()
  pool    <- if (!is.na(p_pool))    safe_read(p_pool)    else tibble()
  
  # Build summary table (never error if any are empty)
  stage5_summary <- tibble::tibble() %>%
    bind_rows(tibble(
      object = c("candidate_pool", "cycle_metrics", "robustness_summary", "final_ranked"),
      file   = c(p_pool, p_metrics, p_robust, p_rank),
      n_rows = c(nrow(pool), nrow(metrics), nrow(robust), nrow(ranked))
    )) %>%
    mutate(file = ifelse(is.na(file), "", file))
  
  # A convenient merged table if the keys exist
  merged <- tibble()
  if (nrow(metrics) > 0) {
    merged <- metrics
    if ("candidate_id" %in% names(metrics) && nrow(robust) > 0 && "candidate_id" %in% names(robust)) {
      merged <- merged %>% left_join(robust, by = "candidate_id")
    }
    if ("candidate_id" %in% names(metrics) && nrow(pool) > 0 && "candidate_id" %in% names(pool)) {
      merged <- merged %>% left_join(pool, by = "candidate_id", suffix = c("", "_pool"))
    }
    if ("candidate_id" %in% names(metrics) && nrow(ranked) > 0 && "candidate_id" %in% names(ranked)) {
      merged <- merged %>% left_join(ranked %>% select(any_of(c("candidate_id","score_stage5"))),
                                     by = "candidate_id")
    }
  }
  
  if (write_files) {
    out_manifest <- file.path(base_dir, out_name)
    readr::write_csv(manifest, out_manifest)
    
    readr::write_csv(stage5_summary, file.path(base_dir, "stage5_summary_index.csv"))
    
    # optional: export merged handle if it exists
    if (nrow(merged) > 0) {
      readr::write_csv(merged, file.path(base_dir, "stage5_merged_handle.csv"))
    }
  }
  
  list(
    manifest = manifest,
    summary  = stage5_summary,
    merged   = merged
  )
}

# ---- run ----
stage5_index <- make_stage5_manifest(
  base_dir = "outputs/wealth_goodwin/grid_search",
  write_files = TRUE
)

# quick peek
print(stage5_index$summary)
print(stage5_index$manifest %>% count(artifact, sort = TRUE))
