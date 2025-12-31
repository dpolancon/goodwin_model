# ============================================================
# analyze_stage5_outputs_to_analysis.R
# Reads Stage 5 outputs from: outputs/wealth_goodwin/grid_search/...
# Writes ALL analysis outputs to: outputs/wealth_goodwin/analysis
#
# Restart-safe:
# - works if robustness_summary is empty/missing
# - works if ranked file missing
# - works if Stage 5 found zero cycles
# - never errors on absent columns (fills with NA/0 as appropriate)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

analyze_stage5_to_analysis <- function(
    grid_base = "outputs/wealth_goodwin/grid_search",
    analysis_base = "outputs/wealth_goodwin/analysis",
    write_files = TRUE
) {
  stopifnot(dir.exists(grid_base))
  dir.create(analysis_base, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------------------------------------
  # Helpers
  # ----------------------------------------------------------
  pick_file <- function(filename) {
    cand <- list.files(grid_base, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
    hit <- cand[basename(cand) == filename]
    if (length(hit) == 0) return(NA_character_)
    hit[[1]]
  }
  
  safe_read <- function(path) {
    if (is.na(path) || !file.exists(path)) return(tibble())
    suppressWarnings(readr::read_csv(path, show_col_types = FALSE))
  }
  
  ensure_cols <- function(df, cols) {
    # If empty, just return (schema will be handled later)
    if (nrow(df) == 0) return(df)
    missing <- setdiff(cols, names(df))
    if (length(missing) > 0) df[missing] <- NA
    df
  }
  
  # Make sure a tibble has these columns even if it's 0 rows
  ensure_schema0 <- function(df, schema_cols) {
    for (nm in schema_cols) if (!(nm %in% names(df))) df[[nm]] <- vector(mode = "logical", length = nrow(df))
    df
  }
  
  # ----------------------------------------------------------
  # Locate Stage 5 canonical files (wherever they are)
  # ----------------------------------------------------------
  p_metrics <- pick_file("stage5B_cycle_metrics.csv")
  p_robust  <- pick_file("robustness_summary.csv")
  p_rank    <- pick_file("stage5_final_ranked.csv")
  p_pool    <- pick_file("stage5_candidate_pool.csv")
  p_hopfA   <- pick_file("stage5A_hopf_certified.csv")
  p_top20   <- pick_file("stage5_top20.csv")
  
  metrics <- safe_read(p_metrics)
  robust  <- safe_read(p_robust)
  ranked  <- safe_read(p_rank)
  pool    <- safe_read(p_pool)
  hopfA   <- safe_read(p_hopfA)
  top20   <- safe_read(p_top20)
  
  # ----------------------------------------------------------
  # Harmonize schemas
  # ----------------------------------------------------------
  metrics <- ensure_cols(metrics, c(
    "candidate_id",
    "sim_ok", "sim_reason",
    "bounded_below","bounded_above",
    "has_cycle_below","has_cycle_above",
    "amp_below","amp_above",
    "period_below","period_above",
    "omega_mean_above",
    "rF_lo","rF_hi"
  ))
  
  robust <- ensure_cols(robust, c(
    "candidate_id", "survival_rate", "bounded_rate", "cycle_rate", "amp_med", "period_med"
  ))
  
  ranked <- ensure_cols(ranked, c("candidate_id", "score_stage5"))
  
  # Coerce keys
  coerce_id <- function(df) {
    if ("candidate_id" %in% names(df)) {
      df <- df %>% mutate(candidate_id = suppressWarnings(as.integer(candidate_id)))
    }
    df
  }
  metrics <- coerce_id(metrics)
  robust  <- coerce_id(robust)
  ranked  <- coerce_id(ranked)
  pool    <- coerce_id(pool)
  hopfA   <- coerce_id(hopfA)
  top20   <- coerce_id(top20)
  
  # ----------------------------------------------------------
  # Build "decision table" (join-safe)
  # ----------------------------------------------------------
  decision <- metrics
  
  # Join robustness if possible
  if (nrow(decision) > 0 && nrow(robust) > 0 && "candidate_id" %in% names(decision) && "candidate_id" %in% names(robust)) {
    decision <- decision %>% left_join(robust, by = "candidate_id")
  } else {
    # guarantee these cols exist
    decision <- ensure_cols(decision, c("survival_rate","bounded_rate","cycle_rate","amp_med","period_med"))
  }
  
  # Join ranked if possible
  if (nrow(decision) > 0 && nrow(ranked) > 0 && "candidate_id" %in% names(decision) && "candidate_id" %in% names(ranked)) {
    decision <- decision %>% left_join(ranked %>% select(any_of(c("candidate_id","score_stage5"))), by = "candidate_id")
  } else {
    decision <- ensure_cols(decision, c("score_stage5"))
  }
  
  # Join pool parameters if possible (useful for knob maps)
  if (nrow(decision) > 0 && nrow(pool) > 0 && "candidate_id" %in% names(decision) && "candidate_id" %in% names(pool)) {
    decision <- decision %>% left_join(pool, by = "candidate_id", suffix = c("", "_pool"))
  }
  
  # Join hopf-certified view if present (optional)
  if (nrow(decision) > 0 && nrow(hopfA) > 0 && "candidate_id" %in% names(decision) && "candidate_id" %in% names(hopfA)) {
    # avoid duplicating too much
    keep_hopfA <- setdiff(names(hopfA), names(decision))
    decision <- decision %>% left_join(hopfA %>% select(any_of(c("candidate_id", keep_hopfA))), by = "candidate_id")
  }
  
  # Defaults for robustness metrics (so ranking logic doesn't break)
  decision <- decision %>%
    mutate(
      survival_rate = coalesce(as.numeric(survival_rate), 0),
      bounded_rate  = coalesce(as.numeric(bounded_rate), 0),
      cycle_rate    = coalesce(as.numeric(cycle_rate), 0),
      score_stage5  = suppressWarnings(as.numeric(score_stage5)),
      sim_ok = as.logical(sim_ok),
      bounded_above = isTRUE(bounded_above),
      bounded_below = isTRUE(bounded_below),
      has_cycle_above = isTRUE(has_cycle_above),
      has_cycle_below = isTRUE(has_cycle_below),
      cycle_flip = has_cycle_above & !has_cycle_below
    ) %>%
    arrange(desc(score_stage5), desc(cycle_rate), desc(survival_rate), desc(bounded_rate))
  
  # ----------------------------------------------------------
  # Files presence report
  # ----------------------------------------------------------
  diag <- tibble(
    object = c("stage5_candidate_pool", "stage5A_hopf_certified", "stage5B_cycle_metrics",
               "stage5C_robustness_summary", "stage5D_final_ranked", "stage5_top20"),
    file = c(p_pool, p_hopfA, p_metrics, p_robust, p_rank, p_top20),
    exists = c(file.exists(p_pool), file.exists(p_hopfA), file.exists(p_metrics),
               file.exists(p_robust), file.exists(p_rank), file.exists(p_top20)),
    n_rows = c(nrow(pool), nrow(hopfA), nrow(metrics), nrow(robust), nrow(ranked), nrow(top20))
  )
  
  counts <- decision %>%
    summarise(
      n_candidates = n(),
      n_sim_ok = sum(sim_ok, na.rm = TRUE),
      n_bounded_above = sum(bounded_above, na.rm = TRUE),
      n_cycle_above = sum(has_cycle_above, na.rm = TRUE),
      n_cycle_flip = sum(cycle_flip, na.rm = TRUE),
      n_robust_nonzero = sum((survival_rate > 0) | (bounded_rate > 0) | (cycle_rate > 0), na.rm = TRUE)
    )
  
  # ----------------------------------------------------------
  # Write outputs to analysis_base
  # ----------------------------------------------------------
  if (isTRUE(write_files)) {
    readr::write_csv(diag,     file.path(analysis_base, "stage5_files_present.csv"))
    readr::write_csv(counts,   file.path(analysis_base, "stage5_counts.csv"))
    readr::write_csv(decision, file.path(analysis_base, "stage5_decision_table.csv"))
    readr::write_csv(decision %>% slice_head(n = 20), file.path(analysis_base, "stage5_top20_decision.csv"))
    
    # Also drop a manifest of the analysis folder itself
    analysis_manifest <- tibble(path = list.files(analysis_base, recursive = TRUE, full.names = TRUE)) %>%
      mutate(
        rel_path = str_replace(path, paste0("^", fixed(analysis_base), "/?"), ""),
        bytes = suppressWarnings(file.size(path)),
        ext = tools::file_ext(basename(path)),
        mtime = suppressWarnings(file.info(path)$mtime)
      ) %>%
      arrange(rel_path)
    readr::write_csv(analysis_manifest, file.path(analysis_base, "manifest_analysis.csv"))
  }
  
  # ----------------------------------------------------------
  # Plots (only if decision has rows)
  # ----------------------------------------------------------
  if (nrow(decision) > 0) {
    p1 <- ggplot(decision, aes(x = score_stage5)) +
      geom_histogram(bins = 40, alpha = 0.80) +
      labs(title = "Stage 5: score_stage5 distribution", x = "score_stage5", y = "count") +
      theme_minimal()
    ggsave(file.path(analysis_base, "plot_score_stage5_hist.png"), p1, width = 10, height = 5, dpi = 160)
    
    p2 <- ggplot(decision, aes(x = cycle_rate, y = survival_rate)) +
      geom_point(alpha = 0.80) +
      labs(title = "Stage 5: robustness tradeoff", x = "cycle_rate", y = "survival_rate") +
      theme_minimal()
    ggsave(file.path(analysis_base, "plot_robust_tradeoff.png"), p2, width = 9, height = 6, dpi = 160)
    
    p3 <- ggplot(decision, aes(x = amp_above, y = period_above)) +
      geom_point(alpha = 0.80) +
      labs(title = "Stage 5: cycle geometry (above root)", x = "amp_above", y = "period_above") +
      theme_minimal()
    ggsave(file.path(analysis_base, "plot_cycle_geometry_above.png"), p3, width = 9, height = 6, dpi = 160)
    
    # Optional knob map if parameters are available
    knob_cols <- c("psi","phi2","kappa_max","rF")
    if (all(knob_cols %in% names(decision))) {
      p4 <- ggplot(decision, aes(x = phi2, y = psi, color = cycle_rate)) +
        geom_point(alpha = 0.85) +
        facet_wrap(~kappa_max) +
        labs(title = "Stage 5: cycle_rate over (phi2, psi) by kappa_max",
             x = "phi2", y = "psi", color = "cycle_rate") +
        theme_minimal()
      ggsave(file.path(analysis_base, "plot_cycle_rate_knobs.png"), p4, width = 12, height = 7, dpi = 160)
    }
  }
  
  # ----------------------------------------------------------
  # Short text report
  # ----------------------------------------------------------
  report_txt <- c(
    "STAGE 5 ANALYSIS REPORT",
    "======================",
    paste0("grid_base: ", grid_base),
    paste0("analysis_base: ", analysis_base),
    "",
    "FILES PRESENT:",
    paste0(diag$object, ": ", ifelse(diag$exists, "OK", "MISSING"),
           " | rows=", diag$n_rows,
           " | file=", ifelse(is.na(diag$file), "", diag$file)),
    "",
    "COUNTS:",
    paste0("n_candidates: ", counts$n_candidates),
    paste0("n_sim_ok: ", counts$n_sim_ok),
    paste0("n_bounded_above: ", counts$n_bounded_above),
    paste0("n_cycle_above: ", counts$n_cycle_above),
    paste0("n_cycle_flip (below no, above yes): ", counts$n_cycle_flip),
    paste0("n_robust_nonzero: ", counts$n_robust_nonzero),
    "",
    "WRITTEN OUTPUTS:",
    file.path(analysis_base, "stage5_decision_table.csv"),
    file.path(analysis_base, "stage5_top20_decision.csv"),
    file.path(analysis_base, "stage5_counts.csv"),
    file.path(analysis_base, "stage5_files_present.csv"),
    file.path(analysis_base, "manifest_analysis.csv")
  )
  
  if (isTRUE(write_files)) writeLines(report_txt, file.path(analysis_base, "stage5_report.txt"))
  
  # Return objects for interactive use
  list(
    diag = diag,
    counts = counts,
    decision = decision,
    analysis_base = analysis_base
  )
}

# ---- run ----
stage5_analysis <- analyze_stage5_to_analysis(
  grid_base = "outputs/wealth_goodwin/grid_search",
  analysis_base = "outputs/wealth_goodwin/analysis",
  write_files = TRUE
)

print(stage5_analysis$counts)
print(stage5_analysis$diag)
