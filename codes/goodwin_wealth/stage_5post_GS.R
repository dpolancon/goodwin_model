# ============================================================
# STAGE 5: SIMULATION (NUKED + SAFEGUARDS)
# - Never marks cycles if unbounded
# - Always writes Stage 5B/5C/5D outputs (schema-safe)
# - Robust to empty targets and missing columns
# - Uses simulate_model hook if available
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
})

# ----------------------------
# Safe CSV writer (schema-friendly)
# ----------------------------
write_csv_safe <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (is.null(df)) df <- tibble::tibble()
  readr::write_csv(df, path)
  invisible(path)
}

# ----------------------------
# Vector-safe logical coercion
# ----------------------------
as_logi <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  if (is.character(x)) return(toupper(x) %in% c("TRUE","T","1","YES","Y"))
  rep(FALSE, length(x))
}

# ----------------------------
# simulate_model wrapper: supports optional x0
# ----------------------------
call_simulate_model <- function(row, rF_sim, t_end, dt, x0 = NULL) {
  fm <- names(formals(simulate_model))
  if (!is.null(x0) && "x0" %in% fm) {
    simulate_model(row, rF_sim = rF_sim, t_end = t_end, dt = dt, x0 = x0)
  } else {
    simulate_model(row, rF_sim = rF_sim, t_end = t_end, dt = dt)
  }
}

# ----------------------------
# Boundedness: return BOTH ok + reason (helps diagnosis)
# ----------------------------
bounded_report <- function(sim,
                           e_rng = c(0, 1.5),
                           omega_rng = c(0, 1.5),
                           d_rng = c(-0.1, 10)) {
  sim <- tibble::as_tibble(sim)
  
  req <- c("time","e","omega","d")
  if (!all(req %in% names(sim))) {
    return(list(ok = FALSE, reason = "missing_required_cols"))
  }
  
  # nonfinite?
  if (any(!is.finite(sim$e)) || any(!is.finite(sim$omega)) || any(!is.finite(sim$d))) {
    return(list(ok = FALSE, reason = "nonfinite_values"))
  }
  
  # range checks
  if (min(sim$e) < e_rng[1] || max(sim$e) > e_rng[2]) {
    return(list(ok = FALSE, reason = "e_out_of_bounds"))
  }
  if (min(sim$omega) < omega_rng[1] || max(sim$omega) > omega_rng[2]) {
    return(list(ok = FALSE, reason = "omega_out_of_bounds"))
  }
  if (min(sim$d) < d_rng[1] || max(sim$d) > d_rng[2]) {
    return(list(ok = FALSE, reason = "d_out_of_bounds"))
  }
  
  list(ok = TRUE, reason = "ok")
}

bounded_ok <- function(sim, ...) {
  bounded_report(sim, ...)$ok
}

# ----------------------------
# Cycle metrics: NEVER cycle if unbounded
# - requires enough peaks for a period estimate
# - constant schema always
# ----------------------------
cycle_metrics <- function(sim,
                          burn_in,
                          var_min,
                          amp_min,
                          bounded = TRUE,
                          min_peaks = 4) {
  sim2 <- tibble::as_tibble(sim) %>%
    dplyr::filter(.data$time >= burn_in)
  
  out <- function(has_cycle, reason,
                  v = NA_real_, a = NA_real_, per = NA_real_, om = NA_real_, n_peaks = NA_integer_) {
    tibble::tibble(
      has_cycle = has_cycle,
      reason = reason,
      var = v,
      amp = a,
      period = per,
      omega_mean = om,
      n_peaks = n_peaks
    )
  }
  
  if (!all(c("time","omega") %in% names(sim2))) return(out(FALSE, "missing_omega_or_time"))
  
  # hard rule
  if (!isTRUE(bounded)) return(out(FALSE, "unbounded_path"))
  
  x <- sim2$omega
  t <- sim2$time
  
  v  <- stats::var(x, na.rm = TRUE)
  a  <- (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) / 2
  om <- mean(x, na.rm = TRUE)
  
  if (!is.finite(v) || !is.finite(a) || !is.finite(om)) return(out(FALSE, "nonfinite", v, a, NA_real_, om))
  
  if (v < var_min || a < amp_min) return(out(FALSE, "low_var_or_amp", v, a, NA_real_, om))
  
  dx <- diff(x)
  peaks <- which(head(dx, -1) > 0 & tail(dx, -1) <= 0) + 1
  
  if (length(peaks) < min_peaks) return(out(FALSE, "too_few_peaks", v, a, NA_real_, om, length(peaks)))
  
  periods <- diff(t[peaks])
  per <- mean(periods, na.rm = TRUE)
  
  if (!is.finite(per) || per <= 0) return(out(FALSE, "bad_period", v, a, per, om, length(peaks)))
  
  out(TRUE, "ok", v, a, per, om, length(peaks))
}

# ============================================================
# run_stage5 (NUKED + SAFEGUARDS)
# ============================================================
run_stage5 <- function(cfg, dirs) {
  cat("\n====================\nStage 5: simulation\n====================\n")
  set.seed(cfg$stage5$seed_stage5)
  
  # Ensure stage folders exist
  dir.create(dirs$stage5A, showWarnings = FALSE, recursive = TRUE)
  dir.create(dirs$stage5B, showWarnings = FALSE, recursive = TRUE)
  dir.create(dirs$stage5C, showWarnings = FALSE, recursive = TRUE)
  dir.create(dirs$stage5D, showWarnings = FALSE, recursive = TRUE)
  
  # Try to source hook if missing
  if (!exists("simulate_model")) {
    if (!is.null(cfg$stage5$hook_path) && file.exists(cfg$stage5$hook_path)) {
      source(cfg$stage5$hook_path)
    }
  }
  
  # If still missing, write empties and exit
  if (!exists("simulate_model")) {
    msg <- "Stage 5 requires simulate_model(row, rF_sim, t_end, dt[, x0]) -> tibble(time,e,omega,d)"
    if (isTRUE(cfg$stage5$require_simulate_hook)) stop(msg)
    
    cat("Skipping Stage 5 (hook missing): ", msg, "\n", sep = "")
    
    write_csv_safe(tibble::tibble(), file.path(dirs$stage5B, "stage5B_cycle_metrics.csv"))
    write_csv_safe(tibble::tibble(
      candidate_id = integer(),
      survival_rate = double(),
      bounded_rate  = double(),
      cycle_rate    = double(),
      amp_med       = double(),
      period_med    = double()
    ), file.path(dirs$stage5C, "robustness_summary.csv"))
    write_csv_safe(tibble::tibble(), file.path(dirs$stage5D, "stage5_final_ranked.csv"))
    write_csv_safe(tibble::tibble(), file.path(dirs$stage5D, "stage5_top20.csv"))
    return(invisible(NULL))
  }
  
  # Candidate pool (prefer Stage 4 relaxed shortlist; fallback progressively)
  cand_pool_file <- file.path(dirs$stage4, "stage4_shortlist_top50_GATED_RELAXED.csv")
  if (!file.exists(cand_pool_file)) cand_pool_file <- file.path(dirs$stage4, "stage4_scored_candidates_GATED_RELAXED.csv")
  if (!file.exists(cand_pool_file)) cand_pool_file <- file.path(dirs$stage4, "stage4_scored_candidates.csv")
  stopifnot(file.exists(cand_pool_file))
  
  cand <- readr::read_csv(cand_pool_file, show_col_types = FALSE) %>%
    dplyr::slice_head(n = cfg$stage5$n_candidates) %>%
    dplyr::mutate(candidate_id = dplyr::row_number())
  
  write_csv_safe(cand, file.path(dirs$stage5A, "stage5_candidate_pool.csv"))
  
  # Root seed selection: rF_root_med if present else rF
  cand_5A <- cand %>%
    dplyr::mutate(
      rF_seed = dplyr::case_when(
        "rF_root_med" %in% names(.) & is.finite(.data$rF_root_med) ~ as.numeric(.data$rF_root_med),
        TRUE ~ as.numeric(.data$rF)
      ),
      rF_root = .data$rF_seed
    )
  
  write_csv_safe(cand_5A, file.path(dirs$stage5A, "stage5A_hopf_certified.csv"))
  
  # -----------------------------------------
  # 5B: simulate below/above root and compute metrics
  # -----------------------------------------
  sim_5B <- cand_5A %>%
    dplyr::mutate(
      sim = purrr::pmap(., function(...) {
        row <- tibble::tibble(...)
        id  <- as.character(row$candidate_id)
        
        rF_root <- suppressWarnings(as.numeric(row$rF_root))
        if (!is.finite(rF_root)) rF_root <- suppressWarnings(as.numeric(row$rF_seed))
        
        rF_lo <- max(1e-6, rF_root - cfg$stage5$eps_root)
        rF_hi <- rF_root + cfg$stage5$eps_root
        
        # initial condition suggestion
        x0 <- c(e = 0.90, omega = 0.65, d = 0.5) + cfg$stage5$perturb
        if (all(c("e_star","omega_star","d_star") %in% names(row))) {
          x0 <- c(e = row$e_star, omega = row$omega_star, d = row$d_star) + cfg$stage5$perturb
        } else if (all(c("omega_star","d_star") %in% names(row))) {
          x0 <- c(e = cfg$targets$e_target, omega = row$omega_star, d = row$d_star) + cfg$stage5$perturb
        }
        
        sim_lo <- tryCatch(
          call_simulate_model(row, rF_sim = rF_lo, t_end = cfg$stage5$t_end, dt = cfg$stage5$dt, x0 = x0),
          error = function(e) NULL
        )
        sim_hi <- tryCatch(
          call_simulate_model(row, rF_sim = rF_hi, t_end = cfg$stage5$t_end, dt = cfg$stage5$dt, x0 = x0),
          error = function(e) NULL
        )
        
        if (is.null(sim_lo) || is.null(sim_hi)) {
          return(tibble::tibble(
            sim_ok = FALSE, sim_reason = "simulate_model_error",
            bounded_below = FALSE, bounded_above = FALSE,
            bounded_reason_below = "simulate_model_error",
            bounded_reason_above = "simulate_model_error",
            has_cycle_below = FALSE, has_cycle_above = FALSE,
            amp_below = NA_real_, amp_above = NA_real_,
            period_below = NA_real_, period_above = NA_real_,
            omega_mean_above = NA_real_,
            rF_lo = rF_lo, rF_hi = rF_hi
          ))
        }
        
        sim_lo <- tibble::as_tibble(sim_lo)
        sim_hi <- tibble::as_tibble(sim_hi)
        
        write_csv_safe(sim_lo, file.path(dirs$stage5B, paste0("traj_", id, "_below.csv")))
        write_csv_safe(sim_hi, file.path(dirs$stage5B, paste0("traj_", id, "_above.csv")))
        
        br_lo <- bounded_report(sim_lo)
        br_hi <- bounded_report(sim_hi)
        
        met_lo <- cycle_metrics(sim_lo,
                                burn_in = cfg$stage5$burn_in,
                                var_min = cfg$stage5$cycle_var_min,
                                amp_min = cfg$stage5$amp_min,
                                bounded = br_lo$ok,
                                min_peaks = cfg$stage5$min_peaks)
        
        met_hi <- cycle_metrics(sim_hi,
                                burn_in = cfg$stage5$burn_in,
                                var_min = cfg$stage5$cycle_var_min,
                                amp_min = cfg$stage5$amp_min,
                                bounded = br_hi$ok,
                                min_peaks = cfg$stage5$min_peaks)
        
        tibble::tibble(
          sim_ok = TRUE, sim_reason = "ok",
          bounded_below = br_lo$ok,
          bounded_above = br_hi$ok,
          bounded_reason_below = br_lo$reason,
          bounded_reason_above = br_hi$reason,
          has_cycle_below = met_lo$has_cycle,
          has_cycle_above = met_hi$has_cycle,
          amp_below = met_lo$amp,
          amp_above = met_hi$amp,
          period_below = met_lo$period,
          period_above = met_hi$period,
          omega_mean_above = met_hi$omega_mean,
          rF_lo = rF_lo, rF_hi = rF_hi
        )
      })
    ) %>%
    tidyr::unnest(sim)
  
  write_csv_safe(sim_5B, file.path(dirs$stage5B, "stage5B_cycle_metrics.csv"))
  
  # -----------------------------------------
  # 5C: robustness ring (only bounded + cycle above)
  # -----------------------------------------
  targets <- sim_5B %>%
    mutate(
      sim_ok = as_logi(sim_ok),
      bounded_above = as_logi(bounded_above),
      has_cycle_above = as_logi(has_cycle_above)
    ) %>%
    filter(sim_ok, bounded_above, has_cycle_above)
  
  robust_schema <- tibble::tibble(
    candidate_id = integer(),
    survival_rate = double(),
    bounded_rate  = double(),
    cycle_rate    = double(),
    amp_med       = double(),
    period_med    = double()
  )
  
  if (nrow(targets) == 0) {
    cat("Stage 5C: no bounded cycle candidates above root. Writing empty robustness_summary.\n")
    robust_5C <- robust_schema
  } else {
    robust_5C <- targets %>%
      mutate(
        robust = purrr::pmap(., function(...) {
          row <- tibble::tibble(...)
          
          rF_test <- suppressWarnings(as.numeric(row$rF_root))
          if (!is.finite(rF_test)) rF_test <- suppressWarnings(as.numeric(row$rF_seed))
          rF_test <- rF_test + cfg$stage5$eps_root
          
          out <- purrr::map_dfr(seq_len(cfg$stage5$robust_reps), function(rep_id) {
            rowj <- row
            
            for (nm in cfg$stage5$jitter_cols) {
              if (nm %in% names(rowj) && is.finite(rowj[[nm]])) {
                rowj[[nm]] <- rowj[[nm]] * exp(rnorm(1, 0, cfg$stage5$jitter_sd))
              }
            }
            
            sim <- tryCatch(
              call_simulate_model(rowj, rF_sim = rF_test, t_end = cfg$stage5$robust_t_end, dt = cfg$stage5$dt),
              error = function(e) NULL
            )
            
            if (is.null(sim)) {
              return(tibble::tibble(rep = rep_id, ok = FALSE, bounded = FALSE, has_cycle = FALSE,
                                    amp = NA_real_, period = NA_real_))
            }
            
            sim <- tibble::as_tibble(sim)
            br  <- bounded_report(sim)
            
            met <- cycle_metrics(sim,
                                 burn_in = min(cfg$stage5$burn_in, cfg$stage5$robust_t_end / 2),
                                 var_min = cfg$stage5$cycle_var_min,
                                 amp_min = cfg$stage5$amp_min,
                                 bounded = br$ok,
                                 min_peaks = cfg$stage5$min_peaks)
            
            tibble::tibble(
              rep = rep_id,
              ok = TRUE,
              bounded = br$ok,
              has_cycle = met$has_cycle,
              amp = met$amp,
              period = met$period
            )
          })
          
          tibble::tibble(
            survival_rate = mean(out$bounded & out$has_cycle, na.rm = TRUE),
            bounded_rate  = mean(out$bounded, na.rm = TRUE),
            cycle_rate    = mean(out$has_cycle, na.rm = TRUE),
            amp_med       = stats::median(out$amp, na.rm = TRUE),
            period_med    = stats::median(out$period, na.rm = TRUE)
          )
        })
      ) %>%
      tidyr::unnest(robust) %>%
      select(candidate_id, survival_rate, bounded_rate, cycle_rate, amp_med, period_med)
    
    # If for some reason it ends empty, enforce schema
    if (nrow(robust_5C) == 0) robust_5C <- robust_schema
  }
  
  write_csv_safe(robust_5C, file.path(dirs$stage5C, "robustness_summary.csv"))
  
  # -----------------------------------------
  # 5D: final ranking (ALWAYS produced)
  # -----------------------------------------
  final_5D <- sim_5B %>%
    left_join(robust_5C, by = "candidate_id") %>%
    mutate(
      survival_rate = dplyr::coalesce(as.numeric(survival_rate), 0),
      bounded_rate  = dplyr::coalesce(as.numeric(bounded_rate), 0),
      cycle_rate    = dplyr::coalesce(as.numeric(cycle_rate), 0),
      amp_above_num = suppressWarnings(as.numeric(amp_above)),
      score_stage5  = 3 * survival_rate + 1 * bounded_rate + 1 * cycle_rate +
        0.2 * if_else(is.finite(amp_above_num), pmin(amp_above_num, 1), 0)
    ) %>%
    arrange(desc(score_stage5))
  
  write_csv_safe(final_5D, file.path(dirs$stage5D, "stage5_final_ranked.csv"))
  write_csv_safe(final_5D %>% slice_head(n = 20), file.path(dirs$stage5D, "stage5_top20.csv"))
  
  cat("\nStage 5 complete.\n")
  cat("  5A: ", file.path(dirs$stage5A, "stage5A_hopf_certified.csv"), "\n", sep = "")
  cat("  5B: ", file.path(dirs$stage5B, "stage5B_cycle_metrics.csv"), "\n", sep = "")
  cat("  5C: ", file.path(dirs$stage5C, "robustness_summary.csv"), "\n", sep = "")
  cat("  5D: ", file.path(dirs$stage5D, "stage5_final_ranked.csv"), "\n", sep = "")
  
  invisible(final_5D)
}
