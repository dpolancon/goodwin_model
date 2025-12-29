suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(fs)
  library(stringr)
  library(tibble)
  library(purrr)  # <- needed for keep()
})

OUT_DIR <- fs::path("outputs", "wealth_goodwin", "grid_search")

# -----------------------------
# Helper: safe read
# -----------------------------
read_csv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

# -----------------------------
# 1) Outputs manifest (review checklist)
# -----------------------------
make_outputs_manifest <- function(out_dir = OUT_DIR) {
  expected <- tribble(
    ~file, ~purpose, ~reviewed,
    "hopf_plane_scores.csv", "Raw plane ranking metrics", FALSE,
    "TABLE_plane_ranking_tidy.csv", "Tidy plane ranking table", FALSE,
    "PLOT_plane_selection_scatter.png", "Plane selection scatter plot", FALSE,
    "FIG_A_hopf_plane_signH_phi1_lambda.png", "Hopf map A (phi1, lambda)", TRUE,
    "FIG_B_hopf_plane_signH_phi1_i.png", "Hopf map B (phi1, i)", TRUE,
    "TABLE_scenarios_tail_diagnostics.csv", "Tail-window diagnostics per scenario", FALSE,
    "TABLE_scenarios_tidy.csv", "Scenario table (lambda_rel, cycle_flag)", FALSE,
    "PLOT_scenarios_amplitudes_vs_lambdaRel.png", "Amplitudes vs lambda_rel", FALSE,
    "PLOT_scenarios_period_vs_lambdaRel.png", "Period vs lambda_rel", FALSE,
    "scenario_phase_end_*", "Phase plots with start/end labels", FALSE
  )
  
  expected %>%
    mutate(
      path = fs::path(out_dir, file),
      exists = if_else(grepl("\\*$", file), TRUE, fs::file_exists(path)), # wildcard: assume generated
      status = case_when(
        reviewed ~ "Reviewed",
        exists & !reviewed ~ "Ready to review",
        !exists ~ "Missing",
        TRUE ~ "Unknown"
      )
    ) %>%
    select(file, exists, reviewed, status, purpose)
}

manifest <- make_outputs_manifest()
write_csv(manifest, fs::path(OUT_DIR, "TABLE_outputs_review_manifest.csv"))
print(manifest)

# -----------------------------
# 2) Plane ranking summary (what to actually look at)
# -----------------------------
planes <- read_csv_safe(fs::path(OUT_DIR, "TABLE_plane_ranking_tidy.csv"))
if (!is.null(planes)) {
  
  plane_summary <- planes %>%
    mutate(
      plane = paste0("(", p1, ", ", p2, ")"),
      rank_rule = case_when(
        sign_change ~ "Good: has Hopf boundary (H changes sign)",
        !sign_change ~ "Weak: no sign-change in H (likely no Hopf in this window)"
      )
    ) %>%
    arrange(desc(sign_change), desc(hopf_band_share), desc(feasible_share)) %>%
    select(plane, rank_rule, feasible_share, hopf_band_share, sign_change)
  
  write_csv(plane_summary, fs::path(OUT_DIR, "TABLE_plane_summary_for_review.csv"))
  print(plane_summary)
}

# -----------------------------
# 3) Scenario integrity check + “paper-grade” scenario summary
# -----------------------------
scen <- read_csv_safe(fs::path(OUT_DIR, "TABLE_scenarios_tidy.csv"))
if (!is.null(scen)) {
  
  scen_review <- scen %>%
    mutate(
      domain_ok = is.finite(e_mean) & is.finite(omega_mean) & is.finite(d_mean) &
        (e_mean > 0 & e_mean < 1) &
        (omega_mean > 0 & omega_mean < 1) &
        (d_mean > 0),
      regime = case_when(
        !domain_ok ~ "Blow-up / out-of-domain",
        domain_ok & cycle_flag ~ "Cycle / sustained oscillations (tail)",
        domain_ok & !cycle_flag ~ "Convergent / damped"
      )
    ) %>%
    arrange(lambda_rel) %>%
    select(run, lambda_star, lambda, lambda_rel, regime,
           e_amp, omega_amp, d_amp, u_amp, xi_amp, period_e)
  
  write_csv(scen_review, fs::path(OUT_DIR, "TABLE_scenarios_review_for_paper.csv"))
  print(scen_review)
}




OUT_DIR <- fs::path("outputs", "wealth_goodwin", "grid_search")

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

# ============================================================
# 1) PLANE AUDIT: inspect H distribution + RH magnitudes
# ============================================================
audit_plane_df <- function(df, plane_name) {
  
  dfF <- df %>% filter(ok, feasible_RH, is.finite(H))
  
  if (nrow(dfF) == 0) {
    return(tibble(
      plane = plane_name,
      n_ok = sum(df$ok, na.rm = TRUE),
      n_feasible = 0,
      H_min = NA_real_, H_q05 = NA_real_, H_med = NA_real_, H_q95 = NA_real_, H_max = NA_real_,
      share_absH_lt_1e6 = NA_real_,
      share_absH_lt_1e4 = NA_real_,
      share_absH_lt_1e3 = NA_real_,
      share_absH_lt_1e2 = NA_real_
    ))
  }
  
  tibble(
    plane = plane_name,
    n_ok = sum(df$ok, na.rm = TRUE),
    n_feasible = nrow(dfF),
    
    H_min = min(dfF$H),
    H_q05 = as.numeric(quantile(dfF$H, 0.05)),
    H_med = median(dfF$H),
    H_q95 = as.numeric(quantile(dfF$H, 0.95)),
    H_max = max(dfF$H),
    
    share_absH_lt_1e6 = mean(abs(dfF$H) < 1e-6),
    share_absH_lt_1e4 = mean(abs(dfF$H) < 1e-4),
    share_absH_lt_1e3 = mean(abs(dfF$H) < 1e-3),
    share_absH_lt_1e2 = mean(abs(dfF$H) < 1e-2)
  )
}

# Try to load plane CSVs that your driver likely created
plane_files <- fs::dir_ls(OUT_DIR, glob = "hopf_plane_*.csv")

plane_audits <- lapply(plane_files, function(f) {
  df <- read_csv_safe(f)
  nm <- fs::path_ext_remove(fs::path_file(f))
  if (is.null(df)) return(NULL)
  audit_plane_df(df, plane_name = nm)
}) %>% bind_rows()

write_csv(plane_audits, fs::path(OUT_DIR, "TABLE_plane_H_audit.csv"))
print(plane_audits)

# Quick plot: H quantiles by plane (helps catch degeneracy)
if (nrow(plane_audits) > 0) {
  pH <- plane_audits %>%
    mutate(plane = factor(plane, levels = plane)) %>%
    ggplot(aes(x = plane, y = H_med)) +
    geom_point() +
    geom_errorbar(aes(ymin = H_q05, ymax = H_q95), width = 0.15) +
    coord_flip() +
    labs(
      title = "Plane audit: H median with 5–95% interval (within feasible RH region)",
      x = NULL, y = "H"
    )
  
  ggsave(fs::path(OUT_DIR, "PLOT_plane_H_quantiles.png"), pH, width = 8.5, height = 5.2, dpi = 160)
}

# ============================================================
# 2) SCENARIO RECLASSIFICATION: use whole-path domain checks + amp thresholds
# ============================================================
# domain checks over entire simulation path
domain_flags <- function(df) {
  tibble(
    ever_bad_e = any(!is.finite(df$e) | df$e <= 0 | df$e >= 1),
    ever_bad_omega = any(!is.finite(df$omega) | df$omega <= 0 | df$omega >= 1),
    ever_bad_d = any(!is.finite(df$d) | df$d <= 0),
    ever_bad_u = any(!is.finite(df$u) | df$u <= 0)
  ) %>%
    mutate(domain_ok_all = !(ever_bad_e | ever_bad_omega | ever_bad_d | ever_bad_u))
}

tail_window <- function(df, frac = 0.25) {
  n <- nrow(df)
  i0 <- max(1, floor((1 - frac) * n))
  df[i0:n, , drop = FALSE]
}

amp <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 20) return(NA_real_)
  as.numeric(quantile(x, 0.95) - quantile(x, 0.05))
}

# classify using amplitude thresholds (tune these once, then keep fixed)
classify_scenario <- function(df, frac_tail = 0.25,
                              amp_tol_e = 0.01,
                              amp_tol_omega = 0.005) {
  
  dom <- domain_flags(df)
  
  dft <- tail_window(df, frac_tail)
  eA <- amp(dft$e)
  wA <- amp(dft$omega)
  
  if (!dom$domain_ok_all) return("Blow-up / out-of-domain")
  
  # Sustained cycle if tail amplitude is nontrivial
  if (is.finite(eA) && is.finite(wA) && (eA >= amp_tol_e || wA >= amp_tol_omega)) {
    return("Cycle / sustained (nontrivial amp)")
  }
  
  "Convergent / damped (tiny amp)"
}

# Build scenario audit from the tidy table + the sim csvs that match each run name
scen_tidy <- read_csv_safe(fs::path(OUT_DIR, "TABLE_scenarios_tidy.csv"))

if (!is.null(scen_tidy)) {
  
  scen_audit <- scen_tidy %>%
    rowwise() %>%
    mutate(
      sim_file = fs::path(OUT_DIR, paste0("scenario_sim_", run, ".csv")),
      sim_exists = fs::file_exists(sim_file),
      regime2 = if (sim_exists) {
        df <- read_csv_safe(sim_file)
        classify_scenario(df)
      } else {
        NA_character_
      }
    ) %>%
    ungroup() %>%
    select(run, lambda_star, lambda, lambda_rel, cycle_flag, regime2,
           e_amp, omega_amp, d_amp, u_amp, xi_amp, period_e, sim_exists)
  
  write_csv(scen_audit, fs::path(OUT_DIR, "TABLE_scenarios_audit_reclassified.csv"))
  print(scen_audit)
  
  # Plot: tail amplitude vs lambda_rel (your paper picture)
  pAmp <- scen_audit %>%
    ggplot(aes(x = lambda_rel, y = e_amp)) +
    geom_point() +
    labs(title = "Scenario: tail amplitude (e) vs lambda_rel", x = "lambda / lambda*", y = "e_amp (tail)")
  
  ggsave(fs::path(OUT_DIR, "PLOT_scenarios_eAmp_vs_lambdaRel_audit.png"),
         pAmp, width = 6.8, height = 4.8, dpi = 160)
}


# --- set OUT_DIR robustly (use here::here if you want, but keeping your style) ---
OUT_DIR <- fs::path("outputs", "wealth_goodwin", "grid_search")

cat("\n[debug] getwd() =", getwd(), "\n")
cat("[debug] OUT_DIR  =", OUT_DIR, "\n")
cat("[debug] OUT_DIR exists? ", fs::dir_exists(OUT_DIR), "\n\n")

# List what's actually inside OUT_DIR
if (fs::dir_exists(OUT_DIR)) {
  cat("[debug] Files in OUT_DIR:\n")
  print(fs::dir_ls(OUT_DIR))
} else {
  stop("OUT_DIR does not exist from this working directory: ", OUT_DIR)
}

read_csv_safe <- function(path) {
  if (!fs::file_exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

# ============================================================
# 1) OUTPUTS MANIFEST (UPDATED: mark what we've reviewed)
# ============================================================
make_outputs_manifest <- function(out_dir = OUT_DIR) {
  expected <- tribble(
    ~file, ~purpose, ~reviewed,
    "hopf_plane_scores.csv", "Raw plane ranking metrics", FALSE,
    "TABLE_plane_ranking_tidy.csv", "Tidy plane ranking table", FALSE,
    "PLOT_plane_selection_scatter.png", "Plane selection scatter plot", FALSE,
    "FIG_A_hopf_plane_signH_phi1_lambda.png", "Hopf map A (phi1, lambda)", TRUE,
    "FIG_B_hopf_plane_signH_phi1_i.png", "Hopf map B (phi1, i)", TRUE,
    "TABLE_scenarios_tail_diagnostics.csv", "Tail-window diagnostics per scenario", FALSE,
    "TABLE_scenarios_tidy.csv", "Scenario table (lambda_rel, cycle_flag)", TRUE,  # you reviewed it via prints
    "PLOT_scenarios_amplitudes_vs_lambdaRel.png", "Amplitudes vs lambda_rel", FALSE,
    "PLOT_scenarios_period_vs_lambdaRel.png", "Period vs lambda_rel", FALSE,
    "scenario_phase_end_*", "Phase plots with start/end labels", FALSE,
    "TABLE_scenarios_audit_reclassified.csv", "Scenario table (reclassified regime2)", TRUE,
    "TABLE_plane_H_audit.csv", "Plane-level H audit (quantiles/shares)", FALSE,
    "PLOT_plane_H_quantiles.png", "Plot of H quantiles by plane", FALSE
  )
  
  expected %>%
    mutate(
      path = fs::path(out_dir, file),
      exists = if_else(str_detect(file, "\\*$"), TRUE, fs::file_exists(path)),
      status = case_when(
        reviewed ~ "Reviewed",
        exists & !reviewed ~ "Ready to review",
        !exists ~ "Missing",
        TRUE ~ "Unknown"
      )
    ) %>%
    select(file, exists, reviewed, status, purpose)
}

manifest <- make_outputs_manifest()
write_csv(manifest, fs::path(OUT_DIR, "TABLE_outputs_review_manifest.csv"))
cat("\n--- UPDATED OUTPUTS REVIEW MANIFEST ---\n")
print(manifest)

# ============================================================
# 2) PLANE AUDIT (FIXED FILE PICKUP)
#    - Find any CSV that looks like a hopf plane grid
# ============================================================
audit_plane_df <- function(df, plane_name) {
  
  # guard: must contain these columns or it's not a plane grid
  needed <- c("ok", "feasible_RH", "H")
  if (!all(needed %in% names(df))) {
    return(tibble(
      plane = plane_name,
      note = "Skipped (not a hopf plane grid file)",
      n_ok = NA_integer_, n_feasible = NA_integer_,
      H_min = NA_real_, H_q05 = NA_real_, H_med = NA_real_, H_q95 = NA_real_, H_max = NA_real_,
      share_absH_lt_1e6 = NA_real_, share_absH_lt_1e4 = NA_real_, share_absH_lt_1e3 = NA_real_, share_absH_lt_1e2 = NA_real_
    ))
  }
  
  dfF <- df %>% filter(ok, feasible_RH, is.finite(H))
  
  if (nrow(dfF) == 0) {
    return(tibble(
      plane = plane_name,
      note = "No feasible RH points",
      n_ok = sum(df$ok, na.rm = TRUE),
      n_feasible = 0,
      H_min = NA_real_, H_q05 = NA_real_, H_med = NA_real_, H_q95 = NA_real_, H_max = NA_real_,
      share_absH_lt_1e6 = NA_real_,
      share_absH_lt_1e4 = NA_real_,
      share_absH_lt_1e3 = NA_real_,
      share_absH_lt_1e2 = NA_real_
    ))
  }
  
  tibble(
    plane = plane_name,
    note = "OK",
    n_ok = sum(df$ok, na.rm = TRUE),
    n_feasible = nrow(dfF),
    
    H_min = min(dfF$H),
    H_q05 = as.numeric(quantile(dfF$H, 0.05)),
    H_med = median(dfF$H),
    H_q95 = as.numeric(quantile(dfF$H, 0.95)),
    H_max = max(dfF$H),
    
    share_absH_lt_1e6 = mean(abs(dfF$H) < 1e-6),
    share_absH_lt_1e4 = mean(abs(dfF$H) < 1e-4),
    share_absH_lt_1e3 = mean(abs(dfF$H) < 1e-3),
    share_absH_lt_1e2 = mean(abs(dfF$H) < 1e-2)
  )
}

# robust: list all csvs and then filter by filename
csv_files <- fs::dir_ls(OUT_DIR, regexp = "\\.csv$")
bn <- basename(csv_files)

plane_files <- csv_files[
  grepl("^hopf_plane_", bn) & !grepl("scores", bn)
]
plane_files


plane_audits <- lapply(plane_files, function(f) {
  df <- read_csv_safe(f)
  nm <- fs::path_ext_remove(fs::path_file(f))
  if (is.null(df)) return(NULL)
  audit_plane_df(df, plane_name = nm)
}) %>% bind_rows()

write_csv(plane_audits, fs::path(OUT_DIR, "TABLE_plane_H_audit.csv"))
cat("\n--- PLANE H AUDIT ---\n")
print(plane_audits)

# ============================================================
# 3) DASHBOARD TABLE: one place to look
# ============================================================
plane_summary <- read_csv_safe(fs::path(OUT_DIR, "TABLE_plane_ranking_tidy.csv")) %>%
  { if (is.null(.)) tibble() else . } %>%
  mutate(plane = paste0("(", p1, ", ", p2, ")")) %>%
  select(plane, feasible_share, hopf_band_share, sign_change)

scen_audit <- read_csv_safe(fs::path(OUT_DIR, "TABLE_scenarios_audit_reclassified.csv")) %>%
  { if (is.null(.)) tibble() else . }

dashboard <- list(
  outputs_manifest = manifest,
  plane_summary = plane_summary,
  plane_H_audit = plane_audits,
  scenarios_reclassified = scen_audit
)

# save each piece in a tidy "dashboard" folder inside OUT_DIR
DASH_DIR <- fs::path(OUT_DIR, "dashboard")
fs::dir_create(DASH_DIR)

write_csv(manifest, fs::path(DASH_DIR, "01_outputs_manifest.csv"))
write_csv(plane_summary, fs::path(DASH_DIR, "02_plane_summary.csv"))
write_csv(plane_audits, fs::path(DASH_DIR, "03_plane_H_audit.csv"))
if (nrow(scen_audit) > 0) write_csv(scen_audit, fs::path(DASH_DIR, "04_scenarios_reclassified.csv"))

cat("\n[done] Dashboard exports written to: ", DASH_DIR, "\n")

