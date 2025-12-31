suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(forcats)
})

# ----------------------------
# Paths
# ----------------------------
base_dir <- "outputs/wealth_goodwin/grid_search"
diag_dir <- file.path(base_dir, "diagnostics")
dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)

res_file <- file.path(base_dir, "grid_results_augmented.csv")
if (!file.exists(res_file)) res_file <- file.path(base_dir, "grid_results.csv")
stopifnot(file.exists(res_file))

hopf_file <- file.path(base_dir, "hopf_roots.csv")

res <- readr::read_csv(res_file, show_col_types = FALSE)
hopf <- if (file.exists(hopf_file)) readr::read_csv(hopf_file, show_col_types = FALSE) else tibble()

# ----------------------------
# Derived diagnostics (robust if augmented flags are missing)
# ----------------------------
if (!"stable" %in% names(res) && "maxReEig" %in% names(res)) {
  res <- res %>% mutate(stable = feasible & !is.na(maxReEig) & (maxReEig < 0))
}
if (!"oscillatory" %in% names(res) && "maxImEig" %in% names(res)) {
  res <- res %>% mutate(oscillatory = feasible & !is.na(maxImEig) & (maxImEig > 1e-6))
}
if (!"rh_ok" %in% names(res) && all(c("a1","a2","a3") %in% names(res))) {
  res <- res %>% mutate(rh_ok = feasible & !is.na(a1) & !is.na(a2) & !is.na(a3) & (a1>0) & (a2>0) & (a3>0))
}

res <- res %>%
  mutate(
    # implied local cycle period (only meaningful if oscillatory)
    period = ifelse(oscillatory, 2*pi/maxImEig, NA_real_),
    # local damping proxy (positive => stable)
    damp = ifelse(feasible, -maxReEig, NA_real_),
    regime = case_when(
      !feasible ~ "infeasible",
      stable & oscillatory ~ "stable_oscillatory",
      stable ~ "stable_real",
      !stable & oscillatory ~ "unstable_oscillatory",
      TRUE ~ "unstable_real"
    )
  )

# ----------------------------
# 1) Quick “is r_star locked?” check (your suspicion)
# ----------------------------
locked_check <- res %>%
  filter(feasible) %>%
  summarise(
    n = n(),
    r_star_min = min(r_star, na.rm = TRUE),
    r_star_max = max(r_star, na.rm = TRUE),
    r_star_range = r_star_max - r_star_min,
    omega_min = min(omega_star, na.rm = TRUE),
    omega_max = max(omega_star, na.rm = TRUE),
    d_min = min(d_star, na.rm = TRUE),
    d_max = max(d_star, na.rm = TRUE)
  )

print(locked_check)
write_csv(locked_check, file.path(diag_dir, "locked_check.csv"))

# ----------------------------
# 2) Regime counts (where the action is)
# ----------------------------
reg_counts <- res %>%
  count(regime, sort = TRUE) %>%
  mutate(share = n / sum(n))

print(reg_counts)
write_csv(reg_counts, file.path(diag_dir, "regime_counts.csv"))

# ----------------------------
# 3) Failure reasons (what kills feasibility)
# ----------------------------
reason_tab <- res %>%
  filter(!feasible) %>%
  count(reason, sort = TRUE) %>%
  mutate(share = n / sum(n))

print(reason_tab)
write_csv(reason_tab, file.path(diag_dir, "failure_reasons.csv"))

# ----------------------------
# 4) Ranges by regime (what moves vs what’s pinned)
# ----------------------------
range_by_regime <- res %>%
  group_by(regime) %>%
  summarise(
    n = n(),
    across(c(r_star, d_star, omega_star, e_star, f_star, maxReEig, maxImEig, H, period),
           list(min = ~suppressWarnings(min(.x, na.rm = TRUE)),
                med = ~suppressWarnings(median(.x, na.rm = TRUE)),
                max = ~suppressWarnings(max(.x, na.rm = TRUE))),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

write_csv(range_by_regime, file.path(diag_dir, "ranges_by_regime.csv"))

# ----------------------------
# 5) Heatmaps (auto-adapt to available grid dimensions)
# ----------------------------
pick_mid <- function(x) {
  ux <- sort(unique(x))
  ux[ceiling(length(ux)/2)]
}

tol <- 1e-12

has_phi2   <- "phi2" %in% names(res)
has_sigma  <- "sigma" %in% names(res)
has_psi    <- "psi" %in% names(res)
has_rFtil  <- "rF_tilde" %in% names(res)

phi2_slice  <- if (has_phi2)  pick_mid(res$phi2)  else NA_real_
sigma_slice <- if (has_sigma) pick_mid(res$sigma) else NA_real_
psi_slice   <- if (has_psi)   pick_mid(res$psi)   else NA_real_

# Case A: OLD grid uses (rF, rF_tilde) with phi2 slice
if (has_rFtil && has_phi2) {
  
  slice_df <- res %>% filter(abs(phi2 - phi2_slice) < tol)
  
  p_stable <- ggplot(slice_df %>% filter(feasible),
                     aes(x = rF, y = rF_tilde, fill = stable)) +
    geom_tile() +
    labs(title = paste0("Stable region (maxReEig<0), phi2 slice = ", phi2_slice),
         x = "rF", y = "rF_tilde", fill = "stable") +
    theme_minimal(base_size = 12)
  
  ggsave(file.path(diag_dir, "heat_stable_phi2slice.png"),
         p_stable, width = 7, height = 5, dpi = 200)
  
  p_period <- ggplot(slice_df %>% filter(stable, oscillatory),
                     aes(x = rF, y = rF_tilde, fill = period)) +
    geom_tile() +
    labs(title = paste0("Period 2π/Im(λ) in stable oscillatory region, phi2 = ", phi2_slice),
         x = "rF", y = "rF_tilde", fill = "period") +
    theme_minimal(base_size = 12)
  
  ggsave(file.path(diag_dir, "heat_period_stable_osc_phi2slice.png"),
         p_period, width = 7, height = 5, dpi = 200)
  
  if ("rh_ok" %in% names(slice_df)) {
    p_rh <- ggplot(slice_df, aes(x = rF, y = rF_tilde, fill = rh_ok)) +
      geom_tile() +
      labs(title = paste0("Routh–Hurwitz OK, phi2 slice = ", phi2_slice),
           x = "rF", y = "rF_tilde", fill = "rh_ok") +
      theme_minimal(base_size = 12)
    
    ggsave(file.path(diag_dir, "heat_rh_ok_phi2slice.png"),
           p_rh, width = 7, height = 5, dpi = 200)
  }
  
  # Case B: NEW grid uses (sigma, psi, phi2) plus rF
} else if (has_sigma && has_psi && has_phi2) {
  
  # slice on phi2 AND sigma, then plot (rF x psi)
  slice_df <- res %>%
    filter(abs(phi2 - phi2_slice) < tol,
           abs(sigma - sigma_slice) < tol)
  
  p_stable <- ggplot(slice_df %>% filter(feasible),
                     aes(x = rF, y = psi, fill = stable)) +
    geom_tile() +
    labs(title = paste0("Stable region (maxReEig<0), phi2=", phi2_slice,
                        ", sigma=", sigma_slice),
         x = "rF", y = "psi", fill = "stable") +
    theme_minimal(base_size = 12)
  
  ggsave(file.path(diag_dir, "heat_stable_rF_psi_phi2_sigma_slice.png"),
         p_stable, width = 7, height = 5, dpi = 200)
  
  p_period <- ggplot(slice_df %>% filter(stable, oscillatory),
                     aes(x = rF, y = psi, fill = period)) +
    geom_tile() +
    labs(title = paste0("Period 2π/Im(λ) in stable oscillatory region, phi2=",
                        phi2_slice, ", sigma=", sigma_slice),
         x = "rF", y = "psi", fill = "period") +
    theme_minimal(base_size = 12)
  
  ggsave(file.path(diag_dir, "heat_period_rF_psi_phi2_sigma_slice.png"),
         p_period, width = 7, height = 5, dpi = 200)
  
  if ("rh_ok" %in% names(slice_df)) {
    p_rh <- ggplot(slice_df, aes(x = rF, y = psi, fill = rh_ok)) +
      geom_tile() +
      labs(title = paste0("Routh–Hurwitz OK, phi2=", phi2_slice,
                          ", sigma=", sigma_slice),
           x = "rF", y = "psi", fill = "rh_ok") +
      theme_minimal(base_size = 12)
    
    ggsave(file.path(diag_dir, "heat_rh_ok_rF_psi_phi2_sigma_slice.png"),
           p_rh, width = 7, height = 5, dpi = 200)
  }
  
} else {
  message("Heatmaps skipped: could not detect expected grid columns. Found: ",
          paste(names(res), collapse = ", "))
}

# ----------------------------
# 6) Hopf plots (auto-adapt to hopf_roots schema)
# ----------------------------

# Base-R quantile binning (no scales dependency)
bin_var <- function(x, n) {
  x <- as.numeric(x)
  ux <- sort(unique(x[is.finite(x)]))
  if (is.infinite(n) || length(ux) <= n) return(factor(x, levels = ux))
  
  probs <- seq(0, 1, length.out = n + 1)
  brks  <- as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, type = 7))
  brks  <- unique(brks)
  if (length(brks) < 3) {
    brks <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE),
                length.out = min(n + 1, length(ux) + 1))
  }
  cut(x, breaks = brks, include.lowest = TRUE, right = TRUE, ordered_result = TRUE)
}

if (nrow(hopf) > 0 && "hopf_rF" %in% names(hopf)) {
  
  # OLD hopf file: (rF_tilde, phi2, hopf_rF)
  if (all(c("rF_tilde","phi2","hopf_rF") %in% names(hopf))) {
    
    p_hopf <- ggplot(hopf, aes(x = rF_tilde, y = phi2, fill = hopf_rF)) +
      geom_tile() +
      labs(title = "Hopf boundary (old schema): hopf_rF over (rF_tilde, phi2)",
           x = "rF_tilde", y = "phi2", fill = "hopf_rF") +
      theme_minimal(base_size = 12)
    
    ggsave(file.path(diag_dir, "heat_hopf_root_surface.png"),
           p_hopf, width = 7, height = 5, dpi = 200)
    
    # NEW hopf file: (sigma, psi, phi2, hopf_rF)
  } else if (all(c("sigma","psi","phi2","hopf_rF") %in% names(hopf))) {
    
    hopf2 <- hopf %>%
      filter(!is.na(hopf_rF)) %>%
      mutate(
        sigma_f = bin_var(sigma, 9),
        psi_f   = bin_var(psi,   9),
        phi2_f  = bin_var(phi2,  9)
      )
    
    # Hopf tiles: (psi x phi2) faceted by sigma
    p_hopfA <- ggplot(hopf2, aes(x = psi, y = phi2, fill = hopf_rF)) +
      geom_tile() +
      facet_wrap(~sigma_f, ncol = 3) +
      labs(title = "Hopf boundary (new schema): hopf_rF over (psi, phi2), faceted by sigma",
           x = "psi", y = "phi2", fill = "hopf_rF") +
      theme_minimal(base_size = 12)
    
    ggsave(file.path(diag_dir, "hopf_tiles_psi_phi2_by_sigma.png"),
           p_hopfA, width = 12, height = 9, dpi = 220)
    
    # Slice views (single page): fix sigma at median
    sigma_slice2 <- pick_mid(hopf2$sigma)
    p_hopfS <- ggplot(hopf2 %>% filter(abs(sigma - sigma_slice2) < tol),
                      aes(x = psi, y = phi2, fill = hopf_rF)) +
      geom_tile() +
      labs(title = paste0("Hopf slice: sigma=", sigma_slice2, " (hopf_rF over psi x phi2)"),
           x = "psi", y = "phi2", fill = "hopf_rF") +
      theme_minimal(base_size = 12)
    
    ggsave(file.path(diag_dir, "hopf_slice_sigma.png"),
           p_hopfS, width = 7, height = 5, dpi = 220)
    
  } else {
    message("Hopf plotting skipped: hopf_roots.csv schema not recognized. Found cols: ",
            paste(names(hopf), collapse = ", "))
  }
  
}

# ----------------------------
# 7) (Optional) Visual check: r_star vs omega_star vs d_star
# ----------------------------
p_r_omega <- ggplot(res %>% filter(feasible),
                    aes(x = r_star, y = omega_star)) +
  geom_point(alpha = 0.4, size = 0.8) +
  labs(title = "Feasible cloud: omega* vs r* (should reveal pinning)", x = "r_star", y = "omega_star") +
  theme_minimal(base_size = 12)

ggsave(file.path(diag_dir, "scatter_omega_vs_rstar.png"), p_r_omega, width = 7, height = 5, dpi = 200)

message("Diagnostics written to: ", diag_dir)
