# ============================================================
# grid_diagnostics.R (Integrated)
# Calibration grid diagnostics for reduced 3D Goodwin–Minsky
# with external risky asset + discipline channel.
#
# KEY MODEL CHOICES
#  (A) Option 1 for finance stock dynamics:
#      dF/dt = I_F  =>  df/dt = iota_F - g f  =>  f* = iotaF*/g_n
#
#  (B) Portfolio tilt uses risk-haircut effective return:
#      rF_eff = (1-psi)*rF,  psi in [0,1]
#      Implemented as rF_tilde = -psi*rF, so rF_eff = rF + rF_tilde.
#
# GRID DIMENSIONS
#   sigma, psi, rF (Hopf scan), phi2
#
# Outputs:
#  - outputs/wealth_goodwin/grid_search/grid_results.csv
#  - outputs/wealth_goodwin/grid_search/grid_results_augmented.csv
#  - outputs/wealth_goodwin/grid_search/hopf_roots.csv
#  - outputs/wealth_goodwin/grid_search/plots/*.png
#  - outputs/wealth_goodwin/grid_search/sessionInfo.txt
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(readr)
})

# ----------------------------
# 0) Repro + paths
# ----------------------------
set.seed(123)
out_dir  <- "outputs/wealth_goodwin/grid_search"
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# GRID EXTENSIONS
# ============================================================
GRID_EXT <- list(
  im_tol = 1e-6,          # treat |Im(eig)| > im_tol as oscillatory
  hopf_H_tol = 0.0,       # 0 uses pure sign change; set small >0 for tolerance band
  write_augmented_csv = TRUE,
  augmented_csv_name = "grid_results_augmented.csv",
  hopf_csv_name = "hopf_roots.csv",
  do_extra_plots = TRUE,
  sigma_slice = NA_real_, # NA = median(sigma)
  psi_slice   = NA_real_, # NA = median(psi)
  phi2_slice  = NA_real_  # NA = median(phi2)
)

gridext_augment_results <- function(res, par_base) {
  gn0 <- par_base$alpha + par_base$beta
  
  res %>%
    mutate(
      gn = gn0,
      f_nonneg = ifelse(is.na(f_star), NA, f_star >= 0),
      
      stable = feasible & !is.na(maxReEig) & (maxReEig < 0),
      oscillatory = feasible & !is.na(maxImEig) & (maxImEig > GRID_EXT$im_tol),
      
      rh_ok = feasible &
        !is.na(a1) & !is.na(a2) & !is.na(a3) &
        (a1 > 0) & (a2 > 0) & (a3 > 0),
      
      hopf_candidate = rh_ok & !is.na(H) & (abs(H) <= max(GRID_EXT$hopf_H_tol, 0))
    )
}

# Hopf boundary in rF for each (sigma, psi, phi2) by sign change in H
gridext_hopf_roots <- function(res_aug) {
  res_aug %>%
    filter(rh_ok) %>%
    arrange(rF) %>%
    group_by(sigma, psi, phi2) %>%
    summarise(
      hopf_rF = {
        Hvec <- H
        rvec <- rF
        
        ok <- which(!is.na(Hvec) & is.finite(Hvec))
        if (length(ok) < 2) {
          NA_real_
        } else {
          H2 <- Hvec[ok]
          r2 <- rvec[ok]
          
          if (GRID_EXT$hopf_H_tol > 0) {
            near <- which(abs(H2) <= GRID_EXT$hopf_H_tol)
            if (length(near) > 0) {
              r2[near[1]]
            } else {
              idx <- which(sign(H2[-1]) != sign(H2[-length(H2)]))
              if (length(idx) == 0) NA_real_ else {
                i <- idx[1]
                r1 <- r2[i]; r2b <- r2[i+1]
                H1 <- H2[i]; H2b <- H2[i+1]
                r1 + (0 - H1) * (r2b - r1) / (H2b - H1)
              }
            }
          } else {
            idx <- which(sign(H2[-1]) != sign(H2[-length(H2)]))
            if (length(idx) == 0) NA_real_ else {
              i <- idx[1]
              r1 <- r2[i]; r2b <- r2[i+1]
              H1 <- H2[i]; H2b <- H2[i+1]
              r1 + (0 - H1) * (r2b - r1) / (H2b - H1)
            }
          }
        }
      },
      .groups = "drop"
    )
}

gridext_plot_heat <- function(df, x, y, fill_var, title = NULL) {
  ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill_var]])) +
    geom_tile() +
    labs(
      title = ifelse(is.null(title),
                     paste("Heatmap:", fill_var, "over", x, "x", y),
                     title),
      x = x, y = y, fill = fill_var
    ) +
    theme_minimal(base_size = 12)
}

gridext_postprocess_and_export <- function(res, grid, par_base, out_dir, plot_dir) {
  res_aug <- gridext_augment_results(res, par_base)
  hopf_roots <- gridext_hopf_roots(res_aug)
  
  if (isTRUE(GRID_EXT$write_augmented_csv)) {
    readr::write_csv(res_aug, file.path(out_dir, GRID_EXT$augmented_csv_name))
    readr::write_csv(hopf_roots, file.path(out_dir, GRID_EXT$hopf_csv_name))
  }
  
  if (isTRUE(GRID_EXT$do_extra_plots)) {
    sig_vals <- sort(unique(grid$sigma))
    psi_vals <- sort(unique(grid$psi))
    phi2_vals <- sort(unique(grid$phi2))
    
    sigma_slice <- GRID_EXT$sigma_slice
    if (is.na(sigma_slice)) sigma_slice <- stats::median(sig_vals)
    
    psi_slice <- GRID_EXT$psi_slice
    if (is.na(psi_slice)) psi_slice <- stats::median(psi_vals)
    
    phi2_slice <- GRID_EXT$phi2_slice
    if (is.na(phi2_slice)) phi2_slice <- stats::median(phi2_vals)
    
    slice_df <- res_aug %>%
      filter(abs(sigma - sigma_slice) < 1e-12,
             abs(psi   - psi_slice)   < 1e-12,
             abs(phi2  - phi2_slice)  < 1e-12)
    
    p_stable <- ggplot(slice_df %>% filter(feasible),
                       aes(x = rF, y = rF_eff, fill = stable)) +
      geom_tile() +
      labs(
        title = paste0("Stability slice: sigma=", sigma_slice, ", psi=", psi_slice, ", phi2=", phi2_slice),
        x = "rF", y = "rF_eff", fill = "stable"
      ) +
      theme_minimal(base_size = 12)
    
    p_rh <- ggplot(slice_df,
                   aes(x = rF, y = rF_eff, fill = rh_ok)) +
      geom_tile() +
      labs(
        title = paste0("RH-ok slice: sigma=", sigma_slice, ", psi=", psi_slice, ", phi2=", phi2_slice),
        x = "rF", y = "rF_eff", fill = "rh_ok"
      ) +
      theme_minimal(base_size = 12)
    
    ggsave(file.path(plot_dir, "ext_heat_stable_rF_vs_rFeff_slice.png"),
           p_stable, width = 7, height = 5, dpi = 200)
    ggsave(file.path(plot_dir, "ext_heat_rhOK_rF_vs_rFeff_slice.png"),
           p_rh, width = 7, height = 5, dpi = 200)
    
    # Hopf surface: for a fixed sigma & psi, map (phi2 -> hopf_rF)
    hopf_slice <- hopf_roots %>%
      filter(abs(sigma - sigma_slice) < 1e-12,
             abs(psi   - psi_slice)   < 1e-12)
    
    if (nrow(hopf_slice) > 0) {
      p_hopf <- ggplot(hopf_slice, aes(x = phi2, y = hopf_rF)) +
        geom_line() +
        geom_point(size = 1) +
        labs(
          title = paste0("Hopf boundary in rF (sigma=", sigma_slice, ", psi=", psi_slice, ")"),
          x = "phi2", y = "hopf_rF"
        ) +
        theme_minimal(base_size = 12)
      
      ggsave(file.path(plot_dir, "ext_hopf_curve_phi2_to_hopf_rF_slice.png"),
             p_hopf, width = 7, height = 5, dpi = 200)
    }
  }
  
  list(res_aug = res_aug, hopf_roots = hopf_roots)
}

# ----------------------------
# 1) Helpers: logistic + inverse
# ----------------------------
logistic_bounded <- function(x, a, b, x0, k) {
  a + (b - a) / (1 + exp(-k * (x - x0)))
}
logistic_bounded_d <- function(x, a, b, x0, k) {
  L <- 1 / (1 + exp(-k * (x - x0)))
  (b - a) * k * L * (1 - L)
}
logistic_bounded_inv <- function(y, a, b, x0, k) {
  stopifnot(all(y > a & y < b))
  x0 - (1 / k) * log((b - a) / (y - a) - 1)
}

# ----------------------------
# 2) Model primitives
# ----------------------------
par_base <- list(
  # macro
  sigma = 3.0,
  delta = 0.05,
  alpha = 0.02,
  beta  = 0.01,
  i     = 0.03,
  
  # wage-share block
  phi0 = 0.00,
  phi1 = 1.00,
  phi2 = 0.50,
  phi3 = 5.0,
  phi4 = 1.0,
  
  # investment share kappa(r): bounded logistic
  kappa_min = 0.01,
  kappa_max = 0.35,
  kappa0    = 0.05,
  kappa1    = 40.0,
  
  # portfolio tilt lambda
  lam0 = 0.00,
  lam1 = 40.0,
  
  # finance returns
  rF = 0.08,
  rF_tilde = 0.00,
  
  # NEW: haircut parameter
  psi = 0.20
)

clamp01 <- function(x) pmin(1, pmax(0, x))
rf_tilde_from_rf <- function(rF, psi) -clamp01(psi) * rF

kappa_fun <- function(r, p) logistic_bounded(r, p$kappa_min, p$kappa_max, p$kappa0, p$kappa1)
kappa_r   <- function(r, p) logistic_bounded_d(r, p$kappa_min, p$kappa_max, p$kappa0, p$kappa1)

g_fun  <- function(r, p) kappa_fun(r, p) / p$sigma - p$delta
g_rfun <- function(r, p) kappa_r(r, p) / p$sigma

lambda_fun <- function(r, p) {
  rF_eff <- p$rF + p$rF_tilde   # = (1-psi)*rF
  x <- (rF_eff - r - p$lam0)
  1 / (1 + exp(-p$lam1 * x))
}
lambda_r <- function(r, p) {
  lam <- lambda_fun(r, p)
  -p$lam1 * lam * (1 - lam)
}

Z_fun <- function(d, f, p) {
  arg <- (d - 1) + p$phi4 * (f - 1)
  1 / (1 + exp(-p$phi3 * arg))
}
Z_d_fun <- function(d, f, p) {
  Z <- Z_fun(d, f, p)
  p$phi3 * Z * (1 - Z)
}
Z_f_fun <- function(d, f, p) {
  Z <- Z_fun(d, f, p)
  p$phi3 * p$phi4 * Z * (1 - Z)
}

# ----------------------------
# 3) Steady state solver (Option 1 for f*)
# ----------------------------
steady_state <- function(p) {
  gn <- p$alpha + p$beta
  if (!is.finite(gn) || gn <= 0) return(list(ok = FALSE, reason = "gn_nonpositive"))
  
  kappa_star <- p$sigma * (gn + p$delta)
  if (!(kappa_star > p$kappa_min && kappa_star < p$kappa_max)) {
    return(list(ok = FALSE, reason = "kappa_star_out_of_bounds"))
  }
  
  r_star <- logistic_bounded_inv(kappa_star, p$kappa_min, p$kappa_max, p$kappa0, p$kappa1)
  if (!is.finite(r_star)) return(list(ok = FALSE, reason = "r_star_nonfinite"))
  
  d_star <- (kappa_star - p$sigma * r_star) / gn
  omega_star <- 1 - p$i * d_star - p$sigma * r_star
  
  lam_star <- lambda_fun(r_star, p)
  if (lam_star >= 0.999) return(list(ok = FALSE, reason = "lambda_too_close_to_1"))
  
  # steady-state funds base per K
  sF_star <- r_star / (1 - lam_star)
  
  # iotaF* = lambda* (RE/K)
  iotaF_star <- lam_star * sF_star
  
  # Option 1: f* = iotaF*/g_n
  f_star <- iotaF_star / gn
  
  Z_star <- Z_fun(d_star, f_star, p)
  e_star <- (p$alpha - p$phi0 + p$phi2 * Z_star) / p$phi1
  
  admiss <- TRUE
  reasons <- character(0)
  
  vals_chk <- c(r_star, d_star, omega_star, f_star, e_star)
  if (!all(is.finite(vals_chk))) { admiss <- FALSE; reasons <- c(reasons, "nonfinite_state") }
  if (d_star < 0)  { admiss <- FALSE; reasons <- c(reasons, "d_star_negative") }
  if (omega_star <= 0 || omega_star >= 1) { admiss <- FALSE; reasons <- c(reasons, "omega_star_outside_0_1") }
  if (e_star <= 0) { admiss <- FALSE; reasons <- c(reasons, "e_star_nonpositive") }
  if (f_star < 0)  { reasons <- c(reasons, "f_star_negative") }
  
  list(
    ok = admiss,
    reason = if (admiss) "OK" else paste(reasons, collapse = ";"),
    r_star = r_star,
    kappa_star = kappa_star,
    d_star = d_star,
    omega_star = omega_star,
    lam_star = lam_star,
    sF_star = sF_star,
    iotaF_star = iotaF_star,
    f_star = f_star,
    Z_star = Z_star,
    e_star = e_star
  )
}

# ----------------------------
# 4) Jacobian at steady state (Option 1 for f_r)
# ----------------------------
jacobian_at_ss <- function(ss, p) {
  gn <- p$alpha + p$beta
  r  <- ss$r_star
  d  <- ss$d_star
  e  <- ss$e_star
  om <- ss$omega_star
  f  <- ss$f_star
  
  r_om <- -1 / p$sigma
  r_d  <- -p$i / p$sigma
  
  kap_r <- kappa_r(r, p)
  g_r   <- kap_r / p$sigma
  
  sF <- ss$sF_star
  lam <- lambda_fun(r, p)
  lam_rv <- lambda_r(r, p)
  
  # Option 1: f(r) = (sF * lam) / g(r)
  denom <- g_fun(r, p)
  if (abs(denom) < 1e-10) denom <- 1e-10
  f_r <- sF * (lam_rv * denom - lam * g_r) / (denom^2)
  
  Zd <- Z_d_fun(d, f, p)
  Zf <- Z_f_fun(d, f, p)
  
  Z_om <- Zf * f_r * r_om
  Z_dv <- Zd + Zf * f_r * r_d
  
  J11 <- 0
  J12 <- e * g_r * r_om
  J13 <- e * g_r * r_d
  
  J21 <- om * p$phi1
  J22 <- -om * p$phi2 * Z_om
  J23 <- -om * p$phi2 * Z_dv
  
  J31 <- 0
  J32 <- 1 + kap_r * r_om - d * g_r * r_om
  J33 <- kap_r * r_d + p$i - gn - d * g_r * r_d
  
  matrix(c(J11,J12,J13,
           J21,J22,J23,
           J31,J32,J33), nrow = 3, byrow = TRUE)
}

# ----------------------------
# 5) Routh–Hurwitz / Hopf diagnostics
# ----------------------------
rh_coeffs <- function(J) {
  tr <- sum(diag(J))
  M12 <- det(J[c(1,2), c(1,2)])
  M13 <- det(J[c(1,3), c(1,3)])
  M23 <- det(J[c(2,3), c(2,3)])
  s2  <- M12 + M13 + M23
  detJ <- det(J)
  
  a1 <- -tr
  a2 <- s2
  a3 <- -detJ
  H  <- a1*a2 - a3
  
  tibble(a1 = a1, a2 = a2, a3 = a3, H = H, tr = tr, det = detJ)
}

# ----------------------------
# 6) Define the calibration grid
# ----------------------------
grid <- tidyr::expand_grid(
  sigma = seq(2.0, 5.0, by = 0.25),
  psi   = seq(0.0, 1.0, by = 0.05),
  rF    = seq(0.005, 0.20, by = 0.01),
  phi2  = seq(0.10, 1.50, by = 0.10)
)

# ----------------------------
# 7) Evaluate grid
# ----------------------------
eval_one <- function(sigma, psi, rF, phi2) {
  p <- par_base
  p$sigma <- sigma
  p$psi   <- clamp01(psi)
  p$rF    <- rF
  p$rF_tilde <- rf_tilde_from_rf(rF, p$psi)
  p$phi2  <- phi2
  
  rF_eff <- p$rF + p$rF_tilde
  
  ss <- steady_state(p)
  if (!isTRUE(ss$ok)) {
    return(tibble(
      sigma = sigma, psi = p$psi, rF = rF, rF_tilde = p$rF_tilde, rF_eff = rF_eff, phi2 = phi2,
      feasible = FALSE, reason = ss$reason,
      r_star = NA_real_, d_star = NA_real_, omega_star = NA_real_, e_star = NA_real_, f_star = NA_real_,
      maxReEig = NA_real_, maxImEig = NA_real_,
      a1 = NA_real_, a2 = NA_real_, a3 = NA_real_, H = NA_real_
    ))
  }
  
  J <- jacobian_at_ss(ss, p)
  eig <- eigen(J)$values
  
  rh <- rh_coeffs(J)
  
  tibble(
    sigma = sigma, psi = p$psi, rF = rF, rF_tilde = p$rF_tilde, rF_eff = rF_eff, phi2 = phi2,
    feasible = TRUE, reason = "OK",
    r_star = ss$r_star,
    d_star = ss$d_star,
    omega_star = ss$omega_star,
    e_star = ss$e_star,
    f_star = ss$f_star,
    maxReEig = max(Re(eig)),
    maxImEig = max(abs(Im(eig))),
    a1 = rh$a1, a2 = rh$a2, a3 = rh$a3, H = rh$H
  )
}

res <- purrr::pmap_dfr(grid, eval_one)

pp <- gridext_postprocess_and_export(res, grid, par_base, out_dir, plot_dir)
res <- pp$res_aug

# Save main results
write_csv(res, file.path(out_dir, "grid_results.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# ----------------------------
# 8) Basic plots
# ----------------------------
theme_set(theme_minimal(base_size = 12))

plot_scatter_feasible <- function(df, x, y) {
  ggplot(df, aes(x = .data[[x]], y = .data[[y]], color = feasible)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "black")) +
    labs(title = paste("Feasibility in (", x, ", ", y, ")", sep=""),
         x = x, y = y, color = "feasible")
}

plot_bivariate_density <- function(df, x, y) {
  dd <- df %>% filter(feasible)
  ggplot(dd, aes(x = .data[[x]], y = .data[[y]])) +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon",
                    alpha = 0.6, color = NA) +
    geom_point(alpha = 0.10, size = 0.5) +
    labs(title = paste("Bivariate density (feasible):", x, "vs", y),
         x = x, y = y, fill = "density") +
    guides(fill = guide_colorbar())
}

p1 <- plot_scatter_feasible(res, "rF", "phi2")
p2 <- plot_scatter_feasible(res, "sigma", "psi")
p3 <- plot_scatter_feasible(res, "rF_eff", "phi2")

ggsave(file.path(plot_dir, "feasible_rF_vs_phi2.png"), p1, width = 7, height = 5, dpi = 200)
ggsave(file.path(plot_dir, "feasible_sigma_vs_psi.png"), p2, width = 7, height = 5, dpi = 200)
ggsave(file.path(plot_dir, "feasible_rFeff_vs_phi2.png"), p3, width = 7, height = 5, dpi = 200)

d1 <- plot_bivariate_density(res, "sigma", "phi2")
d2 <- plot_bivariate_density(res, "psi", "phi2")

ggsave(file.path(plot_dir, "density_sigma_phi2.png"), d1, width = 7, height = 5, dpi = 200)
ggsave(file.path(plot_dir, "density_psi_phi2.png"), d2, width = 7, height = 5, dpi = 200)

message("Done. Wrote: grid_results.csv, grid_results_augmented.csv, hopf_roots.csv, and plots/*.png")


  
# ============================================================
# 10) "3D" Hopf visualizations in ggplot2 (2D + fill + facets)
# ============================================================

hopf_roots <- pp$hopf_roots
hopf_df <- hopf_roots %>% dplyr::filter(!is.na(hopf_rF))

if (nrow(hopf_df) == 0) {
  message("No Hopf roots found (all NA). Skipping Hopf plots.")
} else {
  
  # --- Facet control: bin continuous dimensions to avoid insane panel counts ---
  # If you want FULL facets, set these to Inf (or comment out binning).
  N_SIGMA_FACETS <- 9
  N_PSI_FACETS   <- 9
  N_PHI2_FACETS  <- 9
  
  # Binning helper (keeps labels readable)
  # Binning helper (no scales::cut_number dependency)
  bin_var <- function(x, n) {
    x <- as.numeric(x)
    ux <- sort(unique(x[is.finite(x)]))
    
    # If already small or facet limit is Inf, keep exact levels
    if (is.infinite(n) || length(ux) <= n) {
      return(factor(x, levels = ux))
    }
    
    # Quantile breaks (equal-count-ish bins)
    probs <- seq(0, 1, length.out = n + 1)
    brks  <- as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, type = 7))
    
    # Ensure strictly increasing breaks (ties happen if x has few unique values)
    brks <- unique(brks)
    if (length(brks) < 3) {
      # fallback: equal-width bins if quantiles collapse
      brks <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = min(n + 1, length(ux) + 1))
    }
    
    cut(x, breaks = brks, include.lowest = TRUE, right = TRUE, ordered_result = TRUE)
  }
  
  
  hopf_df2 <- hopf_df %>%
    mutate(
      sigma_f = bin_var(sigma, N_SIGMA_FACETS),
      psi_f   = bin_var(psi,   N_PSI_FACETS),
      phi2_f  = bin_var(phi2,  N_PHI2_FACETS)
    )
  
  theme_set(theme_minimal(base_size = 12))
  
  # ------------------------------------------------------------
  # PLOT A: (psi x phi2) with fill = hopf_rF, faceted by sigma
  # "For each sigma, what rF triggers Hopf across (psi,phi2)?"
  # ------------------------------------------------------------
  pA <- ggplot(hopf_df2, aes(x = psi, y = phi2, fill = hopf_rF)) +
    geom_tile() +
    facet_wrap(~sigma_f, ncol = 3) +
    labs(
      title = "Hopf boundary: hopf_rF over (psi, phi2), faceted by sigma",
      x = "psi (haircut)", y = "phi2 (discipline intensity)", fill = "hopf_rF"
    )
  
  ggsave(file.path(plot_dir, "hopf_tiles_psi_phi2_by_sigma.png"),
         pA, width = 12, height = 9, dpi = 220)
  
  # ------------------------------------------------------------
  # PLOT B: (sigma x phi2) with fill = hopf_rF, faceted by psi
  # "How does Hopf shift as the haircut psi changes?"
  # ------------------------------------------------------------
  pB <- ggplot(hopf_df2, aes(x = sigma, y = phi2, fill = hopf_rF)) +
    geom_tile() +
    facet_wrap(~psi_f, ncol = 3) +
    labs(
      title = "Hopf boundary: hopf_rF over (sigma, phi2), faceted by psi",
      x = "sigma (capital-output ratio)", y = "phi2", fill = "hopf_rF"
    )
  
  ggsave(file.path(plot_dir, "hopf_tiles_sigma_phi2_by_psi.png"),
         pB, width = 12, height = 9, dpi = 220)
  
  # ------------------------------------------------------------
  # PLOT C: (sigma x psi) with fill = hopf_rF, faceted by phi2
  # "For a given conflict intensity, where does Hopf exist and at what rF?"
  # ------------------------------------------------------------
  pC <- ggplot(hopf_df2, aes(x = sigma, y = psi, fill = hopf_rF)) +
    geom_tile() +
    facet_wrap(~phi2_f, ncol = 3) +
    labs(
      title = "Hopf boundary: hopf_rF over (sigma, psi), faceted by phi2",
      x = "sigma", y = "psi", fill = "hopf_rF"
    )
  
  ggsave(file.path(plot_dir, "hopf_tiles_sigma_psi_by_phi2.png"),
         pC, width = 12, height = 9, dpi = 220)
  
  # ------------------------------------------------------------
  # Optional: a cleaner single-page “3D” view using fixed slices
  # (use medians from your actual grid levels)
  # ------------------------------------------------------------
  sigma_slice <- median(sort(unique(grid$sigma)))
  psi_slice   <- median(sort(unique(grid$psi)))
  phi2_slice  <- median(sort(unique(grid$phi2)))
  
  slice_sigma <- hopf_df %>% filter(abs(sigma - sigma_slice) < 1e-12)
  slice_psi   <- hopf_df %>% filter(abs(psi   - psi_slice)   < 1e-12)
  slice_phi2  <- hopf_df %>% filter(abs(phi2  - phi2_slice)  < 1e-12)
  
  pS1 <- ggplot(slice_sigma, aes(x = psi, y = phi2, fill = hopf_rF)) +
    geom_tile() +
    labs(
      title = paste0("Slice: sigma = ", sigma_slice, " (hopf_rF over psi x phi2)"),
      x = "psi", y = "phi2", fill = "hopf_rF"
    )
  ggsave(file.path(plot_dir, "hopf_slice_sigma.png"), pS1, width = 7, height = 5, dpi = 220)
  
  pS2 <- ggplot(slice_psi, aes(x = sigma, y = phi2, fill = hopf_rF)) +
    geom_tile() +
    labs(
      title = paste0("Slice: psi = ", psi_slice, " (hopf_rF over sigma x phi2)"),
      x = "sigma", y = "phi2", fill = "hopf_rF"
    )
  ggsave(file.path(plot_dir, "hopf_slice_psi.png"), pS2, width = 7, height = 5, dpi = 220)
  
  pS3 <- ggplot(slice_phi2, aes(x = sigma, y = psi, fill = hopf_rF)) +
    geom_tile() +
    labs(
      title = paste0("Slice: phi2 = ", phi2_slice, " (hopf_rF over sigma x psi)"),
      x = "sigma", y = "psi", fill = "hopf_rF"
    )
  ggsave(file.path(plot_dir, "hopf_slice_phi2.png"), pS3, width = 7, height = 5, dpi = 220)
  
  message("Saved ggplot Hopf panels to plots/: hopf_tiles_* and hopf_slice_*")
}
