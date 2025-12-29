# ============================================================
# run_staged_grid_search_tidy.R
# Wealth-Goodwin staged grid search (Stages 1–5), refactored
# Goals:
#  - Single source of truth (cfg + dirs)
#  - Libraries loaded once
#  - Model + helpers defined once
#  - Stage outputs are restart-friendly (read/write CSV)
#  - Stage 5 hook (simulate_model) is optional: can skip gracefully or require
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(tibble)
})

# ============================================================
# 0) CONFIG (single source of truth)
# ============================================================
cfg <- list(
  base_dir = "outputs/wealth_goodwin/grid_search",
  
  run = list(stage1 = TRUE, stage2 = TRUE, stage3 = TRUE, stage4 = TRUE, stage5 = TRUE),
  
  # Targets / bands
  targets = list(
    omega_target = 0.65,
    omega_band   = c(0.62, 0.70),
    e_target     = 0.94,
    e_band       = c(0.88, 0.98)
  ),
  
  # Grids (Stage 1 backbone)
  grids = list(
    sigma_grid = c(2.25, 2.35, 2.50),
    gn_grid    = seq(0.05, 0.08, length.out = 4),
    i_grid     = seq(0.025, 0.040, length.out = 4),
    delta_grid = seq(0.03, 0.06,  length.out = 4),
    kappa_max_grid = c(0.25, 0.35, 0.45),
    
    # Stage 2 scans
    rF_grid   = seq(0.02, 0.12, length.out = 9),
    psi_grid  = c(5, 10, 15, 20),
    phi2_grid = seq(0.0, 3.0, length.out = 7),
    
    # Hopf scan around each candidate
    hopf_rF_span = 0.05,
    hopf_n_grid  = 31
  ),
  
  # Stage 3 tolerances / gates
  stage3_tol = list(
    H_bracket_eps   = 1e-10,
    eig_re_tol      = 1e-5,
    eig_im_tol      = 1e-5,
    maxRe_flip_tol  = 1e-6,
    uniroot_tol     = 1e-8,
    uniroot_maxiter = 100
  ),
  
  # Stage 4 scoring
  stage4 = list(
    omega_pref      = c(0.62, 0.65),
    omega_hard_hi   = 0.65,
    omega_hard_lo   = 0.60,
    omega_center    = 0.64,
    gate_omega_max  = 0.67,
    
    sigma_pref_max  = 2.35,
    sigma_soft_max  = 2.50,
    
    lambda_center   = 0.50,
    lambda_soft     = c(0.20, 0.80),
    
    rF_root_pref    = c(0.01, 0.12),
    rF_root_soft    = c(0.00, 0.15),
    
    gates = list(
      gate_RH_required      = TRUE,
      gate_complex_required = TRUE,
      gate_hopf_required    = FALSE
    ),
    
    weights = list(
      w_omega  = 8.0,
      w_sigma  = 1.0,
      w_lambda = 1.2,
      w_hopf   = 2.5,
      w_dist   = 2.0,
      w_RH     = 50.0,
      w_cx     = 25.0
    )
  ),
  
  # Stage 5 simulation & robustness
  stage5 = list(
    require_simulate_hook = TRUE,  # if FALSE, Stage 5 exports empty artifacts instead of stopping
    hook_path = "R/hooks/simulate_model.R",
    
    n_candidates = 40,
    eps_root     = 0.002,
    t_end        = 800,
    dt           = 0.05,
    burn_in      = 300,
    perturb      = c(e = 0.002, omega = -0.002, d = 0.02),
    
    robust_reps  = 25,
    jitter_sd    = 0.02,
    jitter_cols  = c("psi", "phi2"),
    robust_t_end = 400,
    
    cycle_var_min = 1e-5,
    amp_min       = 1e-3,
    
    seed_stage5   = 12345
  )
)

# Base parameters (single place)
par_base <- list(
  # kappa(r): logistic bounds + center/steepness
  kappa_min = 0.02,
  kappa_max = 0.25,  # overridden by scan
  kappa0    = 0.10,
  kappa1    = 30.0,
  
  # Z(d,f)
  phi3 = 8.0,
  phi4 = 1.0,
  
  # omega dynamics
  phi0  = -0.02,
  alpha = 0.02,
  
  # phi1 bounds (endogenized Stage 2/3)
  phi1_min = 0.10,
  phi1_max = 5.00
)

# ============================================================
# 1) DIRECTORIES (one coherent structure)
# ============================================================
make_dirs <- function(base_dir) {
  d <- list(
    base   = base_dir,
    stage1 = file.path(base_dir, "stage1_backbone"),
    stage2 = file.path(base_dir, "stage2_finance"),
    stage3 = file.path(base_dir, "stage3_stability"),
    stage4 = file.path(base_dir, "stage4_scoring"),
    stage5 = file.path(base_dir, "stage5")
  )
  d$stage5A <- file.path(d$stage5, "stage5A_hopf_cert")
  d$stage5B <- file.path(d$stage5, "stage5B_sim")
  d$stage5C <- file.path(d$stage5, "stage5C_robust")
  d$stage5D <- file.path(d$stage5, "stage5D_final")
  
  walk(d, ~dir.create(.x, showWarnings = FALSE, recursive = TRUE))
  d
}
dirs <- make_dirs(cfg$base_dir)

# ============================================================
# 2) GENERIC HELPERS (small, reusable)
# ============================================================
logistic <- function(x) 1 / (1 + exp(-x))

collapse_reasons <- function(reasons) {
  reasons <- unique(reasons[!is.na(reasons) & reasons != ""])
  if (length(reasons) == 0) "ok" else paste(reasons, collapse = "|")
}

to_bool <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), NA, x))
  if (is.numeric(x)) return(ifelse(is.na(x), NA, x != 0))
  if (is.character(x)) {
    y <- tolower(trimws(x))
    return(ifelse(y %in% c("true","t","1","yes","y"), TRUE,
                  ifelse(y %in% c("false","f","0","no","n"), FALSE, NA)))
  }
  rep(NA, length(x))
}

write_csv_safe <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(df, path)
  invisible(path)
}

count_reason <- function(df, col = "reason") {
  if (!col %in% names(df)) return(tibble(reason = "missing_reason_column", n = nrow(df)))
  df %>%
    separate_rows(.data[[col]], sep = "\\|") %>%
    count(.data[[col]], sort = TRUE) %>%
    rename(reason = .data[[col]])
}

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&#_{}$])", "\\\\\\1", x, perl = TRUE)
  x
}

df_to_latex_tabular <- function(df, caption = NULL, label = NULL, align = NULL) {
  stopifnot(is.data.frame(df))
  if (is.null(align)) align <- paste0(rep("l", ncol(df)), collapse = "")
  header <- paste(latex_escape(names(df)), collapse = " & ")
  body <- apply(df, 1, function(row) paste(latex_escape(row), collapse = " & "))
  body <- paste0(body, " \\\\")
  out <- c(
    "\\begin{table}[H]",
    "\\centering",
    if (!is.null(caption)) paste0("\\caption{", caption, "}"),
    if (!is.null(label))   paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{", align, "}"),
    "\\hline",
    paste0(header, " \\\\"),
    "\\hline",
    body,
    "\\hline",
    "\\end{tabular}",
    "\\end{table}"
  )
  paste(out, collapse = "\n")
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ============================================================
# 3) MODEL FUNCTIONS (one place)
# ============================================================
kappa_fun <- function(r, p) {
  p$kappa_min + (p$kappa_max - p$kappa_min) * logistic(p$kappa1 * (r - p$kappa0))
}
kappa_prime <- function(r, p) {
  s <- logistic(p$kappa1 * (r - p$kappa0))
  (p$kappa_max - p$kappa_min) * p$kappa1 * s * (1 - s)
}
kappa_inv <- function(kappa_star, p) {
  y <- (kappa_star - p$kappa_min) / (p$kappa_max - p$kappa_min)
  y <- pmin(pmax(y, 1e-12), 1 - 1e-12)
  p$kappa0 + (1 / p$kappa1) * log(y / (1 - y))
}
Z_fun <- function(d, f, p) {
  logistic(p$phi3 * ((d - 1) + p$phi4 * (f - 1)))
}
lambda_fun <- function(r, rF, psi) logistic(psi * (r - rF))

# ============================================================
# 4) CORE EVALUATOR (reused across Stage 3 + optional Stage 5)
# ============================================================
compute_at_rF <- function(row, rFj, cfg, par_base) {
  # row: list-like with sigma,g_n,i,delta,kappa_max,omega_star,d_star,psi,phi2
  row <- as.list(row)
  p <- par_base
  p$kappa_max <- row$kappa_max
  
  sigma <- row$sigma
  g_n   <- row$g_n
  i     <- row$i
  delta <- row$delta
  psi   <- row$psi
  phi2  <- row$phi2
  
  e_star <- cfg$targets$e_target
  omega  <- row$omega_star
  d      <- row$d_star
  
  # implied r from backbone identity
  r <- (1 - omega - i*d) / sigma
  if (!is.finite(r)) return(list(ok=FALSE, reason="nonfinite_r"))
  
  kappa <- kappa_fun(r, p)
  kp    <- kappa_prime(r, p)
  if (!is.finite(kappa) || !is.finite(kp)) return(list(ok=FALSE, reason="nonfinite_kappa_or_kp"))
  
  g  <- kappa / sigma - delta
  gr <- kp / sigma
  
  lam   <- lambda_fun(r, rFj, psi)
  denom <- 1 - lam
  if (!is.finite(lam) || !is.finite(denom) || abs(denom) < 1e-12) {
    return(list(ok=FALSE, reason="lambda_denominator_near_zero_or_nonfinite"))
  }
  
  q     <- lam / denom
  iotaF <- r * q
  f     <- iotaF / g_n
  if (!is.finite(f)) return(list(ok=FALSE, reason="nonfinite_f_star"))
  
  Z <- Z_fun(d, f, p)
  if (!is.finite(Z)) return(list(ok=FALSE, reason="nonfinite_Z_star"))
  
  # endogenize phi1 to hit e_target
  phi1 <- (p$alpha - p$phi0 + phi2 * Z) / e_star
  if (!is.finite(phi1)) return(list(ok=FALSE, reason="nonfinite_phi1_endog"))
  
  # derivatives for Jacobian
  lam_r   <- psi * lam * (1 - lam)
  q_r     <- lam_r / (denom^2)
  iotaF_r <- q + r*q_r
  f_r     <- iotaF_r / g_n
  if (!is.finite(f_r)) return(list(ok=FALSE, reason="nonfinite_f_r"))
  
  Z_s <- Z * (1 - Z)
  
  dr_domega <- -1 / sigma
  dr_dd     <- -i / sigma
  
  ds_domega <- p$phi3 * (p$phi4 * f_r * dr_domega)
  ds_dd     <- p$phi3 * (1 + p$phi4 * f_r * dr_dd)
  
  Z_omega <- Z_s * ds_domega
  Z_d     <- Z_s * ds_dd
  
  # Jacobian for (e, omega, d)
  J11 <- (g - g_n)
  J12 <- e_star * gr * dr_domega
  J13 <- e_star * gr * dr_dd
  
  J21 <- omega * phi1
  J22 <- omega * (-phi2 * Z_omega)
  J23 <- omega * (-phi2 * Z_d)
  
  J31 <- 0
  J32 <- 1 + kp * dr_domega - d * gr * dr_domega
  J33 <- i - g + kp * dr_dd - d * gr * dr_dd
  
  J <- matrix(c(J11,J12,J13,
                J21,J22,J23,
                J31,J32,J33), nrow=3, byrow=TRUE)
  
  if (any(!is.finite(J))) return(list(ok=FALSE, reason="nonfinite_jacobian"))
  
  ev <- tryCatch(eigen(J)$values, error=function(e) NA_complex_)
  if (all(is.na(ev))) return(list(ok=FALSE, reason="eigen_failed"))
  
  maxRe <- max(Re(ev))
  maxIm <- max(abs(Im(ev)))
  
  # RH coefficients for det(λI - J) = λ^3 + c1 λ^2 + c2 λ + c3
  c1 <- -sum(diag(J))
  c2 <- (J11*J22 - J12*J21) + (J11*J33) + (J22*J33 - J23*J32)  # J31=0 drops out cleanly
  c3 <- -det(J)
  H  <- c1*c2 - c3
  
  RH_ok <- is.finite(c1) && is.finite(c2) && is.finite(c3) && is.finite(H) &&
    (c1 > 0 && c2 > 0 && c3 > 0)
  
  list(
    ok=TRUE,
    r=r, lam=lam, f=f, Z=Z, phi1=phi1,
    J=J,
    c1=c1, c2=c2, c3=c3, H=H, RH_ok=RH_ok,
    maxRe=maxRe, maxIm=maxIm,
    stable = (is.finite(maxRe) && maxRe < 0),
    has_complex = (is.finite(maxIm) && maxIm > cfg$stage3_tol$eig_im_tol)
  )
}

# ============================================================
# 5) STAGE 1: BACKBONE FEASIBILITY (with kappa_max scan)
# ============================================================
run_stage1 <- function(cfg, dirs, par_base) {
  cat("\n====================\nStage 1: Backbone\n====================\n")
  
  g <- cfg$grids
  omega_band <- cfg$targets$omega_band
  
  grid <- tidyr::crossing(
    sigma     = g$sigma_grid,
    g_n       = g$gn_grid,
    i         = g$i_grid,
    delta     = g$delta_grid,
    kappa_max = g$kappa_max_grid
  ) %>% mutate(try_id = row_number())
  
  res <- pmap_dfr(grid, function(sigma, g_n, i, delta, kappa_max, try_id) {
    p <- par_base
    p$kappa_max <- kappa_max
    
    reasons <- character()
    
    kappa_star <- sigma * (g_n + delta)
    if (!is.finite(kappa_star)) reasons <- c(reasons, "nonfinite_kappa_star")
    
    r_star <- NA_real_
    d_star <- NA_real_
    omega_star <- NA_real_
    
    if (is.finite(kappa_star)) {
      if (kappa_star <= p$kappa_min || kappa_star >= p$kappa_max) {
        reasons <- c(reasons, "kappa_star_out_of_bounds")
      } else {
        r_star <- kappa_inv(kappa_star, p)
        if (!is.finite(r_star)) reasons <- c(reasons, "nonfinite_r_star")
        
        if (is.finite(r_star) && is.finite(g_n) && g_n > 0) {
          d_star <- (kappa_star - sigma * r_star) / g_n
        } else {
          reasons <- c(reasons, "bad_gn_for_d_star")
        }
        
        if (!is.finite(d_star)) reasons <- c(reasons, "nonfinite_d_star")
        if (is.finite(d_star) && d_star < 0) reasons <- c(reasons, "d_star_negative")
        
        if (is.finite(d_star) && is.finite(r_star)) {
          omega_star <- 1 - i * d_star - sigma * r_star
        }
        
        if (!is.finite(omega_star)) reasons <- c(reasons, "nonfinite_omega_star")
        if (is.finite(omega_star) && (omega_star <= 0 || omega_star >= 1)) reasons <- c(reasons, "omega_outside_0_1")
        if (is.finite(omega_star) && !(omega_star >= omega_band[1] && omega_star <= omega_band[2])) {
          reasons <- c(reasons, "omega_outside_target_band")
        }
      }
    }
    
    tibble(
      try_id = try_id,
      sigma = sigma, g_n = g_n, i = i, delta = delta,
      kappa_min = p$kappa_min, kappa_max = p$kappa_max,
      kappa_star = kappa_star,
      r_star = r_star,
      d_star = d_star,
      omega_star = omega_star,
      backbone_ok = (length(reasons) == 0),
      reason = collapse_reasons(reasons)
    )
  })
  
  cat("Stage 1 counts:\n")
  cat("  n_try         =", nrow(res), "\n")
  cat("  n_backbone_ok =", sum(res$backbone_ok, na.rm = TRUE), "\n")
  
  write_csv_safe(res, file.path(dirs$stage1, "stage1_backbone.csv"))
  
  rc <- count_reason(res, "reason")
  write_csv_safe(rc, file.path(dirs$stage1, "failure_reasons.csv"))
  
  rc_by_kmax <- res %>%
    separate_rows(reason, sep="\\|") %>%
    count(kappa_max, reason, sort = TRUE)
  write_csv_safe(rc_by_kmax, file.path(dirs$stage1, "failure_reasons_by_kappa_max.csv"))
  
  rate_by_kmax <- res %>%
    group_by(kappa_max) %>%
    summarise(
      n = n(),
      n_ok = sum(backbone_ok, na.rm = TRUE),
      ok_rate = n_ok / n,
      .groups = "drop"
    )
  write_csv_safe(rate_by_kmax, file.path(dirs$stage1, "backbone_rate_by_kappa_max.csv"))
  
  # Plots
  p1 <- ggplot(res, aes(x = r_star, y = omega_star, color = backbone_ok)) +
    geom_point(alpha = 0.7) +
    facet_wrap(~kappa_max) +
    labs(title = "Stage 1: omega* vs r* (by kappa_max)", x = "r*", y = "omega*") +
    theme_minimal()
  
  p2 <- ggplot(res, aes(x = r_star, y = d_star, color = backbone_ok)) +
    geom_point(alpha = 0.7) +
    facet_wrap(~kappa_max) +
    labs(title = "Stage 1: d* vs r* (by kappa_max)", x = "r*", y = "d*") +
    theme_minimal()
  
  p3 <- ggplot(res, aes(x = omega_star, fill = backbone_ok)) +
    geom_histogram(bins = 40, alpha = 0.7) +
    geom_vline(xintercept = cfg$targets$omega_band, linetype = "dashed") +
    facet_wrap(~kappa_max) +
    labs(title = "Stage 1: omega* histogram (by kappa_max)", x = "omega*", y = "count") +
    theme_minimal()
  
  p4 <- ggplot(res, aes(x = d_star, fill = backbone_ok)) +
    geom_histogram(bins = 40, alpha = 0.7) +
    facet_wrap(~kappa_max) +
    labs(title = "Stage 1: d* histogram (by kappa_max)", x = "d*", y = "count") +
    theme_minimal()
  
  p4b <- ggplot(rate_by_kmax, aes(x = factor(kappa_max), y = ok_rate)) +
    geom_col() +
    labs(title = "Stage 1: backbone OK rate by kappa_max", x = "kappa_max", y = "OK rate") +
    theme_minimal()
  
  top_reasons <- rc %>% filter(reason != "ok") %>% slice_head(n = 8) %>% pull(reason)
  p4c <- rc_by_kmax %>%
    filter(reason %in% top_reasons) %>%
    ggplot(aes(x = factor(kappa_max), y = reason, fill = n)) +
    geom_tile() +
    labs(title = "Stage 1: reason counts by kappa_max (top reasons)", x = "kappa_max", y = "reason") +
    theme_minimal()
  
  ggsave(file.path(dirs$stage1, "omega_vs_r.png"), p1, width = 10, height = 6, dpi = 160)
  ggsave(file.path(dirs$stage1, "d_vs_r.png"),     p2, width = 10, height = 6, dpi = 160)
  ggsave(file.path(dirs$stage1, "hist_omega.png"), p3, width = 10, height = 6, dpi = 160)
  ggsave(file.path(dirs$stage1, "hist_d.png"),     p4, width = 10, height = 6, dpi = 160)
  ggsave(file.path(dirs$stage1, "backbone_rate_by_kappa_max.png"), p4b, width = 8, height = 5, dpi = 160)
  ggsave(file.path(dirs$stage1, "reasons_by_kappa_max_tile.png"),  p4c, width = 10, height = 5, dpi = 160)
  
  list(stage1_res = res, stage1_ok = res %>% filter(backbone_ok))
}

# ============================================================
# 6) STAGE 2: FINANCE/DISCIPLINE (conditional on Stage 1 survivors)
# ============================================================
run_stage2 <- function(cfg, dirs, par_base, stage1_ok) {
  cat("\n====================\nStage 2: Finance/Discipline\n====================\n")
  
  if (nrow(stage1_ok) == 0) {
    cat("Stage 2: no backbone survivors. Exporting empty schema.\n")
    empty <- tibble()
    write_csv_safe(empty, file.path(dirs$stage2, "stage2_finance_scan.csv"))
    write_csv_safe(tibble(reason = "no_stage1_survivors", n = 0), file.path(dirs$stage2, "failure_reasons.csv"))
    return(list(stage2_res = empty, stage2_ok = empty))
  }
  
  g <- cfg$grids
  e_target <- cfg$targets$e_target
  
  grid <- tidyr::crossing(
    cand_id = seq_len(nrow(stage1_ok)),
    rF = g$rF_grid,
    psi = g$psi_grid,
    phi2 = g$phi2_grid
  )
  
  res <- pmap_dfr(grid, function(cand_id, rF, psi, phi2) {
    row <- stage1_ok[cand_id, ]
    p <- par_base
    p$kappa_max <- row$kappa_max
    
    reasons <- character()
    
    r_star <- row$r_star
    d_star <- row$d_star
    g_n    <- row$g_n
    
    lambda_star <- lambda_fun(r_star, rF = rF, psi = psi)
    if (!is.finite(lambda_star)) reasons <- c(reasons, "nonfinite_lambda_star")
    if (is.finite(lambda_star) && (lambda_star < 1e-4 || lambda_star > 1 - 1e-4)) reasons <- c(reasons, "lambda_saturation")
    
    f_star <- NA_real_
    if (is.finite(lambda_star) && is.finite(r_star) && is.finite(g_n) && g_n > 0) {
      if (lambda_star >= 1) {
        reasons <- c(reasons, "lambda_ge_1")
      } else {
        iotaF_star <- r_star * (lambda_star / (1 - lambda_star))
        f_star <- iotaF_star / g_n
      }
    } else {
      reasons <- c(reasons, "bad_inputs_for_f_star")
    }
    if (!is.finite(f_star)) reasons <- c(reasons, "nonfinite_f_star")
    if (is.finite(f_star) && f_star < 0) reasons <- c(reasons, "f_star_negative")
    
    Z_star <- NA_real_
    if (is.finite(d_star) && is.finite(f_star)) {
      Z_star <- Z_fun(d_star, f_star, p)
    } else {
      reasons <- c(reasons, "bad_inputs_for_Z")
    }
    if (!is.finite(Z_star)) reasons <- c(reasons, "nonfinite_Z_star")
    
    phi1_endog <- NA_real_
    if (is.finite(Z_star) && e_target > 0) {
      phi1_endog <- (p$alpha - p$phi0 + phi2 * Z_star) / e_target
    } else {
      reasons <- c(reasons, "bad_inputs_for_phi1")
    }
    if (!is.finite(phi1_endog)) reasons <- c(reasons, "nonfinite_phi1")
    if (is.finite(phi1_endog) && (phi1_endog < p$phi1_min || phi1_endog > p$phi1_max)) {
      reasons <- c(reasons, "phi1_out_of_bounds")
    }
    
    tibble(
      try_id = row$try_id,
      cand_id = cand_id,
      sigma = row$sigma, g_n = row$g_n, i = row$i, delta = row$delta,
      kappa_min = row$kappa_min, kappa_max = row$kappa_max,
      kappa_star = row$kappa_star, r_star = row$r_star, d_star = row$d_star, omega_star = row$omega_star,
      rF = rF, psi = psi, phi2 = phi2,
      lambda_star = lambda_star, f_star = f_star, Z_star = Z_star,
      phi1_endog = phi1_endog,
      econ_ok = (length(reasons) == 0),
      reason = collapse_reasons(reasons)
    )
  })
  
  cat("Stage 2 counts:\n")
  cat("  n_try     =", nrow(res), "\n")
  cat("  n_econ_ok =", sum(res$econ_ok, na.rm = TRUE), "\n")
  
  write_csv_safe(res, file.path(dirs$stage2, "stage2_finance_scan.csv"))
  write_csv_safe(count_reason(res, "reason"), file.path(dirs$stage2, "failure_reasons.csv"))
  
  top <- res %>%
    mutate(omega_gap = abs(omega_star - cfg$targets$omega_target),
           lam_mid  = abs(lambda_star - 0.5)) %>%
    arrange(desc(econ_ok), omega_gap, lam_mid) %>%
    slice_head(n = 50)
  write_csv_safe(top, file.path(dirs$stage2, "candidates_top.csv"))
  
  p5 <- ggplot(res, aes(x = rF, y = phi2, color = econ_ok)) +
    geom_point(alpha = 0.6) +
    facet_grid(psi ~ kappa_max) +
    labs(title = "Stage 2: econ_ok cloud in (rF, phi2) by psi and kappa_max", x = "rF", y = "phi2") +
    theme_minimal()
  
  p6 <- ggplot(res, aes(x = lambda_star, fill = econ_ok)) +
    geom_histogram(bins = 40, alpha = 0.7) +
    facet_wrap(~kappa_max) +
    labs(title = "Stage 2: lambda* distribution (by kappa_max)", x = "lambda*", y = "count") +
    theme_minimal()
  
  p7 <- ggplot(res, aes(x = lambda_star, y = f_star, color = econ_ok)) +
    geom_point(alpha = 0.6) +
    facet_wrap(~kappa_max) +
    labs(title = "Stage 2: f* vs lambda* (by kappa_max)", x = "lambda*", y = "f*") +
    theme_minimal()
  
  ggsave(file.path(dirs$stage2, "econ_ok_cloud.png"), p5, width = 12, height = 8, dpi = 160)
  ggsave(file.path(dirs$stage2, "lambda_dist.png"),   p6, width = 10, height = 6, dpi = 160)
  ggsave(file.path(dirs$stage2, "f_vs_lambda.png"),   p7, width = 10, height = 6, dpi = 160)
  
  list(stage2_res = res, stage2_ok = res %>% filter(econ_ok))
}

# ============================================================
# 7) STAGE 3: STABILITY + HOPF (strict boundary)
# ============================================================
run_stage3 <- function(cfg, dirs, par_base, stage2_ok) {
  cat("\n====================\nStage 3: Stability/Hopf\n====================\n")
  
  if (nrow(stage2_ok) == 0) {
    cat("Stage 3: no econ_ok candidates. Exporting empty outputs.\n")
    write_csv_safe(tibble(), file.path(dirs$stage3, "stage3_stability.csv"))
    write_csv_safe(tibble(), file.path(dirs$stage3, "hopf_roots.csv"))
    write_csv_safe(tibble(reason="no_stage2_econ_ok", n=0), file.path(dirs$stage3, "failure_reasons.csv"))
    return(list(stage3_res = tibble(), hopf_out = tibble(), hopf_verified = tibble()))
  }
  
  tol <- cfg$stage3_tol
  g <- cfg$grids
  
  # baseline stability at each candidate's rF
  stage3_res <- map_dfr(split(stage2_ok, seq_len(nrow(stage2_ok))), function(df_row) {
    row <- as.list(df_row)
    out <- compute_at_rF(row, rFj = row$rF, cfg = cfg, par_base = par_base)
    
    if (!isTRUE(out$ok)) {
      return(tibble(
        try_id=row$try_id, cand_id=row$cand_id,
        sigma=row$sigma, g_n=row$g_n, i=row$i, delta=row$delta,
        kappa_min=row$kappa_min, kappa_max=row$kappa_max,
        r_star=row$r_star, d_star=row$d_star, omega_star=row$omega_star,
        rF=row$rF, psi=row$psi, phi2=row$phi2,
        r_implied=NA_real_,
        lambda_star=NA_real_, f_star=NA_real_, Z_star=NA_real_, phi1_endog=NA_real_,
        c1=NA_real_, c2=NA_real_, c3=NA_real_, H=NA_real_,
        maxReEig=NA_real_, maxImEig=NA_real_,
        stable=NA, RH_ok=NA, has_complex=NA,
        reason=out$reason
      ))
    }
    
    tibble(
      try_id=row$try_id, cand_id=row$cand_id,
      sigma=row$sigma, g_n=row$g_n, i=row$i, delta=row$delta,
      kappa_min=row$kappa_min, kappa_max=row$kappa_max,
      r_star=row$r_star, d_star=row$d_star, omega_star=row$omega_star,
      rF=row$rF, psi=row$psi, phi2=row$phi2,
      r_implied=out$r,
      lambda_star=out$lam, f_star=out$f, Z_star=out$Z, phi1_endog=out$phi1,
      c1=out$c1, c2=out$c2, c3=out$c3, H=out$H,
      maxReEig=out$maxRe, maxImEig=out$maxIm,
      stable=out$stable, RH_ok=out$RH_ok, has_complex=out$has_complex,
      reason="ok"
    )
  })
  
  write_csv_safe(stage3_res, file.path(dirs$stage3, "stage3_stability.csv"))
  write_csv_safe(count_reason(stage3_res, "reason"), file.path(dirs$stage3, "failure_reasons.csv"))
  
  cat("Stage 3 counts:\n")
  cat("  n_try        =", nrow(stage3_res), "\n")
  cat("  n_stable     =", sum(stage3_res$stable, na.rm=TRUE), "\n")
  cat("  n_RH_ok      =", sum(stage3_res$RH_ok,  na.rm=TRUE), "\n")
  cat("  n_hasComplex =", sum(stage3_res$has_complex, na.rm=TRUE), "\n")
  
  # Hopf scan (strict bracketing + uniroot(H))
  stage3_scan <- stage3_res %>%
    filter(reason == "ok") %>%
    filter(is.finite(H), is.finite(maxReEig))
  
  hopf_rows <- list()
  
  for (k in seq_len(nrow(stage3_scan))) {
    row_df <- stage3_scan[k, ]
    row <- as.list(row_df)
    
    rF_seq <- seq(row$rF - g$hopf_rF_span, row$rF + g$hopf_rF_span, length.out = g$hopf_n_grid)
    
    H_seq     <- rep(NA_real_, length(rF_seq))
    maxRe_seq <- rep(NA_real_, length(rF_seq))
    maxIm_seq <- rep(NA_real_, length(rF_seq))
    RH_seq    <- rep(NA, length(rF_seq))
    st_seq    <- rep(NA, length(rF_seq))
    cx_seq    <- rep(NA, length(rF_seq))
    
    for (j in seq_along(rF_seq)) {
      tmp <- compute_at_rF(row, rFj = rF_seq[j], cfg = cfg, par_base = par_base)
      if (!isTRUE(tmp$ok)) next
      H_seq[j]     <- tmp$H
      maxRe_seq[j] <- tmp$maxRe
      maxIm_seq[j] <- tmp$maxIm
      RH_seq[j]    <- tmp$RH_ok
      st_seq[j]    <- tmp$stable
      cx_seq[j]    <- tmp$has_complex
    }
    
    for (m in seq_len(length(rF_seq) - 1)) {
      H1 <- H_seq[m];  H2 <- H_seq[m+1]
      if (!is.finite(H1) || !is.finite(H2)) next
      if (abs(H1) < tol$H_bracket_eps || abs(H2) < tol$H_bracket_eps) next
      if (H1 * H2 >= 0) next
      if (!(isTRUE(RH_seq[m]) && isTRUE(RH_seq[m+1]))) next
      
      R1 <- maxRe_seq[m]; R2 <- maxRe_seq[m+1]
      if (!is.finite(R1) || !is.finite(R2)) next
      
      flip_ok <- ((R1 < -tol$maxRe_flip_tol && R2 >  tol$maxRe_flip_tol) ||
                    (R2 < -tol$maxRe_flip_tol && R1 >  tol$maxRe_flip_tol) ||
                    (abs(R1) <= tol$maxRe_flip_tol || abs(R2) <= tol$maxRe_flip_tol))
      if (!flip_ok) next
      
      x1 <- rF_seq[m]; x2 <- rF_seq[m+1]
      
      H_of_rF <- function(x) {
        tmp <- compute_at_rF(row, rFj = x, cfg = cfg, par_base = par_base)
        if (!isTRUE(tmp$ok)) return(NA_real_)
        tmp$H
      }
      
      root <- tryCatch(
        uniroot(H_of_rF, lower = x1, upper = x2, tol = tol$uniroot_tol, maxiter = tol$uniroot_maxiter),
        error = function(e) NULL
      )
      if (is.null(root) || !is.finite(root$root)) next
      rF_root <- root$root
      
      tmp_root <- compute_at_rF(row, rFj = rF_root, cfg = cfg, par_base = par_base)
      if (!isTRUE(tmp_root$ok)) next
      
      hopf_ok <- TRUE
      if (!(is.finite(tmp_root$maxRe) && abs(tmp_root$maxRe) <= tol$eig_re_tol)) hopf_ok <- FALSE
      if (!isTRUE(tmp_root$has_complex)) hopf_ok <- FALSE
      
      h <- 1e-4
      tmp_a <- compute_at_rF(row, rFj = rF_root - h, cfg = cfg, par_base = par_base)
      tmp_b <- compute_at_rF(row, rFj = rF_root + h, cfg = cfg, par_base = par_base)
      dmaxRe_drF <- NA_real_
      if (isTRUE(tmp_a$ok) && isTRUE(tmp_b$ok) && is.finite(tmp_a$maxRe) && is.finite(tmp_b$maxRe)) {
        dmaxRe_drF <- (tmp_b$maxRe - tmp_a$maxRe) / (2*h)
      }
      
      hopf_rows[[length(hopf_rows) + 1]] <- tibble(
        try_id=row$try_id, cand_id=row$cand_id,
        sigma=row$sigma, g_n=row$g_n, i=row$i, delta=row$delta,
        kappa_max=row$kappa_max,
        psi=row$psi, phi2=row$phi2,
        omega_star=row$omega_star, d_star=row$d_star,
        rF0=row$rF,
        rF_left=x1, rF_right=x2,
        H_left=H1, H_right=H2,
        maxRe_left=R1, maxRe_right=R2,
        stable_left=st_seq[m], stable_right=st_seq[m+1],
        RH_left=RH_seq[m], RH_right=RH_seq[m+1],
        rF_root=rF_root,
        H_root=tmp_root$H,
        maxRe_root=tmp_root$maxRe, maxIm_root=tmp_root$maxIm,
        lambda_root=tmp_root$lam, f_root=tmp_root$f, Z_root=tmp_root$Z, phi1_root=tmp_root$phi1,
        dmaxRe_drF=dmaxRe_drF,
        hopf_ok=hopf_ok,
        reason=if (hopf_ok) "ok" else "failed_root_verification"
      )
    }
  }
  
  hopf_out <- if (length(hopf_rows) == 0) tibble() else bind_rows(hopf_rows)
  write_csv_safe(hopf_out, file.path(dirs$stage3, "hopf_roots.csv"))
  
  hopf_verified <- hopf_out %>% filter(is.finite(rF_root), hopf_ok)
  
  cat("\nStage 3 Hopf roots found (strict):\n")
  cat("  total rows    =", nrow(hopf_out), "\n")
  cat("  verified Hopf =", nrow(hopf_verified), "\n")
  
  # Plots
  p8 <- ggplot(stage3_res, aes(x = rF, y = phi2, fill = stable)) +
    geom_tile() +
    facet_grid(psi ~ kappa_max) +
    labs(title = "Stage 3: stability tiles", x = "rF", y = "phi2") +
    theme_minimal()
  ggsave(file.path(dirs$stage3, "stability_tiles.png"), p8, width = 12, height = 8, dpi = 160)
  
  p8b <- ggplot(stage3_res, aes(x = rF, y = phi2, fill = RH_ok)) +
    geom_tile() +
    facet_grid(psi ~ kappa_max) +
    labs(title = "Stage 3: RH positivity tiles", x = "rF", y = "phi2") +
    theme_minimal()
  ggsave(file.path(dirs$stage3, "RH_ok_tiles.png"), p8b, width = 12, height = 8, dpi = 160)
  
  p8c <- ggplot(stage3_res, aes(x = rF, y = phi2, fill = maxReEig)) +
    geom_tile() +
    facet_grid(psi ~ kappa_max) +
    labs(title = "Stage 3: maxRe(eigs) heatmap", x = "rF", y = "phi2") +
    theme_minimal()
  ggsave(file.path(dirs$stage3, "maxRe_tiles.png"), p8c, width = 12, height = 8, dpi = 160)
  
  p10 <- hopf_verified %>%
    ggplot(aes(x = rF_root, y = phi2)) +
    geom_point(alpha = 0.7) +
    facet_grid(psi ~ kappa_max) +
    labs(title = "Stage 3: verified Hopf boundary points", x = "rF_root", y = "phi2") +
    theme_minimal()
  ggsave(file.path(dirs$stage3, "hopf_boundary_rFroot_phi2.png"), p10, width = 12, height = 8, dpi = 160)
  
  # Tables
  tab_stage3 <- stage3_res %>%
    filter(reason == "ok") %>%
    group_by(psi, kappa_max) %>%
    summarise(
      n = n(),
      n_stable = sum(stable, na.rm=TRUE),
      share_stable = n_stable / n,
      n_RH_ok = sum(RH_ok, na.rm=TRUE),
      share_RH_ok = n_RH_ok / n,
      .groups = "drop"
    ) %>%
    arrange(psi, kappa_max)
  
  write_csv_safe(tab_stage3, file.path(dirs$stage3, "table_stage3_summary.csv"))
  tex_stage3 <- df_to_latex_tabular(tab_stage3,
                                    caption = "Stage 3 summary by $(\\psi,\\kappa_{\\max})$.",
                                    label   = "tab:stage3_summary",
                                    align   = "rrrrrr"
  )
  writeLines(tex_stage3, file.path(dirs$stage3, "table_stage3_summary.tex"))
  
  tab_hopf <- hopf_verified %>%
    group_by(psi, kappa_max) %>%
    summarise(
      n_verified = n(),
      rF_root_med = median(rF_root, na.rm=TRUE),
      rF_root_p25 = quantile(rF_root, 0.25, na.rm=TRUE),
      rF_root_p75 = quantile(rF_root, 0.75, na.rm=TRUE),
      slope_med  = median(dmaxRe_drF, na.rm=TRUE),
      .groups = "drop"
    ) %>%
    arrange(psi, kappa_max)
  
  write_csv_safe(tab_hopf, file.path(dirs$stage3, "table_hopf_summary.csv"))
  tex_hopf <- df_to_latex_tabular(tab_hopf,
                                  caption = "Verified Hopf root summary by $(\\psi,\\kappa_{\\max})$ (median/quantiles).",
                                  label   = "tab:hopf_summary",
                                  align   = "rrrrrr"
  )
  writeLines(tex_hopf, file.path(dirs$stage3, "table_hopf_summary.tex"))
  
  # Manifest
  manifest <- tibble(path = list.files(dirs$base, recursive = TRUE, full.names = TRUE)) %>%
    mutate(bytes = file.size(path), ext = tools::file_ext(path)) %>%
    arrange(path)
  write_csv_safe(manifest, file.path(dirs$base, "manifest_outputs.csv"))
  
  list(stage3_res = stage3_res, hopf_out = hopf_out, hopf_verified = hopf_verified)
}

# ============================================================
# 8) STAGE 4: SCORING + SHORTLISTS (reads Stage 3 outputs)
# ============================================================
run_stage4 <- function(cfg, dirs) {
  cat("\n====================\nStage 4: scoring\n====================\n")
  
  stage3_file <- file.path(dirs$stage3, "stage3_stability.csv")
  hopf_file   <- file.path(dirs$stage3, "hopf_roots.csv")
  stopifnot(file.exists(stage3_file), file.exists(hopf_file))
  
  stage3 <- readr::read_csv(stage3_file, show_col_types = FALSE)
  hopf   <- readr::read_csv(hopf_file,   show_col_types = FALSE)
  
  stage3_ok <- stage3 %>%
    filter(reason == "ok") %>%
    mutate(
      RH_ok       = to_bool(RH_ok),
      stable      = to_bool(stable),
      has_complex = to_bool(has_complex)
    ) %>%
    filter(is.finite(try_id), is.finite(cand_id), is.finite(rF)) %>%
    filter(is.finite(omega_star), is.finite(sigma), is.finite(lambda_star))
  
  hopf2 <- hopf %>%
    mutate(
      hopf_ok = to_bool(hopf_ok),
      rF_root = suppressWarnings(as.numeric(rF_root)),
      rF0     = suppressWarnings(as.numeric(rF0))
    )
  
  # tolerate older schema: rF0 might be missing
  if (!("rF0" %in% names(hopf2)) || all(!is.finite(hopf2$rF0))) {
    if ("rF" %in% names(hopf2)) hopf2$rF0 <- suppressWarnings(as.numeric(hopf2$rF))
  }
  
  hopf_verified <- hopf2 %>%
    filter(is.finite(rF_root), is.finite(rF0)) %>%
    filter((hopf_ok == TRUE) | (reason == "ok"))
  
  hopf_by_cand <- hopf_verified %>%
    group_by(try_id, cand_id, rF0) %>%
    summarise(
      n_roots     = n(),
      rF_root_med = median(rF_root, na.rm = TRUE),
      rF_root_p25 = quantile(rF_root, 0.25, na.rm = TRUE),
      rF_root_p75 = quantile(rF_root, 0.75, na.rm = TRUE),
      maxIm_root_med = if ("maxIm_root" %in% names(.)) median(maxIm_root, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      hopf_complex = is.finite(maxIm_root_med) & (maxIm_root_med > 1e-5)
    )
  
  df <- stage3_ok %>%
    left_join(hopf_by_cand, by = c("try_id","cand_id", "rF" = "rF0")) %>%
    mutate(
      has_hopf = is.finite(rF_root_med) & (n_roots > 0)
    )
  
  cat("Input counts:\n")
  cat("  stage3_ok rows     =", nrow(stage3_ok), "\n")
  cat("  hopf_verified rows =", nrow(hopf_verified), "\n")
  cat("  merged rows        =", nrow(df), "\n")
  
  s4 <- cfg$stage4
  w  <- s4$weights
  gates <- s4$gates
  
  df_scored <- df %>%
    mutate(
      pen_omega = case_when(
        omega_star >= s4$omega_pref[1] & omega_star <= s4$omega_pref[2] ~ 0,
        omega_star > s4$omega_hard_hi ~ 5 + 50*(omega_star - s4$omega_hard_hi),
        omega_star < s4$omega_hard_lo ~ 2 + 20*(s4$omega_hard_lo - omega_star),
        TRUE ~ abs(omega_star - s4$omega_center)
      ),
      pen_sigma = case_when(
        sigma <= s4$sigma_pref_max ~ 0,
        sigma <= s4$sigma_soft_max ~ (sigma - s4$sigma_pref_max),
        TRUE ~ 2 + 5*(sigma - s4$sigma_soft_max)
      ),
      pen_lambda = case_when(
        lambda_star >= s4$lambda_soft[1] & lambda_star <= s4$lambda_soft[2] ~ abs(lambda_star - s4$lambda_center),
        TRUE ~ 2 + 5*abs(lambda_star - s4$lambda_center)
      ),
      pen_hopf = case_when(
        has_hopf & rF_root_med >= s4$rF_root_pref[1] & rF_root_med <= s4$rF_root_pref[2] ~ 0,
        has_hopf & rF_root_med >= s4$rF_root_soft[1] & rF_root_med <= s4$rF_root_soft[2] ~ abs(rF_root_med - median(s4$rF_root_pref)),
        has_hopf ~ 2 + 5*abs(rF_root_med - median(s4$rF_root_pref)),
        TRUE ~ 10
      ),
      pen_dist = case_when(
        has_hopf ~ abs(rF - rF_root_med),
        TRUE ~ 1
      ),
      pen_RH = ifelse(RH_ok == TRUE, 0, 1),
      pen_cx = ifelse(has_complex == TRUE, 0, 1),
      score = w$w_omega*pen_omega + w$w_sigma*pen_sigma + w$w_lambda*pen_lambda +
        w$w_hopf*pen_hopf + w$w_dist*pen_dist + w$w_RH*pen_RH + w$w_cx*pen_cx
    ) %>%
    arrange(score)
  
  write_csv_safe(df_scored, file.path(dirs$stage4, "stage4_scored_candidates.csv"))
  
  gate_diag <- tibble(
    n_total = nrow(df_scored),
    n_RH_ok = sum(df_scored$RH_ok == TRUE, na.rm = TRUE),
    n_complex = sum(df_scored$has_complex == TRUE, na.rm = TRUE),
    n_hopf = sum(df_scored$has_hopf == TRUE, na.rm = TRUE),
    n_omega_le_cap = sum(df_scored$omega_star <= s4$gate_omega_max, na.rm = TRUE)
  )
  write_csv_safe(gate_diag, file.path(dirs$stage4, "gate_diagnostics.csv"))
  print(gate_diag)
  
  df_gated_strict <- df_scored %>%
    filter(if (gates$gate_RH_required) RH_ok == TRUE else TRUE) %>%
    filter(if (gates$gate_complex_required) has_complex == TRUE else TRUE) %>%
    filter(if (gates$gate_hopf_required) has_hopf == TRUE else TRUE) %>%
    filter(is.finite(omega_star) & omega_star <= s4$gate_omega_max)
  
  df_gated_relaxed <- df_scored %>%
    filter(is.finite(omega_star) & omega_star <= s4$gate_omega_max)
  
  write_csv_safe(df_gated_strict,  file.path(dirs$stage4, "stage4_scored_candidates_GATED_STRICT.csv"))
  write_csv_safe(df_gated_relaxed, file.path(dirs$stage4, "stage4_scored_candidates_GATED_RELAXED.csv"))
  
  cat("Stage 4 gated sizes:\n")
  cat("  strict  =", nrow(df_gated_strict), "\n")
  cat("  relaxed =", nrow(df_gated_relaxed), "\n")
  
  # Top50 exports
  pick_cols <- c(
    "try_id","cand_id",
    "sigma","g_n","i","delta","kappa_max",
    "psi","phi2","rF",
    "omega_star","d_star",
    "lambda_star","f_star","Z_star",
    "stable","RH_ok","has_complex",
    "has_hopf","n_roots","rF_root_med","rF_root_p25","rF_root_p75",
    "pen_omega","pen_sigma","pen_lambda","pen_hopf","pen_dist",
    "score"
  )
  write_csv_safe(df_scored        %>% select(any_of(pick_cols)) %>% slice_head(n = 50),
                 file.path(dirs$stage4, "stage4_shortlist_top50_ALL.csv"))
  write_csv_safe(df_gated_relaxed %>% select(any_of(pick_cols)) %>% slice_head(n = 50),
                 file.path(dirs$stage4, "stage4_shortlist_top50_GATED_RELAXED.csv"))
  write_csv_safe(df_gated_strict  %>% select(any_of(pick_cols)) %>% slice_head(n = 50),
                 file.path(dirs$stage4, "stage4_shortlist_top50_GATED_STRICT.csv"))
  
  # Plots
  df_plot <- df_scored %>% mutate(score_decile = ntile(score, 10) %>% as.factor())
  
  p1 <- ggplot(df_plot, aes(x = sigma, y = omega_star, color = score_decile)) +
    geom_point(alpha = 0.75) +
    geom_hline(yintercept = s4$omega_pref, linetype = "dashed") +
    labs(title = "Stage 4: omega* vs sigma (score deciles)", x = "sigma", y = "omega*") +
    theme_minimal()
  ggsave(file.path(dirs$stage4, "score_scatter_omega_sigma.png"), p1, width = 10, height = 6, dpi = 160)
  
  p2 <- ggplot(df_plot %>% filter(has_hopf), aes(x = rF_root_med, y = phi2, color = score_decile)) +
    geom_point(alpha = 0.75) +
    facet_grid(psi ~ kappa_max) +
    geom_vline(xintercept = s4$rF_root_pref, linetype = "dashed") +
    labs(title = "Stage 4: Hopf boundary locations (rF_root_med vs phi2)",
         x = "rF_root (median)", y = "phi2") +
    theme_minimal()
  ggsave(file.path(dirs$stage4, "score_hopf_boundary_phi2.png"), p2, width = 12, height = 8, dpi = 160)
  
  p3 <- ggplot(df_plot, aes(x = score)) +
    geom_histogram(bins = 50, alpha = 0.75) +
    labs(title = "Stage 4: score distribution", x = "score (lower=better)", y = "count") +
    theme_minimal()
  ggsave(file.path(dirs$stage4, "score_hist.png"), p3, width = 10, height = 5, dpi = 160)
  
  p4 <- ggplot(df_plot, aes(x = omega_star)) +
    geom_histogram(bins = 40, alpha = 0.75) +
    geom_vline(xintercept = s4$omega_pref, linetype = "dashed") +
    labs(title = "Stage 4: omega* distribution", x = "omega*", y = "count") +
    theme_minimal()
  ggsave(file.path(dirs$stage4, "omega_hist.png"), p4, width = 10, height = 5, dpi = 160)
  
  # LaTeX shortlist (top 20 relaxed)
  short_tex <- df_gated_relaxed %>%
    select(try_id, cand_id, psi, kappa_max, sigma, omega_star, rF, has_hopf, rF_root_med, score) %>%
    slice_head(n = 20) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    mutate(has_hopf = ifelse(has_hopf == TRUE, "T", "F"))
  
  cols  <- names(short_tex)
  align <- paste(rep("l", length(cols)), collapse = "")
  header <- paste(cols, collapse = " & ")
  rows <- apply(short_tex, 1, function(r) paste(r, collapse = " & "))
  rows <- paste0(rows, " \\\\")
  
  tex_lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{Stage 4 shortlist (top 20 by score; rounded to 2 decimals).}",
    "\\label{tab:stage4_shortlist_top20}",
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(tex_lines, file.path(dirs$stage4, "shortlist_top20.tex"))
  
  list(df_scored = df_scored, df_gated_relaxed = df_gated_relaxed, df_gated_strict = df_gated_strict)
}

# ============================================================
# 9) STAGE 5: SIMULATION (hook-based), graceful or strict
# ============================================================

# ============================================================
# R/hooks/simulate_model.R
# Hook for Stage 5
# Must define: simulate_model(row, rF_sim, t_end, dt, x0=NULL)
# Returns a tibble with columns: time, e, omega, d
# ============================================================

simulate_model <- function(row, rF_sim, t_end, dt, x0 = NULL) {
  if (!requireNamespace("deSolve", quietly = TRUE)) {
    stop("simulate_model() requires the 'deSolve' package. Install it via install.packages('deSolve').")
  }
  
  # row arrives as a 1-row tibble/data.frame from pmap; coerce safely
  if (is.data.frame(row)) row <- as.list(row[1, ])
  if (!is.list(row)) stop("simulate_model(): 'row' must be list-like (a 1-row tibble/data.frame is fine).")
  
  # ---- Helpers (self-contained) ----
  logistic <- function(x) 1 / (1 + exp(-x))
  
  # Pull global defaults if you have them, else fallback
  # These match the par_base you showed earlier.
  p <- if (exists("par_base", inherits = TRUE)) get("par_base", inherits = TRUE) else list(
    kappa_min = 0.02,
    kappa_max = 0.25,
    kappa0    = 0.10,
    kappa1    = 30.0,
    phi3      = 8.0,
    phi4      = 1.0,
    phi0      = -0.02,
    alpha     = 0.02,
    phi1_min  = 0.10,
    phi1_max  = 5.00
  )
  
  # Candidate overrides
  if (!is.null(row$kappa_max) && is.finite(row$kappa_max)) p$kappa_max <- as.numeric(row$kappa_max)
  if (!is.null(row$kappa_min) && is.finite(row$kappa_min)) p$kappa_min <- as.numeric(row$kappa_min)
  
  # Core candidate parameters (must exist)
  sigma <- as.numeric(row$sigma)
  g_n   <- as.numeric(row$g_n)
  i     <- as.numeric(row$i)
  delta <- as.numeric(row$delta)
  
  psi  <- as.numeric(row$psi)
  phi2 <- as.numeric(row$phi2)
  
  stopifnot(is.finite(sigma), is.finite(g_n), is.finite(i), is.finite(delta),
            is.finite(psi), is.finite(phi2), g_n > 0, sigma > 0)
  
  kappa_fun <- function(r) {
    p$kappa_min + (p$kappa_max - p$kappa_min) * logistic(p$kappa1 * (r - p$kappa0))
  }
  
  Z_fun <- function(d, f) {
    logistic(p$phi3 * ((d - 1) + p$phi4 * (f - 1)))
  }
  
  lambda_fun <- function(r, rF, psi) logistic(psi * (r - rF))
  
  # ---- Choose phi1 (constant) ----
  # Prefer Stage2's phi1_endog if present, else reconstruct from stored Z_star if possible.
  e_target <- if (exists("cfg", inherits = TRUE) && !is.null(get("cfg", inherits=TRUE)$targets$e_target)) {
    get("cfg", inherits = TRUE)$targets$e_target
  } else {
    0.94
  }
  
  phi1 <- NA_real_
  if (!is.null(row$phi1_endog) && is.finite(row$phi1_endog)) {
    phi1 <- as.numeric(row$phi1_endog)
  } else if (!is.null(row$Z_star) && is.finite(row$Z_star)) {
    phi1 <- (p$alpha - p$phi0 + phi2 * as.numeric(row$Z_star)) / e_target
  } else {
    # last resort: use a sane mid value; dynamics still run
    phi1 <- 1.0
  }
  if (!is.finite(phi1)) phi1 <- 1.0
  
  # ---- Initial conditions ----
  if (is.null(x0)) {
    # Use steady state if available, else defaults
    e0 <- if (!is.null(row$e_star) && is.finite(row$e_star)) as.numeric(row$e_star) else e_target
    w0 <- if (!is.null(row$omega_star) && is.finite(row$omega_star)) as.numeric(row$omega_star) else 0.65
    d0 <- if (!is.null(row$d_star) && is.finite(row$d_star)) as.numeric(row$d_star) else 0.5
    x0 <- c(e = e0, omega = w0, d = d0)
  } else {
    x0 <- as.numeric(x0)
    names(x0) <- names(x0) %||% c("e", "omega", "d")
    if (!all(c("e","omega","d") %in% names(x0))) stop("x0 must be a named vector with e, omega, d.")
  }
  
  # ---- ODE system (reduced 3D) ----
  rhs <- function(t, state, parms) {
    e     <- state[["e"]]
    omega <- state[["omega"]]
    d     <- state[["d"]]
    
    # implied profit rate
    r <- (1 - omega - i * d) / sigma
    
    # investment / growth
    kappa <- kappa_fun(r)
    g     <- kappa / sigma - delta
    
    # financialization block
    lam <- lambda_fun(r, rF = rF_sim, psi = psi)
    # avoid division blow-ups if lam ~ 1
    lam <- min(max(lam, 1e-10), 1 - 1e-10)
    
    iotaF <- r * (lam / (1 - lam))
    f     <- iotaF / g_n
    
    # discipline
    Z <- Z_fun(d, f)
    
    de     <- (g - g_n) * e
    domega <- omega * (p$phi0 + phi1 * e - p$alpha - phi2 * Z)
    dd     <- kappa - (1 - omega) + i * d - d * g
    
    list(c(de, domega, dd))
  }
  
  times <- seq(0, t_end, by = dt)
  
  sol <- deSolve::ode(
    y = x0,
    times = times,
    func = rhs,
    parms = NULL,
    method = "lsoda"
  )
  
  sol <- as.data.frame(sol)
  tibble::as_tibble(sol) %>%
    dplyr::rename(time = .data$time) %>%
    dplyr::select(time, e, omega, d)
}





call_simulate_model <- function(row, rF_sim, t_end, dt, x0 = NULL) {
  # supports either simulate_model(row, rF_sim, t_end, dt)
  # or simulate_model(row, rF_sim, t_end, dt, x0)
  fm <- names(formals(simulate_model))
  if (!is.null(x0) && "x0" %in% fm) {
    simulate_model(row, rF_sim = rF_sim, t_end = t_end, dt = dt, x0 = x0)
  } else {
    simulate_model(row, rF_sim = rF_sim, t_end = t_end, dt = dt)
  }
}

cycle_metrics <- function(sim, burn_in, var_min, amp_min) {
  sim2 <- as_tibble(sim) %>% filter(time >= burn_in)
  if (!("omega" %in% names(sim2))) return(tibble(has_cycle = FALSE, reason = "missing_omega"))
  
  x <- sim2$omega
  t <- sim2$time
  v <- var(x, na.rm = TRUE)
  a <- (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) / 2
  
  if (!is.finite(v) || !is.finite(a)) return(tibble(has_cycle = FALSE, reason = "nonfinite", var=v, amp=a))
  if (v < var_min || a < amp_min) return(tibble(has_cycle = FALSE, reason = "low_var_or_amp", var=v, amp=a))
  
  dx <- diff(x)
  peaks <- which(head(dx, -1) > 0 & tail(dx, -1) <= 0) + 1
  if (length(peaks) < 3) return(tibble(has_cycle = TRUE, reason = "few_peaks", var=v, amp=a, period=NA_real_))
  
  periods <- diff(t[peaks])
  tibble(has_cycle = TRUE, reason = "ok", var=v, amp=a,
         period = mean(periods, na.rm = TRUE),
         omega_mean = mean(x, na.rm = TRUE))
}

bounded_ok <- function(sim) {
  sim <- as_tibble(sim)
  ok <- TRUE
  if ("e" %in% names(sim))     ok <- ok && all(is.finite(sim$e))     && all(sim$e >= 0)     && all(sim$e <= 1.5)
  if ("omega" %in% names(sim)) ok <- ok && all(is.finite(sim$omega)) && all(sim$omega >= 0) && all(sim$omega <= 1.5)
  if ("d" %in% names(sim))     ok <- ok && all(is.finite(sim$d))     && all(sim$d >= -0.1)  && all(sim$d <= 10)
  ok
}

run_stage5 <- function(cfg, dirs) {
  cat("\n====================\nStage 5: simulation\n====================\n")
  set.seed(cfg$stage5$seed_stage5)
  
  # try to source hook if missing
  if (!exists("simulate_model")) {
    if (file.exists(cfg$stage5$hook_path)) source(cfg$stage5$hook_path)
  }
  
  if (!exists("simulate_model")) {
    msg <- "Stage 5 requires simulate_model(row, rF_sim, t_end, dt) -> columns: time,e,omega,d"
    if (isTRUE(cfg$stage5$require_simulate_hook)) stop(msg)
    cat("Skipping Stage 5 (hook missing): ", msg, "\n", sep="")
    write_csv_safe(tibble(), file.path(dirs$stage5B, "stage5B_cycle_metrics.csv"))
    write_csv_safe(tibble(), file.path(dirs$stage5C, "robustness_summary.csv"))
    write_csv_safe(tibble(), file.path(dirs$stage5D, "stage5_final_ranked.csv"))
    return(invisible(NULL))
  }
  
  # Candidate pool (prefer Stage 4 relaxed shortlist)
  cand_pool_file <- file.path(dirs$stage4, "stage4_shortlist_top50_GATED_RELAXED.csv")
  if (!file.exists(cand_pool_file)) cand_pool_file <- file.path(dirs$stage4, "stage4_scored_candidates_GATED_RELAXED.csv")
  if (!file.exists(cand_pool_file)) cand_pool_file <- file.path(dirs$stage4, "stage4_scored_candidates.csv")
  stopifnot(file.exists(cand_pool_file))
  
  cand <- read_csv(cand_pool_file, show_col_types = FALSE) %>%
    slice_head(n = cfg$stage5$n_candidates) %>%
    mutate(candidate_id = row_number())
  
  write_csv_safe(cand, file.path(dirs$stage5A, "stage5_candidate_pool.csv"))
  
  # Root seed selection: rF_root_med if present else rF
  cand_5A <- cand %>%
    mutate(
      rF_seed = case_when(
        "rF_root_med" %in% names(.) & is.finite(rF_root_med) ~ as.numeric(rF_root_med),
        TRUE ~ as.numeric(rF)
      ),
      rF_root = rF_seed
    )
  write_csv_safe(cand_5A, file.path(dirs$stage5A, "stage5A_hopf_certified.csv"))
  
  # Sims below/above root
  sim_5B <- cand_5A %>%
    mutate(
      sim = pmap(., function(...) {
        row <- tibble(...)
        id <- as.character(row$candidate_id)
        
        rF_root <- as.numeric(row$rF_root)
        if (!is.finite(rF_root)) rF_root <- as.numeric(row$rF_seed)
        
        rF_lo <- max(1e-6, rF_root - cfg$stage5$eps_root)
        rF_hi <- rF_root + cfg$stage5$eps_root
        
        # initial condition suggestion (optional)
        x0 <- c(e = 0.90, omega = 0.65, d = 0.5) + cfg$stage5$perturb
        if (all(c("e_star","omega_star","d_star") %in% names(row))) {
          x0 <- c(e = row$e_star, omega = row$omega_star, d = row$d_star) + cfg$stage5$perturb
        } else if (all(c("omega_star","d_star") %in% names(row))) {
          x0 <- c(e = cfg$targets$e_target, omega = row$omega_star, d = row$d_star) + cfg$stage5$perturb
        }
        
        sim_lo <- tryCatch(call_simulate_model(row, rF_sim = rF_lo, t_end = cfg$stage5$t_end, dt = cfg$stage5$dt, x0 = x0),
                           error = function(e) NULL)
        sim_hi <- tryCatch(call_simulate_model(row, rF_sim = rF_hi, t_end = cfg$stage5$t_end, dt = cfg$stage5$dt, x0 = x0),
                           error = function(e) NULL)
        
        if (is.null(sim_lo) || is.null(sim_hi)) {
          return(tibble(
            sim_ok = FALSE, sim_reason = "simulate_model_error",
            bounded_below = NA, bounded_above = NA,
            has_cycle_below = NA, has_cycle_above = NA,
            amp_below = NA_real_, amp_above = NA_real_,
            period_below = NA_real_, period_above = NA_real_,
            omega_mean_above = NA_real_,
            rF_lo = rF_lo, rF_hi = rF_hi
          ))
        }
        
        sim_lo <- as_tibble(sim_lo)
        sim_hi <- as_tibble(sim_hi)
        
        write_csv_safe(sim_lo, file.path(dirs$stage5B, paste0("traj_", id, "_below.csv")))
        write_csv_safe(sim_hi, file.path(dirs$stage5B, paste0("traj_", id, "_above.csv")))
        
        met_lo <- cycle_metrics(sim_lo, burn_in = cfg$stage5$burn_in,
                                var_min = cfg$stage5$cycle_var_min, amp_min = cfg$stage5$amp_min)
        met_hi <- cycle_metrics(sim_hi, burn_in = cfg$stage5$burn_in,
                                var_min = cfg$stage5$cycle_var_min, amp_min = cfg$stage5$amp_min)
        
        tibble(
          sim_ok = TRUE, sim_reason = "ok",
          bounded_below = bounded_ok(sim_lo),
          bounded_above = bounded_ok(sim_hi),
          has_cycle_below = met_lo$has_cycle %||% FALSE,
          has_cycle_above = met_hi$has_cycle %||% FALSE,
          amp_below = met_lo$amp %||% NA_real_,
          amp_above = met_hi$amp %||% NA_real_,
          period_below = met_lo$period %||% NA_real_,
          period_above = met_hi$period %||% NA_real_,
          omega_mean_above = met_hi$omega_mean %||% NA_real_,
          rF_lo = rF_lo, rF_hi = rF_hi
        )
      })
    ) %>%
    unnest(sim)
  
  write_csv_safe(sim_5B, file.path(dirs$stage5B, "stage5B_cycle_metrics.csv"))
  
  # Robustness ring
  targets <- sim_5B %>%
    filter(sim_ok == TRUE, bounded_above == TRUE, has_cycle_above == TRUE)
  
  robust_5C <- targets %>%
    mutate(
      robust = pmap(., function(...) {
        row <- tibble(...)
        
        rF_test <- as.numeric(row$rF_root)
        if (!is.finite(rF_test)) rF_test <- as.numeric(row$rF_seed)
        rF_test <- rF_test + cfg$stage5$eps_root
        
        out <- map_dfr(seq_len(cfg$stage5$robust_reps), function(rep_id) {
          rowj <- row
          
          for (nm in cfg$stage5$jitter_cols) {
            if (nm %in% names(rowj) && is.finite(rowj[[nm]])) {
              rowj[[nm]] <- rowj[[nm]] * exp(rnorm(1, 0, cfg$stage5$jitter_sd))
            }
          }
          
          sim <- tryCatch(call_simulate_model(rowj, rF_sim = rF_test, t_end = cfg$stage5$robust_t_end, dt = cfg$stage5$dt),
                          error=function(e) NULL)
          if (is.null(sim)) {
            return(tibble(rep = rep_id, ok = FALSE, bounded = FALSE, has_cycle = FALSE, amp = NA_real_, period = NA_real_))
          }
          
          sim <- as_tibble(sim)
          met <- cycle_metrics(sim, burn_in = min(cfg$stage5$burn_in, cfg$stage5$robust_t_end/2),
                               var_min = cfg$stage5$cycle_var_min, amp_min = cfg$stage5$amp_min)
          
          tibble(
            rep = rep_id,
            ok = TRUE,
            bounded = bounded_ok(sim),
            has_cycle = met$has_cycle %||% FALSE,
            amp = met$amp %||% NA_real_,
            period = met$period %||% NA_real_
          )
        })
        
        tibble(
          survival_rate = mean(out$bounded & out$has_cycle, na.rm = TRUE),
          bounded_rate  = mean(out$bounded, na.rm = TRUE),
          cycle_rate    = mean(out$has_cycle, na.rm = TRUE),
          amp_med       = median(out$amp, na.rm = TRUE),
          period_med    = median(out$period, na.rm = TRUE)
        )
      })
    ) %>%
    unnest(robust)
  
  write_csv_safe(robust_5C, file.path(dirs$stage5C, "robustness_summary.csv"))
  
  # Final ranking
  final_5D <- sim_5B %>%
    left_join(robust_5C %>% select(candidate_id, survival_rate, bounded_rate, cycle_rate, amp_med, period_med),
              by = "candidate_id") %>%
    mutate(
      survival_rate = coalesce(survival_rate, 0),
      bounded_rate  = coalesce(bounded_rate, 0),
      cycle_rate    = coalesce(cycle_rate, 0),
      score_stage5 = 3*survival_rate + 1*bounded_rate + 1*cycle_rate +
        0.2*if_else(is.finite(amp_above), pmin(amp_above, 1), 0)
    ) %>%
    arrange(desc(score_stage5))
  
  write_csv_safe(final_5D, file.path(dirs$stage5D, "stage5_final_ranked.csv"))
  write_csv_safe(final_5D %>% slice_head(n = 20), file.path(dirs$stage5D, "stage5_top20.csv"))
  
  cat("\nStage 5 complete.\n")
  cat("  5A: ", file.path(dirs$stage5A, "stage5A_hopf_certified.csv"), "\n", sep="")
  cat("  5B: ", file.path(dirs$stage5B, "stage5B_cycle_metrics.csv"), "\n", sep="")
  cat("  5C: ", file.path(dirs$stage5C, "robustness_summary.csv"), "\n", sep="")
  cat("  5D: ", file.path(dirs$stage5D, "stage5_final_ranked.csv"), "\n", sep="")
  
  invisible(final_5D)
}

# ============================================================
# 10) RUN PIPELINE (sequential, restart-friendly)
# ============================================================
out1 <- NULL; out2 <- NULL; out3 <- NULL; out4 <- NULL

if (isTRUE(cfg$run$stage1)) {
  out1 <- run_stage1(cfg, dirs, par_base)
}

if (isTRUE(cfg$run$stage2)) {
  stage1_ok <- if (!is.null(out1)) out1$stage1_ok else {
    read_csv(file.path(dirs$stage1, "stage1_backbone.csv"), show_col_types = FALSE) %>% filter(backbone_ok)
  }
  out2 <- run_stage2(cfg, dirs, par_base, stage1_ok)
}

if (isTRUE(cfg$run$stage3)) {
  stage2_ok <- if (!is.null(out2)) out2$stage2_ok else {
    read_csv(file.path(dirs$stage2, "stage2_finance_scan.csv"), show_col_types = FALSE) %>% filter(econ_ok)
  }
  out3 <- run_stage3(cfg, dirs, par_base, stage2_ok)
}

if (isTRUE(cfg$run$stage4)) {
  out4 <- run_stage4(cfg, dirs)
}

if (isTRUE(cfg$run$stage5)) {
  run_stage5(cfg, dirs)
}
