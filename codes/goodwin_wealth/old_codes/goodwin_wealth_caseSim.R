############################################################
## wealth_goodwin_wealth_inequality_RUNALL.R
##
## Purpose:
##  - Reduced-system Hopf mapping (SS -> J -> RH -> H) across planes
##  - Rank planes + export Hopf maps with sign(H)
##  - Find lambda* on a chosen slice and simulate scenarios:
##      (i) fixed lambda below/near/above lambda*
##      (ii) moving-Hopffiness “patriarch” scenario (lambda drift)
##  - Wealth block: asset price growth + inequality ratchet
##  - NUMERICAL STABILITY FIXES:
##      * evolve log_pA instead of pA
##      * cap debt_stress used in p_g
##      * cap p_g (pos/neg)
##      * stop integration via rootfun (|d| and log_pA caps)
##  - Exports: CSV tables + plots + audit bundle + manifest
############################################################

suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(fs)
  library(tibble)
  library(grid)   # unit() for arrows
  library(readr)
})

## ============================================================
## GLOBAL OUTPUT DIR
## ============================================================
OUT_DIR <- fs::path("outputs", "wealth_goodwin", "wealth_inequality")
fs::dir_create(OUT_DIR)

write_csv_quiet <- function(x, path) readr::write_csv(x, path)

## ============================================================
## 0) Helpers
## ============================================================
inv_logit <- function(x) 1 / (1 + exp(-x))

require_params <- function(par, keys) {
  miss <- setdiff(keys, names(par))
  if (length(miss) > 0) stop("Missing parameters: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || !is.finite(a[1])) b else a

clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

pos_part <- function(x) pmax(x, 0)

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

period_peaks <- function(t, x) {
  x <- as.numeric(x); t <- as.numeric(t)
  ok <- is.finite(x) & is.finite(t)
  x <- x[ok]; t <- t[ok]
  if (length(x) < 200) return(NA_real_)
  dx1 <- diff(x)
  peaks <- which(dx1[-1] < 0 & dx1[-length(dx1)] > 0) + 1
  if (length(peaks) < 3) return(NA_real_)
  median(diff(t[peaks]))
}

## ============================================================
## 1) Core reduced model pieces (u fixed at u_star for Hopf map)
## ============================================================
Phi_fun   <- function(e, par) par$phi0 + par$phi1 * e
Phi_e_fun <- function(e, par) par$phi1
pi_fun    <- function(omega) 1 - omega

r_net_reduced <- function(omega, d, par) {
  (par$u_star / par$sigma) * (pi_fun(omega) - par$i * d)
}

kappa_fun <- function(r, par) {
  par$kappa_min + (par$kappa_max - par$kappa_min) /
    (1 + exp(-par$lambda * (r - par$r0)))
}

kappa_r_fun <- function(r, par) {
  s <- inv_logit(par$lambda * (r - par$r0))
  (par$kappa_max - par$kappa_min) * par$lambda * s * (1 - s)
}

## Optional leakage psi(d)
psi_fun <- function(d, par) {
  if (isTRUE(par$psi_zero)) return(0)
  inv_logit(par$psi_lambda * (d - par$psi_mid)) * par$psi_max
}

## Owners saving (here usually constant)
sA_fun <- function(chiA, par) {
  if (isTRUE(par$sA_const)) return(par$sA0)
  s <- inv_logit(par$sA_lambda * (chiA - par$sA_mid))
  par$sA_min + (par$sA_max - par$sA_min) * s
}

## Non-owners induced consumption “share”
cN_ind_fun <- function(yN, par) {
  if (yN <= 0) return(par$cN_minus)
  par$cN_minus + (par$cN_plus - par$cN_minus) * (yN / (yN + par$cN_eta))
}

## Autonomous baseline consumption share
cN_bar_fun <- function(t, par) par$varsigma0 * exp(par$varsigma1 * t)

## Demand wedge Delta evaluated at u = u_star
Delta_fun_uStar <- function(t, omega, d, par) {
  u <- par$u_star
  r   <- (u / par$sigma) * (pi_fun(omega) - par$i * d)
  kap <- kappa_fun(r, par)
  
  rho  <- par$sigma / u
  chiA <- par$p * rho + d
  sA   <- sA_fun(chiA, par)
  
  yN     <- omega - par$i * d
  cN_ind <- cN_ind_fun(yN, par)
  
  psi <- psi_fun(d, par)
  cA  <- (1 - sA) * (1 - omega)
  
  1 - psi - cN_ind - cA - kap
}

u_desired_uStar <- function(t, omega, d, par) {
  Delta <- Delta_fun_uStar(t, omega, d, par)
  if (!is.finite(Delta) || Delta <= 1e-10) return(NA_real_)
  cbar <- cN_bar_fun(t, par)
  cbar / Delta
}

## ============================================================
## 2) Reduced steady state (u = u_star)
## ============================================================
steady_state_reduced <- function(par, eps = 1e-12) {
  
  require_params(par, c(
    "alpha","beta","delta","sigma","u_star","i",
    "kappa_min","kappa_max","lambda","r0",
    "phi0","phi1"
  ))
  
  g_n <- par$alpha + par$beta
  if (g_n <= 0) return(list(ok = FALSE, reason = "g_n <= 0 (alpha+beta must be >0)"))
  
  A <- par$u_star / par$sigma
  if (A <= 0) return(list(ok = FALSE, reason = "A=u_star/sigma must be > 0"))
  
  ## Balanced growth at u=u_star:
  ## gY = A*kappa - delta = g_n  => kappa_target
  k_target <- (g_n + par$delta) / A
  
  if (k_target <= par$kappa_min + eps || k_target >= par$kappa_max - eps) {
    return(list(ok = FALSE, reason = "kappa_target outside (kappa_min, kappa_max)"))
  }
  
  s <- (k_target - par$kappa_min) / (par$kappa_max - par$kappa_min)
  if (s <= eps || s >= 1 - eps) return(list(ok = FALSE, reason = "logit inversion s not in (0,1)"))
  
  r_star <- par$r0 + (1 / par$lambda) * log(s / (1 - s))
  
  ## d* (your identity)
  d_star <- (k_target - (1 / A) * r_star) / g_n
  if (!is.finite(d_star)) return(list(ok = FALSE, reason = "d* non-finite"))
  ## NOTE: do NOT force d*>0 (you explicitly allow sign); but reject absurdly close-to-zero denominators
  if (abs(d_star) < eps) d_star <- 0
  
  ## omega* from dd=0: 1-omega = kappa + (i-g_n)d
  pi_star <- k_target + (par$i - g_n) * d_star
  omega_star <- 1 - pi_star
  if (omega_star <= eps || omega_star >= 1 - eps) return(list(ok = FALSE, reason = "omega* not in (0,1)"))
  
  ## e* from Phillips: Phi(e*) = alpha
  if (par$phi1 <= 0) return(list(ok = FALSE, reason = "phi1 must be >0"))
  e_star <- (par$alpha - par$phi0) / par$phi1
  if (e_star <= eps || e_star >= 1 - eps) return(list(ok = FALSE, reason = "e* not in (0,1)"))
  
  list(
    ok = TRUE,
    e = e_star,
    omega = omega_star,
    d = d_star,
    r = r_star,
    kappa = k_target,
    g_n = g_n,
    A = A
  )
}

## ============================================================
## 3) Enforce u_d = u_star at SS by calibrating varsigma0
## ============================================================
calibrate_varsigma0_at_ss <- function(par, t0 = 0) {
  require_params(par, c("u_star","varsigma1","varsigma0"))
  
  if (abs(par$varsigma1) > 0) {
    warning("varsigma1 != 0: calibration enforces u_d=u_star only at t0.")
  }
  
  ss <- steady_state_reduced(par)
  if (!isTRUE(ss$ok)) return(list(ok = FALSE, reason = ss$reason, par = par))
  
  Delta_star <- Delta_fun_uStar(t0, ss$omega, ss$d, par)
  if (!is.finite(Delta_star) || Delta_star <= 1e-10) {
    return(list(ok = FALSE, reason = "Delta* nonpositive at SS; cannot calibrate varsigma0", par = par))
  }
  
  par2 <- par
  par2$varsigma0 <- par$u_star * Delta_star
  
  list(ok = TRUE, par = par2, ss = ss, Delta_star = Delta_star)
}

## ============================================================
## 4) Reduced Jacobian (e, omega, d) at SS + RH/Hopf functional
## ============================================================
jacobian_reduced_at_ss <- function(ss, par) {
  e <- ss$e; omega <- ss$omega; d <- ss$d
  g_n <- ss$g_n
  A <- ss$A
  r <- ss$r
  
  kap_r <- kappa_r_fun(r, par)
  Phi_e <- Phi_e_fun(e, par)
  
  J <- matrix(0, 3, 3)
  
  ## de
  J[1,2] <- - e * (A^2) * kap_r
  J[1,3] <- - e * (A^2) * par$i * kap_r
  
  ## domega
  J[2,1] <- omega * Phi_e
  
  ## dd
  J[3,2] <- 1 - A * kap_r + d * (A^2) * kap_r
  J[3,3] <- (par$i - g_n) - A * par$i * kap_r + d * (A^2) * par$i * kap_r
  
  list(J = J, kappa_r = kap_r, Phi_e = Phi_e)
}

rh_hopf_reduced <- function(J) {
  trJ  <- sum(diag(J))
  M2 <- (J[1,1]*J[2,2] - J[1,2]*J[2,1]) +
    (J[1,1]*J[3,3] - J[1,3]*J[3,1]) +
    (J[2,2]*J[3,3] - J[2,3]*J[3,2])
  detJ <- det(J)
  
  a1 <- -trJ
  a2 <- M2
  a3 <- -detJ
  
  H <- a1 * a2 - a3
  feasible_RH <- (a1 > 0) && (a2 > 0) && (a3 > 0)
  
  list(a1 = a1, a2 = a2, a3 = a3, H = H, feasible_RH = feasible_RH)
}

analyze_reduced <- function(par, calibrate_varsigma0 = TRUE) {
  if (calibrate_varsigma0) {
    cal <- calibrate_varsigma0_at_ss(par)
    if (!isTRUE(cal$ok)) return(list(ok = FALSE, reason = cal$reason))
    par <- cal$par
  }
  
  ss <- steady_state_reduced(par)
  if (!isTRUE(ss$ok)) return(list(ok = FALSE, reason = ss$reason))
  
  Jp <- jacobian_reduced_at_ss(ss, par)
  rh <- rh_hopf_reduced(Jp$J)
  
  list(ok = TRUE, par = par, steady_state = ss, jacobian = Jp$J, rh = rh)
}

hopf_gap_reduced <- function(par, calibrate_varsigma0 = TRUE) {
  out <- analyze_reduced(par, calibrate_varsigma0 = calibrate_varsigma0)
  if (!isTRUE(out$ok)) return(NA_real_)
  out$rh$H
}

find_hopf_mu <- function(par, mu, lo, hi,
                         n_grid = 350, tol = 1e-8, maxit = 200,
                         calibrate_varsigma0 = TRUE, verbose = TRUE) {
  stopifnot(mu %in% names(par))
  
  f <- function(x) {
    par2 <- par
    par2[[mu]] <- x
    hopf_gap_reduced(par2, calibrate_varsigma0 = calibrate_varsigma0)
  }
  
  xs <- seq(lo, hi, length.out = n_grid)
  Hs <- vapply(xs, f, numeric(1))
  
  ok <- is.finite(Hs)
  if (!any(ok)) stop("No admissible points for H on [lo,hi].")
  
  xs_ok <- xs[ok]; Hs_ok <- Hs[ok]
  
  if (verbose) {
    cat(sprintf("\n[scan] admissible points: %d / %d\n", length(xs_ok), length(xs)))
    cat(sprintf("[scan] H range: [%.4g, %.4g]\n", min(Hs_ok), max(Hs_ok)))
  }
  
  sgn <- sign(Hs_ok)
  idx <- which(sgn[-1] * sgn[-length(sgn)] < 0)
  if (length(idx) == 0) stop("No sign change in H on admissible subset. Expand bracket or change plane.")
  
  a <- xs_ok[idx[1]]
  b <- xs_ok[idx[1] + 1]
  if (verbose) cat(sprintf("[scan] bracketing H=0 on [%g, %g]\n", a, b))
  
  uniroot(f, interval = c(a, b), tol = tol, maxiter = maxit)$root
}

## ============================================================
## 5) Hopf mapping on a plane + scoring + plotting
## ============================================================
grid_values <- function(lo, hi, n) seq(lo, hi, length.out = n)

hopf_plane <- function(par_base, p1, p2, p1_vals, p2_vals,
                       calibrate_varsigma0 = TRUE, H_band = 1e-3) {
  stopifnot(p1 %in% names(par_base), p2 %in% names(par_base))
  
  res <- vector("list", length(p1_vals) * length(p2_vals))
  k <- 0
  
  for (v1 in p1_vals) {
    for (v2 in p2_vals) {
      k <- k + 1
      par <- par_base
      par[[p1]] <- v1
      par[[p2]] <- v2
      
      out <- analyze_reduced(par, calibrate_varsigma0 = calibrate_varsigma0)
      
      if (!isTRUE(out$ok)) {
        res[[k]] <- tibble(
          !!p1 := v1, !!p2 := v2,
          ok = FALSE,
          reason = out$reason,
          a1 = NA_real_, a2 = NA_real_, a3 = NA_real_, H = NA_real_,
          feasible_RH = FALSE,
          hopf_band = FALSE
        )
      } else {
        rh <- out$rh
        res[[k]] <- tibble(
          !!p1 := v1, !!p2 := v2,
          ok = TRUE,
          reason = NA_character_,
          a1 = rh$a1, a2 = rh$a2, a3 = rh$a3, H = rh$H,
          feasible_RH = rh$feasible_RH,
          hopf_band = isTRUE(rh$feasible_RH) && is.finite(rh$H) && (abs(rh$H) < H_band)
        )
      }
    }
  }
  bind_rows(res)
}

score_hopf_plane <- function(df, H_band = 1e-3, tol = 1e-6) {
  dfF <- df %>% filter(ok, feasible_RH, is.finite(H))
  if (nrow(dfF) == 0) {
    return(tibble(
      n_ok = sum(df$ok, na.rm = TRUE),
      n_feasible = 0,
      feasible_share = 0,
      share_H_pos = NA_real_, share_H_neg = NA_real_,
      H_min = NA_real_, H_q05 = NA_real_, H_med = NA_real_, H_q95 = NA_real_, H_max = NA_real_,
      sd_H = NA_real_,
      sign_change = FALSE,
      has_robust_sign_change = FALSE,
      hopf_band_share = NA_real_,
      degeneracy_flag = TRUE
    ))
  }
  
  Hmin <- min(dfF$H)
  Hmax <- max(dfF$H)
  share_pos <- mean(dfF$H >  tol)
  share_neg <- mean(dfF$H < -tol)
  hopf_band_share <- mean(abs(dfF$H) < H_band)
  sdH <- sd(dfF$H)
  
  sign_change <- (Hmin < 0) && (Hmax > 0)
  has_robust_sign_change <- (Hmin < -tol) && (Hmax > tol)
  degeneracy_flag <- (sdH < 1e-8) || (!has_robust_sign_change)
  
  tibble(
    n_ok = sum(df$ok, na.rm = TRUE),
    n_feasible = nrow(dfF),
    feasible_share = nrow(dfF) / max(1, sum(df$ok, na.rm = TRUE)),
    hopf_band_share = hopf_band_share,
    share_H_pos = share_pos,
    share_H_neg = share_neg,
    H_min = Hmin,
    H_q05 = as.numeric(quantile(dfF$H, 0.05)),
    H_med = median(dfF$H),
    H_q95 = as.numeric(quantile(dfF$H, 0.95)),
    H_max = Hmax,
    sd_H = sdH,
    sign_change = sign_change,
    has_robust_sign_change = has_robust_sign_change,
    degeneracy_flag = degeneracy_flag
  )
}

plot_hopf_plane_signH <- function(df, p1, p2, title = NULL) {
  df2 <- df %>%
    mutate(
      status = case_when(
        !ok ~ "SS infeasible",
        ok & !feasible_RH ~ "RH prereq fail",
        ok & feasible_RH & is.finite(H) & H > 0 ~ "H>0 (stable)",
        ok & feasible_RH & is.finite(H) & H < 0 ~ "H<0 (unstable)",
        ok & feasible_RH & is.finite(H) & abs(H) <= 1e-10 ~ "H≈0",
        TRUE ~ "Other"
      ),
      status = factor(status, levels = c("SS infeasible","RH prereq fail","H>0 (stable)","H<0 (unstable)","H≈0","Other"))
    )
  
  ggplot(df2, aes(x = .data[[p1]], y = .data[[p2]])) +
    geom_tile(aes(fill = status), alpha = 0.92) +
    geom_contour(aes(z = H), breaks = 0, linewidth = 0.55, na.rm = TRUE) +
    labs(
      title = title %||% paste0("Hopf plane: ", p2, " vs ", p1),
      subtitle = "Tiles: SS feasibility + RH prereqs + sign(H). Contour: H=0 (Hopf boundary candidate).",
      x = p1, y = p2, fill = NULL
    )
}

compare_planes <- function(par_base, planes, ranges, n = 80,
                           H_band = 1e-3, tol = 1e-6, calibrate_varsigma0 = TRUE) {
  out <- vector("list", length(planes))
  for (j in seq_along(planes)) {
    p1 <- planes[[j]][1]
    p2 <- planes[[j]][2]
    p1_rng <- ranges[[p1]]
    p2_rng <- ranges[[p2]]
    
    df <- hopf_plane(
      par_base, p1, p2,
      p1_vals = grid_values(p1_rng[1], p1_rng[2], n),
      p2_vals = grid_values(p2_rng[1], p2_rng[2], n),
      calibrate_varsigma0 = calibrate_varsigma0,
      H_band = H_band
    )
    
    sc <- score_hopf_plane(df, H_band = H_band, tol = tol) %>%
      mutate(p1 = p1, p2 = p2)
    out[[j]] <- sc
  }
  
  bind_rows(out) %>%
    arrange(desc(has_robust_sign_change), degeneracy_flag, desc(feasible_share), desc(sd_H))
}

make_plane_decision_outputs <- function(scores, out_dir = OUT_DIR) {
  plane_tbl <- scores %>%
    mutate(
      plane = paste0("(", p1, ", ", p2, ")"),
      score = 100*feasible_share + 25*hopf_band_share + 10*as.numeric(has_robust_sign_change)
    ) %>%
    arrange(desc(has_robust_sign_change), desc(feasible_share), desc(hopf_band_share)) %>%
    select(plane, p1, p2, has_robust_sign_change, degeneracy_flag, sign_change,
           feasible_share, hopf_band_share, sd_H, H_min, H_max, n_feasible, n_ok, score)
  
  write_csv_quiet(plane_tbl, fs::path(out_dir, "TABLE_plane_ranking_tidy.csv"))
  
  p_sel <- plane_tbl %>%
    mutate(has_robust_sign_change = as.factor(has_robust_sign_change)) %>%
    ggplot(aes(x = feasible_share, y = hopf_band_share, label = plane, shape = has_robust_sign_change)) +
    geom_point(size = 2.7) +
    geom_text(check_overlap = TRUE, nudge_y = 0.002, size = 3) +
    labs(
      title = "Plane selection: feasibility vs near-Hopf thickness",
      subtitle = "Prefer robust sign change + high feasible_share. hopf_band_share proxies boundary density.",
      x = "feasible_share (RH prereqs satisfied among ok SS points)",
      y = "hopf_band_share (|H| < band, within feasible RH)",
      shape = "robust sign-change?"
    )
  
  ggsave(fs::path(out_dir, "PLOT_plane_selection_scatter.png"),
         p_sel, width = 9.2, height = 6.0, dpi = 160)
  
  plane_tbl
}

## ============================================================
## 6) Plane H audit (distribution + “degeneracy smell test”)
## ============================================================
audit_plane_df <- function(df, plane_name) {
  dfF <- df %>% filter(ok, feasible_RH, is.finite(H))
  if (nrow(dfF) == 0) {
    return(tibble(
      plane = plane_name,
      n_ok = sum(df$ok, na.rm = TRUE),
      n_feasible = 0,
      H_min = NA_real_, H_q05 = NA_real_, H_med = NA_real_, H_q95 = NA_real_, H_max = NA_real_,
      share_absH_lt_1e_06 = NA_real_,
      share_absH_lt_1e_04 = NA_real_,
      share_absH_lt_1e_03 = NA_real_,
      share_absH_lt_1e_02 = NA_real_
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
    share_absH_lt_1e_06 = mean(abs(dfF$H) < 1e-6),
    share_absH_lt_1e_04 = mean(abs(dfF$H) < 1e-4),
    share_absH_lt_1e_03 = mean(abs(dfF$H) < 1e-3),
    share_absH_lt_1e_02 = mean(abs(dfF$H) < 1e-2)
  )
}

plot_plane_H_quantiles <- function(plane_audits, out_dir = OUT_DIR) {
  if (nrow(plane_audits) == 0) return(invisible(NULL))
  p <- plane_audits %>%
    mutate(plane = factor(plane, levels = plane)) %>%
    ggplot(aes(x = plane, y = H_med)) +
    geom_point() +
    geom_errorbar(aes(ymin = H_q05, ymax = H_q95), width = 0.15) +
    coord_flip() +
    labs(
      title = "Plane audit: H median with 5–95% interval (feasible RH region)",
      x = NULL, y = "H"
    )
  ggsave(fs::path(out_dir, "PLOT_plane_H_quantiles.png"), p, width = 8.5, height = 5.2, dpi = 160)
  p
}

## ============================================================
## 7) Regime functions for simulation (moving Hopf “patriarch”)
## ============================================================
smoothstep <- function(x) {
  x <- clamp(x, 0, 1)
  x*x*(3 - 2*x)
}

## time-varying parameter interpolator
get_par_t <- function(t, par) {
  t0 <- par$t_switch %||% 350
  tw <- par$t_width  %||% 40
  s  <- smoothstep((t - t0) / tw)
  
  lambda_t <- (par$lambda_pre %||% par$lambda) * (1 - s) + (par$lambda_post %||% par$lambda) * s
  r0_t     <- (par$r0_pre %||% par$r0)       * (1 - s) + (par$r0_post %||% par$r0)       * s
  psi_max  <- (par$psi_max_pre %||% 0)       * (1 - s) + (par$psi_max_post %||% 0)       * s
  
  ## optional: shift autonomous consumption after switch (as multiplier)
  varsigma0_t <- (par$varsigma0_pre %||% par$varsigma0) * (1 - s) +
    ((par$varsigma0_post %||% NA_real_) %||% (par$varsigma0_pre %||% par$varsigma0)) * s
  
  par_t <- par
  par_t$lambda <- lambda_t
  par_t$r0     <- r0_t
  par_t$psi_max <- psi_max
  par_t$varsigma0 <- varsigma0_t
  par_t$switch_s <- s
  par_t
}

debt_stress_fun <- function(d, d_bar = 1.05, s0 = 0) {
  ## simple: stress when d/d_bar above 1 (works if d positive); if d can be negative, use abs for stress
  z <- abs(d) / max(1e-9, d_bar) - 1
  pos_part(z - s0)
}

p_growth_fun <- function(R, debt_stress, r, r_star,
                         g_p0 = 0.00, g_pR = 0.015, g_pS = 0.10, g_pRgap = 0.30) {
  ## R is optional regime indicator you may pass; keep it harmless:
  ## asset price growth increases with: (i) regime (financialization), (ii) stress (speculation),
  ## and can be damped if profit rate is below a reference.
  g_p0 + g_pR * R + g_pS * debt_stress + g_pRgap * (r - r_star)
}

## ============================================================
## 8) RHS: full simulation with residual u + log_pA + xi ratchet
## ============================================================
rhs_regime_patriarch <- function(t, state, par) {
  
  par_t <- get_par_t(t, par)
  
  e      <- state[["e"]]
  w      <- state[["omega"]]
  d      <- state[["d"]]
  u      <- state[["u"]]
  log_pA <- state[["log_pA"]]
  xi     <- state[["xi"]]
  xi_floor <- state[["xi_floor"]]
  
  ## guardrails (but do NOT force d sign)
  e <- clamp(e, 1e-9, 1 - 1e-9)
  w <- clamp(w, 1e-9, 1 - 1e-9)
  u <- max(u, 1e-9)
  log_pA <- ifelse(is.finite(log_pA), log_pA, 0)
  xi <- max(xi, 1e-12)
  xi_floor <- max(xi_floor, 1e-12)
  
  ## core evaluated at u_star for e/omega/d block (your “residual u” convention)
  r   <- r_net_reduced(w, d, par_t)
  kap <- kappa_fun(r, par_t)
  
  A   <- par_t$u_star / par_t$sigma
  gK  <- A * kap - par_t$delta
  gY  <- gK
  g_n <- par_t$alpha + par_t$beta
  
  Phi <- Phi_fun(e, par_t)
  de  <- (gY - g_n) * e
  dw  <- (Phi - par_t$alpha) * w
  dd  <- kap - (1 - w) + par_t$i * d - d * gY
  
  ## residual follower u
  u_d <- u_desired_uStar(t, w, d, par_t)
  du <- if (is.finite(u_d)) (par_t$lambda_u %||% 2.0) * (u_d - u) else 0
  
  ## --- asset price growth block (NUMERIC FIXES)
  ## compute stress (raw + clipped)
  debt_stress_raw <- debt_stress_fun(d, d_bar = par_t$d_bar %||% 1.05, s0 = par_t$stress_s0 %||% 0)
  debt_stress <- pmin(debt_stress_raw, par_t$debt_stress_cap %||% 5)
  
  ## regime indicator R can be mapped to the switch itself (or external knob)
  R <- (par_t$R_level_pre %||% 0) * (1 - (par_t$switch_s %||% 0)) +
    (par_t$R_level_post %||% 1) * (par_t$switch_s %||% 0)
  
  p_g_raw <- p_growth_fun(
    R = R,
    debt_stress = debt_stress,
    r = r,
    r_star = par_t$r_star %||% (par_t$r0),
    g_p0 = par_t$g_p0 %||% 0.00,
    g_pR = par_t$g_pR %||% 0.015,
    g_pS = par_t$g_pS %||% 0.10,
    g_pRgap = par_t$g_pRgap %||% 0.30
  )
  
  ## cap p_g (positive + negative caps can differ)
  p_g <- pmin(pmax(p_g_raw, par_t$p_g_cap_neg %||% -0.15), par_t$p_g_cap_pos %||% 0.08)
  
  dlog_pA <- p_g
  
  ## --- inequality block (xi dynamics)
  ## xi grows with positive asset price growth and stress, but has a ratchet to a moving floor
  k_xi_p  <- par_t$k_xi_p  %||% 0.30
  k_xi_d  <- par_t$k_xi_d  %||% 0.20
  k_xi_mr <- par_t$k_xi_mr %||% 0.08
  
  xi_dot <- xi * (
    k_xi_p * pos_part(p_g) +
      k_xi_d * pos_part(debt_stress) -
      k_xi_mr * (xi - xi_floor)
  )
  
  ## xi_floor “ratchet”: rises when regime switches / after clash
  xi_floor_base <- par_t$xi_floor_base %||% 0.60
  xi_floor_jump <- par_t$xi_floor_jump %||% 0.60
  xi_floor_speed <- par_t$xi_floor_speed %||% 0.05
  xi_floor_target <- xi_floor_base + xi_floor_jump * (par_t$switch_s %||% 0)
  
  dxi_floor <- xi_floor_speed * (xi_floor_target - xi_floor)
  
  ## return derivatives + diagnostics
  pA_diag <- exp(pmin(log_pA, 50))
  
  list(
    c(de, dw, dd, du, dlog_pA, xi_dot, dxi_floor),
    c(
      r = r, kappa = kap, gY = gY,
      u_d = u_d, lambda_t = par_t$lambda, r0_t = par_t$r0, psi_max = par_t$psi_max,
      debt_stress = debt_stress_raw, debt_stress_clip = debt_stress,
      R = R, p_g = p_g,
      pA = pA_diag
    )
  )
}

simulate_regime <- function(par, state0 = NULL,
                            t_end = 900, dt = 0.05,
                            perturb = c(e = 0.01, omega = -0.01, d = 0.05, u = -0.05),
                            calibrate_varsigma0 = TRUE) {
  
  ## calibrate varsigma0 to enforce u_d=u_star at SS in pre-regime parameters
  par_use <- par
  if (calibrate_varsigma0) {
    ## use “pre” parameters for calibration if provided
    par_cal <- par_use
    par_cal$lambda <- par_use$lambda_pre %||% par_use$lambda
    par_cal$r0     <- par_use$r0_pre %||% par_use$r0
    par_cal$psi_max <- par_use$psi_max_pre %||% (par_use$psi_max %||% 0)
    par_cal$varsigma0 <- par_use$varsigma0_pre %||% par_use$varsigma0
    
    cal <- calibrate_varsigma0_at_ss(par_cal)
    if (!isTRUE(cal$ok)) stop("Calibration failed: ", cal$reason)
    
    ## store calibrated baseline into the regime container
    par_use$varsigma0_pre <- cal$par$varsigma0
    if (!is.finite(par_use$varsigma0_post %||% NA_real_)) {
      ## default: keep same after switch unless user overrides
      par_use$varsigma0_post <- par_use$varsigma0_pre
    }
  }
  
  ## SS for initial conditions (using pre parameters)
  par_ss <- par_use
  par_ss$lambda <- par_use$lambda_pre %||% par_use$lambda
  par_ss$r0     <- par_use$r0_pre %||% par_use$r0
  par_ss$psi_max <- par_use$psi_max_pre %||% (par_use$psi_max %||% 0)
  par_ss$varsigma0 <- par_use$varsigma0_pre %||% par_use$varsigma0
  
  ss <- steady_state_reduced(par_ss)
  if (!isTRUE(ss$ok)) stop("No admissible reduced SS for initialization: ", ss$reason)
  
  if (is.null(state0)) {
    state0 <- c(
      e        = ss$e + perturb[["e"]],
      omega    = ss$omega + perturb[["omega"]],
      d        = ss$d + perturb[["d"]],
      u        = par_use$u_star + perturb[["u"]],
      log_pA   = log(par_use$pA0 %||% 1.0),
      xi       = par_use$xi0 %||% 1.0,
      xi_floor = par_use$xi_floor0 %||% (par_use$xi_floor_base %||% 0.60)
    )
  }
  
  times <- seq(0, t_end, by = dt)
  
  ## stop conditions (so scenarios don't crash everything)
  rootfun <- function(t, y, parms) {
    c(
      (parms$d_cap_stop %||% 50) - abs(y["d"]),
      (parms$log_pA_cap_stop %||% 40) - y["log_pA"]
    )
  }
  
  sol <- ode(
    y = state0, times = times,
    func = rhs_regime_patriarch, parms = par_use,
    method = "lsoda",
    rootfunc = rootfun
  )
  
  
  ## deSolve returns matrix-like with attributes; convert safely
  df <- as.data.frame(sol)
  
  ## ensure unique names (deSolve can sometimes collide if diagnostics match state names)
  names(df) <- make.unique(names(df))
  
  as_tibble(df) %>%
    mutate(pA = exp(pmin(log_pA, 50)))
}

## ============================================================
## 9) Tail summary + plots
## ============================================================
summarize_tail <- function(df, name, frac_tail = 0.25) {
  dft <- tail_window(df, frac_tail)
  pA_series <- exp(pmin(dft$log_pA, 50))
  
  tibble(
    run = name,
    time_end = max(df$time, na.rm = TRUE),
    
    e_mean = mean(dft$e, na.rm = TRUE),
    e_sd   = sd(dft$e, na.rm = TRUE),
    e_amp  = amp(dft$e),
    
    omega_mean = mean(dft$omega, na.rm = TRUE),
    omega_sd   = sd(dft$omega, na.rm = TRUE),
    omega_amp  = amp(dft$omega),
    
    d_mean = mean(dft$d, na.rm = TRUE),
    d_sd   = sd(dft$d, na.rm = TRUE),
    d_amp  = amp(dft$d),
    
    u_mean = mean(dft$u, na.rm = TRUE),
    u_sd   = sd(dft$u, na.rm = TRUE),
    u_amp  = amp(dft$u),
    
    xi_mean = mean(dft$xi, na.rm = TRUE),
    xi_sd   = sd(dft$xi, na.rm = TRUE),
    xi_amp  = amp(log(pmax(dft$xi, 1e-12))),
    
    gY_mean = mean(dft$gY, na.rm = TRUE),
    
    p_mean = mean(pA_series, na.rm = TRUE),
    p_sd   = sd(pA_series, na.rm = TRUE),
    p_amp  = amp(dft$log_pA),
    
    period_e = period_peaks(dft$time, dft$e)
  )
}

plot_scenario_timeseries <- function(df, tag, out_dir = OUT_DIR) {
  p <- df %>%
    select(time, e, omega, d, u, xi, log_pA, gY, p_g, debt_stress) %>%
    pivot_longer(-time, names_to = "var", values_to = "value") %>%
    ggplot(aes(x = time, y = value)) +
    geom_line() +
    facet_wrap(~var, scales = "free_y", ncol = 1) +
    labs(title = paste0("Scenario time paths (", tag, ")"), x = "t", y = NULL)
  
  ggsave(fs::path(out_dir, paste0("scenario_timeseries_", tag, ".png")),
         p, width = 7.2, height = 10.0, dpi = 160)
  p
}

phase_end_plot <- function(dat, x, y, tag,
                           out_dir = OUT_DIR,
                           tail_frac = 0.25,
                           arrow_every = 40) {
  
  dat <- arrange(dat, time)
  
  n <- nrow(dat)
  i0 <- max(1, floor((1 - tail_frac) * n))
  dat <- dat %>%
    mutate(idx = row_number(),
           is_tail = idx >= i0)
  
  p_start <- dat[1, ]
  p_end   <- dat[n, ]
  
  tail_dat <- dat %>% filter(is_tail) %>%
    mutate(
      x0 = .data[[x]], y0 = .data[[y]],
      x1 = lead(.data[[x]]), y1 = lead(.data[[y]])
    )
  
  arrows <- tail_dat %>% filter(!is.na(x1), idx %% arrow_every == 0)
  
  p <- ggplot(dat, aes(x = .data[[x]], y = .data[[y]])) +
    geom_path(aes(alpha = is_tail)) +
    scale_alpha_manual(values = c(`FALSE` = 0.25, `TRUE` = 1.0), guide = "none") +
    geom_point(data = p_start, size = 2.3) +
    geom_point(data = p_end,   size = 2.3) +
    geom_text(data = p_start, aes(label = "start"), vjust = -0.9, size = 3.2) +
    geom_text(data = p_end,   aes(label = "end"),   vjust = -0.9, size = 3.2) +
    geom_segment(
      data = arrows,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      arrow = arrow(length = unit(0.12, "inches")),
      inherit.aes = FALSE
    ) +
    labs(
      title = paste0("Phase (tail emphasized + direction): ", y, " vs ", x, " (", tag, ")"),
      subtitle = paste0("Tail = last ", round(100*tail_frac), "% | labeled start/end + arrows on tail"),
      x = x, y = y
    )
  
  ggsave(fs::path(out_dir, paste0("phase_end_", y, "_vs_", x, "_", tag, ".png")),
         p, width = 6.6, height = 5.6, dpi = 160)
  p
}

plot_scenario_phase_end <- function(df, tag, out_dir = OUT_DIR) {
  phase_end_plot(df, x = "e", y = "omega", tag = tag, out_dir = out_dir)
  phase_end_plot(df, x = "d", y = "xi",    tag = tag, out_dir = out_dir)
  invisible(TRUE)
}

## ============================================================
## 10) Scenario audit + tidy outputs
## ============================================================
domain_flags_allpath <- function(df) {
  tibble(
    ever_bad_e     = any(!is.finite(df$e) | df$e <= 0 | df$e >= 1),
    ever_bad_omega = any(!is.finite(df$omega) | df$omega <= 0 | df$omega >= 1),
    ever_bad_u     = any(!is.finite(df$u) | df$u <= 0),
    ever_bad_logp  = any(!is.finite(df$log_pA))
  ) %>% mutate(domain_ok_all = !(ever_bad_e | ever_bad_omega | ever_bad_u | ever_bad_logp))
}

classify_scenario <- function(df, frac_tail = 0.25,
                              amp_tol_e = 0.01,
                              amp_tol_omega = 0.005) {
  
  dom <- domain_flags_allpath(df)
  dft <- tail_window(df, frac_tail)
  
  eA <- amp(dft$e)
  wA <- amp(dft$omega)
  
  if (!isTRUE(dom$domain_ok_all)) return("Blow-up / out-of-domain")
  if (is.finite(eA) && is.finite(wA) && (eA >= amp_tol_e || wA >= amp_tol_omega)) {
    return("Cycle / sustained (nontrivial amp)")
  }
  "Convergent / damped (tiny amp)"
}

make_scenarios_tidy <- function(scen_tbl, out_dir = OUT_DIR) {
  scen_tbl2 <- scen_tbl %>%
    mutate(
      lambda_rel = lambda_pre / lambda_star,
      cycle_flag = is.finite(period_e) & period_e > 0 & is.finite(e_amp) & e_amp > 1e-3
    )
  write_csv_quiet(scen_tbl2, fs::path(out_dir, "TABLE_scenarios_tail_diagnostics.csv"))
  
  scen_long <- scen_tbl2 %>%
    select(run, lambda_rel, e_amp, omega_amp, xi_amp, p_amp, period_e, cycle_flag) %>%
    pivot_longer(cols = c(e_amp, omega_amp, xi_amp, p_amp), names_to = "metric", values_to = "value")
  
  p_amp <- ggplot(scen_long, aes(x = lambda_rel, y = value, label = run, shape = cycle_flag)) +
    geom_point(size = 3) +
    geom_text(check_overlap = TRUE, nudge_y = 0.002, size = 3) +
    facet_wrap(~metric, scales = "free_y", ncol = 1) +
    labs(
      title = "Diagnostics: tail amplitudes vs λ/λ*",
      x = "lambda_rel = λ / λ*",
      y = "tail amplitude",
      shape = "cycle_flag"
    )
  ggsave(fs::path(out_dir, "PLOT_scenarios_amplitudes_vs_lambdaRel.png"),
         p_amp, width = 7.6, height = 10.0, dpi = 160)
  
  p_per <- ggplot(scen_tbl2, aes(x = lambda_rel, y = period_e, label = run, shape = cycle_flag)) +
    geom_point(size = 3) +
    geom_text(check_overlap = TRUE, nudge_y = 0.02, size = 3) +
    labs(
      title = "Period proxy: peak-to-peak period of e(t) in tail window",
      x = "lambda_rel = λ / λ*",
      y = "period_e",
      shape = "cycle_flag"
    )
  ggsave(fs::path(out_dir, "PLOT_scenarios_period_vs_lambdaRel.png"),
         p_per, width = 8.2, height = 5.5, dpi = 160)
  
  write_csv_quiet(scen_tbl2, fs::path(out_dir, "TABLE_scenarios_tidy.csv"))
  
  list(scen_tbl = scen_tbl2, p_amp = p_amp, p_per = p_per)
}

## ============================================================
## 11) Outputs manifest
## ============================================================
make_outputs_manifest <- function(out_dir = OUT_DIR) {
  expected <- tribble(
    ~file, ~purpose, ~reviewed,
    "hopf_plane_scores.csv", "Raw plane ranking metrics", FALSE,
    "TABLE_plane_ranking_tidy.csv", "Tidy plane ranking table", FALSE,
    "PLOT_plane_selection_scatter.png", "Plane selection scatter plot", FALSE,
    "hopf_plane_phi1_lambda.csv", "Plane grid data: (phi1, lambda)", FALSE,
    "hopf_plane_r0_lambda.csv", "Plane grid data: (r0, lambda)", FALSE,
    "hopf_plane_phi1_i.csv", "Plane grid data: (phi1, i)", FALSE,
    "FIG_A_hopf_plane_signH_phi1_lambda.png", "Hopf map A", FALSE,
    "FIG_B_hopf_plane_signH_r0_lambda.png", "Hopf map B", FALSE,
    "FIG_C_hopf_plane_signH_phi1_i.png", "Hopf map C", FALSE,
    "TABLE_plane_H_audit.csv", "Plane H distribution audit", FALSE,
    "PLOT_plane_H_quantiles.png", "Plane H quantiles plot", FALSE,
    "TABLE_scenarios_tail_diagnostics.csv", "Tail-window diagnostics per scenario", FALSE,
    "TABLE_scenarios_tidy.csv", "Scenario tidy table", FALSE,
    "PLOT_scenarios_amplitudes_vs_lambdaRel.png", "Amplitudes vs lambda_rel", FALSE,
    "PLOT_scenarios_period_vs_lambdaRel.png", "Period vs lambda_rel", FALSE,
    "scenario_sim_*.csv", "Scenario simulation outputs", FALSE,
    "scenario_timeseries_*.png", "Scenario time-series panels", FALSE,
    "phase_end_*.png", "Phase plots with start/end labels", FALSE
  )
  
  expected %>%
    mutate(
      path = fs::path(out_dir, file),
      exists = if_else(grepl("\\*$", file), TRUE, fs::file_exists(path)),
      status = case_when(
        reviewed ~ "Reviewed",
        exists & !reviewed ~ "Ready to review",
        !exists ~ "Missing",
        TRUE ~ "Unknown"
      )
    ) %>%
    select(file, exists, reviewed, status, purpose)
}

## ============================================================
## 12) DRIVER: baseline parameters + planes + figures + scenarios
## ============================================================
cat("\n[INIT] Setting parameters...\n")

par_base <- list(
  ## growth / tech
  alpha = 0.02, beta = 0.01,
  delta = 0.02,
  sigma = 3.0,
  u_star = 1.0,
  
  ## finance
  i = 0.04,
  
  ## investment (logistic)
  kappa_min = 0.00, kappa_max = 0.30,
  lambda = 20,
  r0 = 0.04,
  
  ## Phillips
  phi0 = -0.06,
  phi1 = 0.10,
  
  ## utilization adjustment speed (residual follower)
  lambda_u = 2.0,
  
  ## price (wealth accounting)
  p = 1.0,
  
  ## autonomous consumption (will be calibrated)
  varsigma0 = 0.55,
  varsigma1 = 0.00,
  
  ## induced non-owner consumption
  cN_minus = 0.02,
  cN_plus  = 0.20,
  cN_eta   = 0.10,
  
  ## owners saving
  sA_const  = TRUE,
  sA0       = 0.40,
  sA_min    = 0.20, sA_max = 0.80, sA_lambda = 2.0, sA_mid = 2.0,
  
  ## leakage psi(d)
  psi_zero  = TRUE,
  psi_max   = 0.00, psi_lambda = 5.0, psi_mid = 1.0,
  
  ## --- stability / caps (THE FIXES)
  debt_stress_cap <- 5,
  p_g_cap_pos     <- 0.08,
  p_g_cap_neg     <- -0.15,
  d_cap_stop      <- 50,
  log_pA_cap_stop <- 40,
  
  ## --- wealth block defaults
  pA0 = 1.0,
  xi0 = 1.0,
  xi_floor_base = 0.60,
  xi_floor0 = 0.60,
  xi_floor_jump = 0.60,
  xi_floor_speed = 0.05,
  
  k_xi_p  = 0.30,
  k_xi_d  = 0.20,
  k_xi_mr = 0.08,
  
  ## asset price growth knobs
  g_p0 = 0.00,
  g_pR = 0.015,
  g_pS = 0.10,
  g_pRgap = 0.30,
  r_star = 0.04,
  
  ## regime switch schedule
  t_switch = 350,
  t_width  = 40,
  R_level_pre  = 0,
  R_level_post = 1,
  d_bar = 1.05,
  stress_s0 = 0
)

planes <- list(
  c("phi1","lambda"),
  c("i","lambda"),
  c("phi1","i"),
  c("phi0","phi1"),
  c("r0","lambda")
)

ranges <- list(
  phi1   = c(0.03, 0.25),
  lambda = c(5, 150),
  i      = c(0.00, 0.10),
  phi0   = c(-0.12, -0.01),
  r0     = c(0.01, 0.08)
)

## -----------------------------
## A) Rank planes + decision plot
## -----------------------------
cat("\n[A] Ranking planes...\n")
scores <- compare_planes(
  par_base,
  planes = planes,
  ranges = ranges,
  n = 80,
  H_band = 1e-3,
  tol = 1e-6,
  calibrate_varsigma0 = TRUE
)
write_csv_quiet(scores, fs::path(OUT_DIR, "hopf_plane_scores.csv"))
plane_tbl <- make_plane_decision_outputs(scores, out_dir = OUT_DIR)
cat("[A] Saved hopf_plane_scores.csv + TABLE_plane_ranking_tidy.csv + scatter plot\n")

## -----------------------------
## B) Hopf maps for selected planes
## -----------------------------
cat("\n[B] Building Hopf maps...\n")

## FIG A: (phi1, lambda)
df_A <- hopf_plane(
  par_base, "phi1", "lambda",
  p1_vals = grid_values(ranges$phi1[1], ranges$phi1[2], 140),
  p2_vals = grid_values(ranges$lambda[1], ranges$lambda[2], 140),
  calibrate_varsigma0 = TRUE,
  H_band = 1e-3
)
write_csv_quiet(df_A, fs::path(OUT_DIR, "hopf_plane_phi1_lambda.csv"))
pA <- plot_hopf_plane_signH(df_A, "phi1","lambda","FIG A: Hopf map (lambda vs phi1)")
ggsave(fs::path(OUT_DIR, "FIG_A_hopf_plane_signH_phi1_lambda.png"), pA, width = 7.2, height = 6.0, dpi = 160)

## FIG B: (r0, lambda)
df_B <- hopf_plane(
  par_base, "r0", "lambda",
  p1_vals = grid_values(ranges$r0[1], ranges$r0[2], 140),
  p2_vals = grid_values(ranges$lambda[1], ranges$lambda[2], 140),
  calibrate_varsigma0 = TRUE,
  H_band = 1e-3
)
write_csv_quiet(df_B, fs::path(OUT_DIR, "hopf_plane_r0_lambda.csv"))
pB <- plot_hopf_plane_signH(df_B, "r0","lambda","FIG B: Hopf map (lambda vs r0)")
ggsave(fs::path(OUT_DIR, "FIG_B_hopf_plane_signH_r0_lambda.png"), pB, width = 7.2, height = 6.0, dpi = 160)

## FIG C: (phi1, i)
df_C <- hopf_plane(
  par_base, "phi1", "i",
  p1_vals = grid_values(ranges$phi1[1], ranges$phi1[2], 140),
  p2_vals = grid_values(ranges$i[1], ranges$i[2], 140),
  calibrate_varsigma0 = TRUE,
  H_band = 1e-3
)
write_csv_quiet(df_C, fs::path(OUT_DIR, "hopf_plane_phi1_i.csv"))
pC <- plot_hopf_plane_signH(df_C, "phi1","i","FIG C: Hopf map (i vs phi1)")
ggsave(fs::path(OUT_DIR, "FIG_C_hopf_plane_signH_phi1_i.png"), pC, width = 7.2, height = 6.0, dpi = 160)

cat("[B] Saved Hopf maps A/B/C and plane CSVs\n")

## -----------------------------
## C) Plane audit bundle
## -----------------------------
cat("\n[C] Auditing H distributions...\n")
plane_audits <- bind_rows(
  audit_plane_df(df_A, "hopf_plane_phi1_lambda"),
  audit_plane_df(df_B, "hopf_plane_r0_lambda"),
  audit_plane_df(df_C, "hopf_plane_phi1_i")
)
write_csv_quiet(plane_audits, fs::path(OUT_DIR, "TABLE_plane_H_audit.csv"))
plot_plane_H_quantiles(plane_audits, out_dir = OUT_DIR)
cat("[C] Saved TABLE_plane_H_audit.csv + PLOT_plane_H_quantiles.png\n")

## -----------------------------
## D) Find lambda* on the (phi1, lambda) slice at phi1_pick
## -----------------------------
cat("\n[D] Finding lambda*...\n")
phi1_pick <- par_base$phi1
par_root <- par_base
par_root$phi1 <- phi1_pick

lambda_star <- find_hopf_mu(
  par_root, mu = "lambda",
  lo = ranges$lambda[1], hi = ranges$lambda[2],
  n_grid = 350,
  calibrate_varsigma0 = TRUE,
  verbose = TRUE
)

cat(sprintf("\n[D] Hopf approx lambda* = %.6f at phi1 = %.4f\n", lambda_star, phi1_pick))

## -----------------------------
## E) Scenarios: fixed lambda below/near/above + moving-Hopffiness flagship
## -----------------------------
cat("\n[E] Running scenarios...\n")
scenario_tbl <- list()

## E1: fixed lambda
lambda_vals <- c(below = 0.70*lambda_star, near = 0.98*lambda_star, above = 1.20*lambda_star)
lambda_vals <- pmin(pmax(lambda_vals, ranges$lambda[1]), ranges$lambda[2])

for (nm in names(lambda_vals)) {
  
  par_s <- par_base
  par_s$phi1 <- phi1_pick
  
  ## fixed regime: pre=post
  par_s$lambda_pre  <- as.numeric(lambda_vals[[nm]])
  par_s$lambda_post <- as.numeric(lambda_vals[[nm]])
  
  par_s$r0_pre  <- par_base$r0
  par_s$r0_post <- par_base$r0
  
  par_s$psi_zero <- TRUE
  par_s$psi_max_pre  <- 0.00
  par_s$psi_max_post <- 0.00
  
  tag <- paste0("fixed_phi1_", sprintf("%.3f", phi1_pick),
                "_lambda_", nm, "_", sprintf("%.2f", par_s$lambda_pre))
  
  df_s <- simulate_regime(par_s, t_end = 900, dt = 0.05)
  
  write_csv_quiet(df_s, fs::path(OUT_DIR, paste0("scenario_sim_", tag, ".csv")))
  plot_scenario_timeseries(df_s, tag, out_dir = OUT_DIR)
  plot_scenario_phase_end(df_s, tag, out_dir = OUT_DIR)
  
  scenario_tbl[[tag]] <- summarize_tail(df_s, name = tag, frac_tail = 0.25) %>%
    mutate(
      phi1 = phi1_pick,
      lambda_pre = par_s$lambda_pre,
      lambda_post = par_s$lambda_post,
      lambda_star = lambda_star,
      scenario_type = "fixed_lambda",
      sim_file = fs::path(OUT_DIR, paste0("scenario_sim_", tag, ".csv")),
      sim_exists = TRUE
    )
}

## E2: flagship moving-Hopffiness (patriarch)
## pick a pre below lambda* and a post above lambda*
par_m <- par_base
par_m$phi1 <- phi1_pick
par_m$lambda_pre  <- pmin(pmax(0.85*lambda_star, ranges$lambda[1]), ranges$lambda[2])
par_m$lambda_post <- pmin(pmax(1.20*lambda_star, ranges$lambda[1]), ranges$lambda[2])

## drift r0 upward post-switch (tightens investment logistic threshold)
par_m$r0_pre  <- par_base$r0
par_m$r0_post <- 0.065

## turn on leakage after switch
par_m$psi_zero <- FALSE
par_m$psi_max_pre  <- 0.00
par_m$psi_max_post <- 0.08
par_m$psi_lambda <- 5.0
par_m$psi_mid    <- 1.0

## stronger inequality ratchet after switch
par_m$xi_floor_jump <- 0.90
par_m$k_xi_p  <- 0.35
par_m$k_xi_d  <- 0.25
par_m$k_xi_mr <- 0.06

tag_m <- paste0("movingHopf_phi1_", sprintf("%.3f", phi1_pick),
                "_lamPre_", sprintf("%.2f", par_m$lambda_pre),
                "_lamPost_", sprintf("%.2f", par_m$lambda_post),
                "_r0Drift")

df_m <- simulate_regime(par_m, t_end = 900, dt = 0.05)
write_csv_quiet(df_m, fs::path(OUT_DIR, paste0("scenario_sim_", tag_m, ".csv")))
plot_scenario_timeseries(df_m, tag_m, out_dir = OUT_DIR)
plot_scenario_phase_end(df_m, tag_m, out_dir = OUT_DIR)

scenario_tbl[[tag_m]] <- summarize_tail(df_m, name = tag_m, frac_tail = 0.25) %>%
  mutate(
    phi1 = phi1_pick,
    lambda_pre = par_m$lambda_pre,
    lambda_post = par_m$lambda_post,
    lambda_star = lambda_star,
    scenario_type = "movingHopf_patriarch",
    sim_file = fs::path(OUT_DIR, paste0("scenario_sim_", tag_m, ".csv")),
    sim_exists = TRUE
  )

scenario_tbl <- bind_rows(scenario_tbl)

## classify scenarios using whole-path checks + tail amplitude threshold
scen_audit <- scenario_tbl %>%
  rowwise() %>%
  mutate(
    sim_exists = fs::file_exists(sim_file),
    regime2 = if (sim_exists) {
      df <- readr::read_csv(sim_file, show_col_types = FALSE)
      classify_scenario(df)
    } else NA_character_
  ) %>%
  ungroup()

write_csv_quiet(scen_audit, fs::path(OUT_DIR, "TABLE_scenarios_audit_reclassified.csv"))

## tidy scenario outputs + plots
scen_outputs <- make_scenarios_tidy(scen_audit %>% mutate(lambda_star = lambda_star), out_dir = OUT_DIR)

cat("[E] Saved scenarios: sims + timeseries + phase plots + tables + amp/period plots\n")

## -----------------------------
## F) Manifest
## -----------------------------
manifest <- make_outputs_manifest(OUT_DIR)
write_csv_quiet(manifest, fs::path(OUT_DIR, "TABLE_outputs_review_manifest.csv"))
cat("[F] Saved outputs manifest\n")

cat("\nDONE. Exports in:\n", OUT_DIR, "\n")
############################################################
