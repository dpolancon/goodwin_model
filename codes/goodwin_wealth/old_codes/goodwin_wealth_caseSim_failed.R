############################################################
## wealth_goodwin_wealth_inequality_RUNALL.R
##
## Purpose:
##  - Reduced-system Hopf mapping (SS -> Jacobian -> RH -> H) on parameter planes
##  - Rank planes, export Hopf sign(H) maps + H audits
##  - Find lambda* (Hopf root) for a chosen phi1
##  - Full simulation with residual u follower AND regime patriarch (moving Hopf)
##  - Wealth inequality block: asset-price process pA + xi dynamics with ratcheting floor
##  - Export: scenario sims, timeseries panels, phase plots with start/end + arrows,
##            tail diagnostics tables, scenario audit tables, and a manifest
############################################################

suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(fs)
  library(tibble)
  library(grid)   # unit() in arrows
  library(readr)
  library(stringr)
})

## ============================================================
## OUTPUT DIR
## ============================================================
OUT_DIR <- fs::path("outputs", "wealth_goodwin", "wealth_inequality")
fs::dir_create(OUT_DIR)

write_csv_quiet <- function(df, path) readr::write_csv(df, path)

## ============================================================
## 0) Helpers
## ============================================================
inv_logit <- function(x) 1 / (1 + exp(-x))
pos_part  <- function(x) pmax(0, x)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || !is.finite(a[1])) b else a

require_params <- function(par, keys) {
  miss <- setdiff(keys, names(par))
  if (length(miss) > 0) stop("Missing parameters: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

grid_values <- function(lo, hi, n) seq(lo, hi, length.out = n)

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
## 1) Core reduced model pieces (for Hopf mapping)
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

## owners saving (optional)
sA_fun <- function(chiA, par) {
  if (isTRUE(par$sA_const)) return(par$sA0)
  s <- inv_logit(par$sA_lambda * (chiA - par$sA_mid))
  par$sA_min + (par$sA_max - par$sA_min) * s
}

## leakage psi(d)
psi_fun <- function(d, par) {
  if (isTRUE(par$psi_zero)) return(0)
  inv_logit(par$psi_lambda * (d - par$psi_mid)) * par$psi_max
}

## induced non-owner consumption "share"
cN_ind_fun <- function(yN, par) {
  if (yN <= 0) return(par$cN_minus)
  par$cN_minus + (par$cN_plus - par$cN_minus) * (yN / (yN + par$cN_eta))
}

## autonomous baseline consumption share
cN_bar_fun <- function(t, par) par$varsigma0 * exp(par$varsigma1 * t)

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
## 2) Reduced steady state (u=u_star), and calibrate varsigma0 at SS
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
  if (A <= 0) return(list(ok = FALSE, reason = "A=u_star/sigma must be >0"))
  
  k_target <- (g_n + par$delta) / A
  if (k_target <= par$kappa_min + eps || k_target >= par$kappa_max - eps) {
    return(list(ok = FALSE, reason = "kappa_target outside (kappa_min,kappa_max)"))
  }
  
  s <- (k_target - par$kappa_min) / (par$kappa_max - par$kappa_min)
  if (s <= eps || s >= 1 - eps) return(list(ok = FALSE, reason = "logit inversion s not in (0,1)"))
  
  r_star <- par$r0 + (1 / par$lambda) * log(s / (1 - s))
  
  d_star <- (k_target - (1 / A) * r_star) / g_n
  if (d_star <= eps) return(list(ok = FALSE, reason = "d* <= 0 (no interior debt steady state)"))
  
  pi_star <- k_target + (par$i - g_n) * d_star
  omega_star <- 1 - pi_star
  if (omega_star <= eps || omega_star >= 1 - eps) return(list(ok = FALSE, reason = "omega* not in (0,1)"))
  
  if (par$phi1 <= 0) return(list(ok = FALSE, reason = "phi1 must be > 0"))
  e_star <- (par$alpha - par$phi0) / par$phi1
  if (e_star <= eps || e_star >= 1 - eps) return(list(ok = FALSE, reason = "e* not in (0,1)"))
  
  list(ok = TRUE, e = e_star, omega = omega_star, d = d_star,
       r = r_star, kappa = k_target, g_n = g_n, A = A)
}

calibrate_varsigma0_at_ss <- function(par, t0 = 0) {
  require_params(par, c("u_star","varsigma1","varsigma0"))
  if (abs(par$varsigma1) > 0) warning("varsigma1 != 0: calibration enforces u_d=u_star only at t0.")
  
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
## 3) Reduced Jacobian + RH/Hopf functional
## ============================================================
jacobian_reduced_at_ss <- function(ss, par) {
  e <- ss$e; omega <- ss$omega; d <- ss$d
  g_n <- ss$g_n
  A   <- ss$A
  r   <- ss$r
  
  kap_r <- kappa_r_fun(r, par)
  Phi_e <- Phi_e_fun(e, par)
  
  J <- matrix(0, 3, 3)
  
  # de
  J[1,2] <- - e * (A^2) * kap_r
  J[1,3] <- - e * (A^2) * par$i * kap_r
  
  # domega
  J[2,1] <- omega * Phi_e
  
  # dd
  J[3,2] <- 1 - A * kap_r + d * (A^2) * kap_r
  J[3,3] <- (par$i - g_n) - A * par$i * kap_r + d * (A^2) * par$i * kap_r
  
  list(J = J, kappa_r = kap_r, Phi_e = Phi_e)
}

rh_hopf_reduced <- function(J) {
  trJ <- sum(diag(J))
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
  eig <- eigen(Jp$J)$values
  
  list(ok = TRUE, par = par, steady_state = ss, jacobian = Jp$J,
       eigenvalues = eig, rh = rh,
       jac_parts = list(kappa_r = Jp$kappa_r, Phi_e = Jp$Phi_e))
}

hopf_gap_reduced <- function(par, calibrate_varsigma0 = TRUE) {
  out <- analyze_reduced(par, calibrate_varsigma0 = calibrate_varsigma0)
  if (!isTRUE(out$ok)) return(NA_real_)
  out$rh$H
}

find_hopf_mu <- function(par, mu, lo, hi,
                         n_grid = 300, tol = 1e-8, maxit = 200,
                         calibrate_varsigma0 = TRUE,
                         verbose = TRUE) {
  
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
## 4) Hopf plane grid + scoring + plots + audits
## ============================================================
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
      share_H_pos = NA_real_,
      share_H_neg = NA_real_,
      H_min = NA_real_, H_max = NA_real_,
      H_med = NA_real_, H_q05 = NA_real_, H_q95 = NA_real_,
      sd_H = NA_real_,
      sign_change = FALSE,
      has_robust_sign_change = FALSE,
      hopf_band_share = NA_real_,
      degeneracy_flag = TRUE
    ))
  }
  
  Hmin <- min(dfF$H); Hmax <- max(dfF$H)
  share_pos <- mean(dfF$H >  tol)
  share_neg <- mean(dfF$H < -tol)
  
  sdH <- sd(dfF$H)
  
  has_robust_sign_change <- (Hmin < -tol) && (Hmax > tol)
  hopf_band_share <- mean(abs(dfF$H) < H_band)
  
  tibble(
    n_ok = sum(df$ok, na.rm = TRUE),
    n_feasible = nrow(dfF),
    feasible_share = nrow(dfF) / max(1, sum(df$ok, na.rm = TRUE)),
    hopf_band_share = hopf_band_share,
    share_H_pos = share_pos,
    share_H_neg = share_neg,
    H_min = Hmin,
    H_max = Hmax,
    H_med = median(dfF$H),
    H_q05 = as.numeric(quantile(dfF$H, 0.05)),
    H_q95 = as.numeric(quantile(dfF$H, 0.95)),
    sd_H = sdH,
    sign_change = (Hmin < 0 && Hmax > 0),
    has_robust_sign_change = has_robust_sign_change,
    degeneracy_flag = (!has_robust_sign_change) || (sdH < 1e-8)
  )
}

plot_hopf_plane_signH <- function(df, p1, p2, title = NULL) {
  
  df2 <- df %>%
    mutate(
      status = case_when(
        !ok ~ "SS infeasible",
        ok & !feasible_RH ~ "RH prereq fail",
        ok & feasible_RH & is.finite(H) & H > 0 ~ "H>0",
        ok & feasible_RH & is.finite(H) & H < 0 ~ "H<0",
        ok & feasible_RH & is.finite(H) & abs(H) <= 1e-10 ~ "H≈0",
        TRUE ~ "Other"
      ),
      status = factor(status, levels = c("SS infeasible","RH prereq fail","H>0","H<0","H≈0","Other"))
    )
  
  ggplot(df2, aes(x = .data[[p1]], y = .data[[p2]])) +
    geom_tile(aes(fill = status), alpha = 0.92) +
    geom_contour(aes(z = H), breaks = 0, linewidth = 0.55, na.rm = TRUE) +
    labs(
      title = title %||% paste0("Hopf plane: ", p2, " vs ", p1),
      subtitle = "Tiles: SS feasibility + RH prereqs + sign(H). Contour: H=0 (candidate Hopf boundary).",
      x = p1, y = p2, fill = NULL
    )
}

plane_H_audit <- function(df, plane_name) {
  dfF <- df %>% filter(ok, feasible_RH, is.finite(H))
  if (nrow(dfF) == 0) {
    return(tibble(
      plane = plane_name, n_ok = sum(df$ok, na.rm = TRUE), n_feasible = 0,
      H_min = NA_real_, H_q05 = NA_real_, H_med = NA_real_, H_q95 = NA_real_, H_max = NA_real_,
      share_absH_lt_1e_06 = NA_real_, share_absH_lt_1e_04 = NA_real_,
      share_absH_lt_1e_03 = NA_real_, share_absH_lt_1e_02 = NA_real_
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

make_plane_ranking_outputs <- function(scores, out_dir = OUT_DIR) {
  plane_tbl <- scores %>%
    mutate(
      plane = paste0("(", p1, ", ", p2, ")"),
      score = 100*feasible_share + 25*hopf_band_share + 10*as.numeric(has_robust_sign_change) - 5*as.numeric(degeneracy_flag)
    ) %>%
    arrange(desc(has_robust_sign_change), degeneracy_flag, desc(feasible_share), desc(sd_H)) %>%
    select(plane, p1, p2,
           has_robust_sign_change, degeneracy_flag, sign_change,
           feasible_share, hopf_band_share, sd_H, H_min, H_max,
           n_feasible, n_ok, score)
  
  write_csv_quiet(plane_tbl, fs::path(out_dir, "TABLE_plane_ranking_tidy.csv"))
  
  p_sel <- plane_tbl %>%
    ggplot(aes(x = feasible_share, y = hopf_band_share, label = plane,
               shape = as.factor(has_robust_sign_change))) +
    geom_point(size = 2.8) +
    geom_text(check_overlap = TRUE, nudge_y = 0.002, size = 3) +
    labs(
      title = "Plane selection: feasibility vs near-Hopf thickness",
      subtitle = "Prefer robust sign-change + high feasible_share. hopf_band_share = density of |H|<band within feasible RH.",
      x = "feasible_share (RH prereqs satisfied among ok SS points)",
      y = "hopf_band_share (|H| < band, within feasible RH)",
      shape = "robust sign-change?"
    )
  
  ggsave(fs::path(out_dir, "PLOT_plane_selection_scatter.png"),
         p_sel, width = 9.2, height = 6.0, dpi = 160)
  
  list(plane_tbl = plane_tbl, plot = p_sel)
}

## ============================================================
## 5) Patriarch regime module: moving Hopf + asset prices + xi
## ============================================================
R_time <- function(t, tc, k) inv_logit(k * (t - tc))
R_mu <- function(mu_t, mu_star, k) inv_logit(k * (mu_t - mu_star))
interp_regime <- function(R, pre, post) pre + (post - pre) * R

debt_stress_fun <- function(d, d_bar, s0 = 0) pos_part(d - d_bar - s0)

p_growth_fun <- function(R, debt_stress, r, r_star,
                         g_p0 = 0.00, g_pR = 0.02,
                         g_pS = 0.10, g_pRgap = 0.30) {
  g_p0 + g_pR * R + g_pS * debt_stress - g_pRgap * pos_part(r - r_star)
}

xi_dot_fun <- function(xi, xi_floor, p_g, debt_stress, R,
                       k_xi_p_pre = 0.00, k_xi_p_post = 0.20,
                       k_xi_d_pre = 0.02, k_xi_d_post = 0.10,
                       k_xi_mr_pre = 0.06, k_xi_mr_post = 0.01) {
  
  k_xi_p  <- interp_regime(R, k_xi_p_pre,  k_xi_p_post)
  k_xi_d  <- interp_regime(R, k_xi_d_pre,  k_xi_d_post)
  k_xi_mr <- interp_regime(R, k_xi_mr_pre, k_xi_mr_post)
  
  xi * (k_xi_p * pos_part(p_g) + k_xi_d * pos_part(debt_stress)) -
    k_xi_mr * pos_part(xi - xi_floor)
}

xi_floor_dot_fun <- function(xi_floor, R, debt_stress,
                             eta_rat = 0.08,
                             s0 = 0.00,
                             delta_floor = 0.005,
                             xi_floor_base = 0.50) {
  ratchet <- eta_rat * R * pos_part(debt_stress - s0)
  relax   <- delta_floor * (xi_floor - xi_floor_base)
  ratchet - relax
}

## ============================================================
## 6) Full simulation RHS: residual u follower + regime patriarch
##     States: e, omega, d, u, pA, xi, xi_floor
## ============================================================
rhs_regime_patriarch <- function(t, state, par) {
  
  e        <- state[["e"]]
  omega    <- state[["omega"]]
  d        <- state[["d"]]
  u        <- state[["u"]]
  pA       <- state[["pA"]]
  xi       <- state[["xi"]]
  xi_floor <- state[["xi_floor"]]
  
  # guardrails
  e <- min(max(e, 1e-9), 1 - 1e-9)
  omega_floor <- par$omega_floor %||% 0.05
  omega <- min(max(omega, omega_floor), 1 - 1e-9)
  if (!is.finite(d)) d <- 0
  u  <- max(u, 1e-9)
  pA <- max(pA, 1e-9)
  xi <- max(xi, 1e-9)
  xi_floor <- max(xi_floor, 1e-9)
  
  # patriarch drift: lambda(t) moves across lambda_star
  R_t  <- R_time(t, par$t_crash %||% 350, par$k_R %||% 0.03)
  lambda_raw <- interp_regime(R_t, par$lambda_pre, par$lambda_post)
  R <- R_mu(lambda_raw, par$lambda_star, par$k_mu %||% 80)
  
  # time-varying params (subordinates)
  lambda_t <- interp_regime(R, par$lambda_pre,  par$lambda_post)
  r0_t     <- interp_regime(R, par$r0_pre,      par$r0_post)
  psi_max  <- interp_regime(R, par$psi_max_pre, par$psi_max_post)
  
  par_t <- par
  par_t$lambda <- lambda_t
  par_t$r0     <- r0_t
  par_t$psi_max <- psi_max
  par_t$psi_zero <- FALSE  # allow psi_max to matter
  
  # core dynamics (still uses u_star in r and gY)
  r   <- r_net_reduced(omega, d, par_t)
  kap <- kappa_fun(r, par_t)
  
  A   <- (par_t$u_star / par_t$sigma)
  g_n <- par_t$alpha + par_t$beta
  
  gY <- A * kap - par_t$delta  # NOTE: not forced to 0
  de <- (gY - g_n) * e
  dw <- (Phi_fun(e, par_t) - par_t$alpha) * omega
  dd <- kap - (1 - omega) + par_t$i * d - d * gY  # d can go +/- in your interpretation
  
  # residual follower u
  u_d <- u_desired_uStar(t, omega, d, par_t)
  du  <- if (is.finite(u_d)) par_t$lambda_u * (u_d - u) else 0
  
  # stress + asset price drift
  debt_stress <- debt_stress_fun(d, d_bar = par_t$d_bar %||% 1.05, s0 = par_t$stress_s0 %||% 0)
  p_g <- p_growth_fun(R, debt_stress, r, r_star = par_t$r_star %||% r0_t,
                      g_p0 = par_t$g_p0 %||% 0.00,
                      g_pR = par_t$g_pR %||% 0.015,
                      g_pS = par_t$g_pS %||% 0.10,
                      g_pRgap = par_t$g_pRgap %||% 0.30)
  
  dpA <- pA * p_g
  
  # inequality + ratchet floor
  dxi <- xi_dot_fun(xi, xi_floor, p_g, debt_stress, R,
                    k_xi_p_pre  = par_t$k_xi_p_pre  %||% 0.00,
                    k_xi_p_post = par_t$k_xi_p_post %||% 0.20,
                    k_xi_d_pre  = par_t$k_xi_d_pre  %||% 0.02,
                    k_xi_d_post = par_t$k_xi_d_post %||% 0.10,
                    k_xi_mr_pre = par_t$k_xi_mr_pre %||% 0.06,
                    k_xi_mr_post= par_t$k_xi_mr_post%||% 0.01)
  
  dxi_floor <- xi_floor_dot_fun(xi_floor, R, debt_stress,
                                eta_rat = par_t$eta_rat %||% 0.08,
                                s0 = par_t$rat_s0 %||% 0.00,
                                delta_floor = par_t$delta_floor %||% 0.005,
                                xi_floor_base = par_t$xi_floor_base %||% 0.60)
  
  list(
    c(de, dw, dd, du, dpA, dxi, dxi_floor),
    c(R = R, lambda_t = lambda_t, r0_t = r0_t, psi_max = psi_max,
      r = r, kappa = kap, gY = gY,
      debt_stress = debt_stress, p_g = p_g,
      u_d = u_d, u = u)
  )
}

rhs_regime_guarded <- function(t, y, parms) {
  # state sanity
  if (any(!is.finite(y))) {
    msg <- paste0(
      "Non-finite STATE at t=", signif(t, 6), "\n",
      paste(names(y), signif(y, 6), sep="=", collapse=", ")
    )
    stop(msg, call. = FALSE)
  }
  
  out <- rhs_regime_patriarch(t, y, parms)
  
  # deSolve convention: out[[1]] is dy vector
  dy <- out[[1]]
  if (any(!is.finite(dy))) {
    msg <- paste0(
      "Non-finite DERIVATIVE at t=", signif(t, 6), "\n",
      paste(names(y), signif(y, 6), sep="=", collapse=", ")
    )
    stop(msg, call. = FALSE)
  }
  
  out
}


simulate_regime <- function(par, state0 = NULL,
                            t_end = 900, dt = 0.05,
                            perturb = c(e = 0.01, omega = -0.01, d = 0.05, u = -0.05),
                            calibrate_varsigma0 = TRUE) {
  
  # calibration at SS uses PRE regime parameters (R=0) by design
  par_pre <- par
  par_pre$lambda <- par$lambda_pre
  par_pre$r0     <- par$r0_pre
  par_pre$psi_max <- par$psi_max_pre
  par_pre$psi_zero <- TRUE
  
  if (calibrate_varsigma0) {
    cal <- calibrate_varsigma0_at_ss(par_pre)
    if (!isTRUE(cal$ok)) stop("Calibration failed: ", cal$reason)
    par_pre <- cal$par
    ss <- cal$ss
  } else {
    ss <- steady_state_reduced(par_pre)
    if (!isTRUE(ss$ok)) stop("No admissible SS: ", ss$reason)
  }
  
  # carry calibrated varsigma0 into par
  par$varsigma0 <- par_pre$varsigma0
  
  if (is.null(state0)) {
    state0 <- c(
      e        = ss$e + perturb[["e"]],
      omega    = ss$omega + perturb[["omega"]],
      d        = ss$d + perturb[["d"]],
      u        = par$u_star + perturb[["u"]],
      pA       = par$pA0 %||% 1.0,
      xi       = par$xi0 %||% 1.0,
      xi_floor = par$xi_floor0 %||% (par$xi_floor_base %||% 0.60)
    )
  }
  
  times <- seq(0, t_end, by = dt)
  
  sol <- ode(
    y = state0,
    times = times,
    func = rhs_regime_guarded,
    parms = par,
    method = "lsoda"
  )
  
  as_tibble(as.data.frame(sol))
}

## ============================================================
## 7) Diagnostics + plots (scenarios)
## ============================================================
summarize_tail <- function(df, name, frac_tail = 0.25) {
  dft <- tail_window(df, frac_tail)
  
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
    
    pA_mean = mean(dft$pA, na.rm = TRUE),
    pA_sd   = sd(dft$pA, na.rm = TRUE),
    pA_amp  = amp(log(pmax(dft$pA, 1e-12))),
    
    xi_mean = mean(dft$xi, na.rm = TRUE),
    xi_sd   = sd(dft$xi, na.rm = TRUE),
    xi_amp  = amp(log(pmax(dft$xi, 1e-12))),
    
    gY_mean = mean(dft$gY, na.rm = TRUE),
    
    period_e = period_peaks(dft$time, dft$e)
  )
}

domain_flags <- function(df) {
  tibble(
    ever_bad_e = any(!is.finite(df$e) | df$e <= 0 | df$e >= 1),
    ever_bad_omega = any(!is.finite(df$omega) | df$omega <= 0 | df$omega >= 1),
    ever_bad_u = any(!is.finite(df$u) | df$u <= 0),
    ever_bad_pA = any(!is.finite(df$pA) | df$pA <= 0),
    ever_bad_xi = any(!is.finite(df$xi) | df$xi <= 0),
    domain_ok_all = TRUE
  ) %>%
    mutate(domain_ok_all = !(ever_bad_e | ever_bad_omega | ever_bad_u | ever_bad_pA | ever_bad_xi))
}

classify_scenario <- function(df, frac_tail = 0.25,
                              amp_tol_e = 0.01, amp_tol_omega = 0.005) {
  
  dom <- domain_flags(df)
  dft <- tail_window(df, frac_tail)
  
  eA <- amp(dft$e)
  wA <- amp(dft$omega)
  
  if (!dom$domain_ok_all) return("Blow-up / out-of-domain")
  if (is.finite(eA) && is.finite(wA) && (eA >= amp_tol_e || wA >= amp_tol_omega)) {
    return("Cycle / sustained (nontrivial amp)")
  }
  "Convergent / damped (tiny amp)"
}

plot_scenario_timeseries <- function(df, tag, out_dir = OUT_DIR) {
  vars <- c("e","omega","d","u","pA","xi","xi_floor","gY","R","lambda_t","r0_t","debt_stress","p_g")
  keep <- intersect(vars, names(df))
  
  p <- df %>%
    select(time, all_of(keep)) %>%
    pivot_longer(-time, names_to = "var", values_to = "value") %>%
    ggplot(aes(x = time, y = value)) +
    geom_line() +
    facet_wrap(~var, scales = "free_y", ncol = 1) +
    labs(title = paste0("Scenario time paths (", tag, ")"), x = "t", y = NULL)
  
  ggsave(fs::path(out_dir, paste0("scenario_timeseries_", tag, ".png")),
         p, width = 7.6, height = 12.0, dpi = 160)
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
         p, width = 6.8, height = 5.8, dpi = 160)
  p
}

plot_scenario_phase_end <- function(df, tag, out_dir = OUT_DIR) {
  phase_end_plot(df, x = "e",  y = "omega", tag = tag, out_dir = out_dir)
  phase_end_plot(df, x = "d",  y = "xi",    tag = tag, out_dir = out_dir)
  phase_end_plot(df, x = "pA", y = "xi",    tag = tag, out_dir = out_dir)
  invisible(TRUE)
}

## ============================================================
## 8) Audit manifest
## ============================================================
make_outputs_manifest <- function(out_dir = OUT_DIR) {
  expected <- tribble(
    ~file, ~purpose, ~reviewed,
    "hopf_plane_scores.csv", "Raw plane scoring metrics", FALSE,
    "TABLE_plane_ranking_tidy.csv", "Ranked plane table", FALSE,
    "PLOT_plane_selection_scatter.png", "Plane selection plot", FALSE,
    
    "hopf_plane_phi1_lambda.csv", "Plane grid data (phi1,lambda)", FALSE,
    "FIG_A_hopf_plane_signH_phi1_lambda.png", "Hopf map A", FALSE,
    
    "hopf_plane_r0_lambda.csv", "Plane grid data (r0,lambda)", FALSE,
    "FIG_B_hopf_plane_signH_r0_lambda.png", "Hopf map B", FALSE,
    
    "hopf_plane_phi1_i.csv", "Plane grid data (phi1,i)", FALSE,
    "FIG_C_hopf_plane_signH_phi1_i.png", "Hopf map C", FALSE,
    
    "TABLE_plane_H_audit.csv", "H distribution audit across planes", FALSE,
    "PLOT_plane_H_quantiles.png", "H quantile intervals across planes", FALSE,
    
    "TABLE_scenarios_tail_diagnostics.csv", "Tail diagnostics per scenario", FALSE,
    "TABLE_scenarios_audit_reclassified.csv", "Scenario reclassification + domain audit", FALSE,
    "TABLE_scenarios_tidy.csv", "Scenario tidy summary", FALSE,
    
    "scenario_sim_*.csv", "Scenario raw sim outputs", FALSE,
    "scenario_timeseries_*.png", "Scenario multi-panel timeseries", FALSE,
    "phase_end_*.png", "Phase plots with start/end + arrows", FALSE
  )
  
  expected %>%
    mutate(
      path = fs::path(out_dir, file),
      exists = if_else(str_detect(file, "\\*"), TRUE, fs::file_exists(path)),
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
## 9) BASELINE PARAMS
## ============================================================
par_base <- list(
  # growth / tech
  alpha = 0.02, beta = 0.01,
  delta = 0.02,
  sigma = 3.0,
  u_star = 1.0,
  
  # finance
  i = 0.04,
  
  # investment logistic placeholders (time-varying uses *_pre/post)
  kappa_min = 0.00, kappa_max = 0.30,
  lambda = 20,   # overwritten in reduced analysis by plane values
  r0 = 0.04,     # overwritten similarly
  
  # Phillips
  phi0 = -0.06,
  phi1 = 0.10,
  
  # residual u adjustment
  lambda_u = 2.0,
  
  # wealth accounting param (DO NOT use "p" as a state)
  p = 1.0,
  
  # autonomous consumption (calibrated)
  varsigma0 = 0.55,
  varsigma1 = 0.00,
  
  # induced non-owner consumption
  cN_minus = 0.02,
  cN_plus  = 0.20,
  cN_eta   = 0.10,
  
  # owners saving
  sA_const  = TRUE,
  sA0       = 0.40,
  sA_min    = 0.20, sA_max = 0.80, sA_lambda = 2.0, sA_mid = 2.0,
  
  # leakage psi(d) will be turned on in regime
  psi_zero  = TRUE,
  psi_max   = 0.00, psi_lambda = 5.0, psi_mid = 1.0,
  
  # wage share floor
  omega_floor = 0.05,
  
  # regime patriarch controls
  lambda_pre  = 51.0,
  lambda_post = 72.0,
  r0_pre      = 0.040,
  r0_post     = 0.065,
  psi_max_pre  = 0.00,
  psi_max_post = 0.08,
  
  t_crash = 350,
  k_R     = 0.03,
  k_mu    = 80,
  
  # stress threshold
  d_bar = 1.05,
  stress_s0 = 0.00,
  rat_s0    = 0.00,
  
  # asset price process
  g_p0 = 0.00,
  g_pR = 0.015,
  g_pS = 0.10,
  g_pRgap = 0.30,
  
  # inequality channels + ratchet
  k_xi_p_pre  = 0.00,
  k_xi_p_post = 0.20,
  k_xi_d_pre  = 0.02,
  k_xi_d_post = 0.10,
  k_xi_mr_pre = 0.06,
  k_xi_mr_post= 0.01,
  
  eta_rat = 0.08,
  delta_floor = 0.005,
  xi_floor_base = 0.60,
  
  # initial conditions for new states
  pA0 = 1.0,
  xi0 = 1.0,
  xi_floor0 = 0.60
)

ranges <- list(
  phi1   = c(0.03, 0.25),
  lambda = c(5, 150),
  i      = c(0.00, 0.10),
  phi0   = c(-0.12, -0.01),
  r0     = c(0.01, 0.08)
)

## ============================================================
## 10) DRIVER A: Plane scores (same planes you already like)
## ============================================================
cat("[A] Scoring candidate Hopf planes...\n")

planes <- list(
  c("phi1","lambda"),
  c("i","lambda"),
  c("r0","lambda"),
  c("phi0","phi1"),
  c("phi1","i")
)

compare_planes <- function(par_base, planes, ranges, n = 80,
                           H_band = 1e-3, tol = 1e-6,
                           calibrate_varsigma0 = TRUE) {
  
  out <- vector("list", length(planes))
  
  for (j in seq_along(planes)) {
    p1 <- planes[[j]][1]
    p2 <- planes[[j]][2]
    
    p1_rng <- ranges[[p1]]
    p2_rng <- ranges[[p2]]
    
    df <- hopf_plane(
      par_base,
      p1 = p1, p2 = p2,
      p1_vals = grid_values(p1_rng[1], p1_rng[2], n),
      p2_vals = grid_values(p2_rng[1], p2_rng[2], n),
      calibrate_varsigma0 = calibrate_varsigma0,
      H_band = H_band
    )
    
    sc <- score_hopf_plane(df, H_band = H_band, tol = tol) %>% mutate(p1 = p1, p2 = p2)
    out[[j]] <- sc
  }
  
  bind_rows(out) %>%
    arrange(desc(has_robust_sign_change), degeneracy_flag, desc(feasible_share), desc(sd_H))
}

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
plane_outputs <- make_plane_ranking_outputs(scores, out_dir = OUT_DIR)

cat("[A] saved: hopf_plane_scores.csv + TABLE_plane_ranking_tidy.csv + PLOT_plane_selection_scatter.png\n")

## ============================================================
## 11) DRIVER B/C/D: Hopf plane figures
## ============================================================
cat("[B] Exporting Hopf plane maps...\n")

# FIG A: (phi1, lambda)
df_A <- hopf_plane(
  par_base,
  p1 = "phi1", p2 = "lambda",
  p1_vals = grid_values(ranges$phi1[1], ranges$phi1[2], 140),
  p2_vals = grid_values(ranges$lambda[1], ranges$lambda[2], 140),
  calibrate_varsigma0 = TRUE,
  H_band = 1e-3
)
write_csv_quiet(df_A, fs::path(OUT_DIR, "hopf_plane_phi1_lambda.csv"))
pA <- plot_hopf_plane_signH(df_A, "phi1", "lambda", title = "FIG A: Hopf map (lambda vs phi1) with sign(H)")
ggsave(fs::path(OUT_DIR, "FIG_A_hopf_plane_signH_phi1_lambda.png"), pA, width = 7.2, height = 6.0, dpi = 160)

# FIG B: (r0, lambda)
df_B <- hopf_plane(
  par_base,
  p1 = "r0", p2 = "lambda",
  p1_vals = grid_values(ranges$r0[1], ranges$r0[2], 140),
  p2_vals = grid_values(ranges$lambda[1], ranges$lambda[2], 140),
  calibrate_varsigma0 = TRUE,
  H_band = 1e-3
)
write_csv_quiet(df_B, fs::path(OUT_DIR, "hopf_plane_r0_lambda.csv"))
pB <- plot_hopf_plane_signH(df_B, "r0", "lambda", title = "FIG B: Hopf map (lambda vs r0) with sign(H)")
ggsave(fs::path(OUT_DIR, "FIG_B_hopf_plane_signH_r0_lambda.png"), pB, width = 7.2, height = 6.0, dpi = 160)

# FIG C: (phi1, i)
df_C <- hopf_plane(
  par_base,
  p1 = "phi1", p2 = "i",
  p1_vals = grid_values(ranges$phi1[1], ranges$phi1[2], 140),
  p2_vals = grid_values(ranges$i[1], ranges$i[2], 140),
  calibrate_varsigma0 = TRUE,
  H_band = 1e-3
)
write_csv_quiet(df_C, fs::path(OUT_DIR, "hopf_plane_phi1_i.csv"))
pC <- plot_hopf_plane_signH(df_C, "phi1", "i", title = "FIG C: Hopf map (i vs phi1) with sign(H)")
ggsave(fs::path(OUT_DIR, "FIG_C_hopf_plane_signH_phi1_i.png"), pC, width = 7.2, height = 6.0, dpi = 160)

cat("[B] saved Hopf planes A/B/C (CSVs + PNGs)\n")

## ============================================================
## 12) DRIVER E: H audit table + H quantile plot
## ============================================================
cat("[C] Auditing H distributions across exported planes...\n")

plane_df_list <- list(
  hopf_plane_phi1_lambda = df_A,
  hopf_plane_r0_lambda   = df_B,
  hopf_plane_phi1_i      = df_C
)

plane_audit <- bind_rows(lapply(names(plane_df_list), function(nm) {
  plane_H_audit(plane_df_list[[nm]], nm)
}))

write_csv_quiet(plane_audit, fs::path(OUT_DIR, "TABLE_plane_H_audit.csv"))

pH <- plane_audit %>%
  mutate(plane = factor(plane, levels = plane)) %>%
  ggplot(aes(x = plane, y = H_med)) +
  geom_point() +
  geom_errorbar(aes(ymin = H_q05, ymax = H_q95), width = 0.15) +
  coord_flip() +
  labs(title = "Plane audit: H median with 5–95% interval (within feasible RH region)",
       x = NULL, y = "H")

ggsave(fs::path(OUT_DIR, "PLOT_plane_H_quantiles.png"), pH, width = 8.5, height = 5.2, dpi = 160)

cat("[C] saved: TABLE_plane_H_audit.csv + PLOT_plane_H_quantiles.png\n")

## ============================================================
## 13) DRIVER F: Find lambda* for a chosen phi1
## ============================================================
cat("[D] Finding lambda* (Hopf root) for phi1_pick...\n")

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

cat(sprintf("\n[D] Hopf approx: lambda* = %.6f at phi1 = %.4f\n", lambda_star, phi1_pick))

# Store in par_base for regime simulation
par_base$lambda_star <- lambda_star

## ============================================================
## 14) DRIVER G: Scenarios (fixed and moving Hopf patriarch)
## ============================================================
cat("[E] Running scenarios...\n")

scenario_tbl <- list()

# ---- G1) fixed lambda below/near/above lambda*
lambda_vals <- c(
  below = 0.70 * lambda_star,
  near  = 0.98 * lambda_star,
  above = 1.20 * lambda_star
)
lambda_vals <- pmin(pmax(lambda_vals, ranges$lambda[1]), ranges$lambda[2])

for (nm in names(lambda_vals)) {
  
  par_s <- par_base
  
  # fixed lambda: pre=post, so moving-Hopffiness is turned off
  par_s$lambda_pre  <- as.numeric(lambda_vals[[nm]])
  par_s$lambda_post <- as.numeric(lambda_vals[[nm]])
  
  # keep r0 drift OFF for fixed baseline comparison
  par_s$r0_pre  <- par_base$r0_pre
  par_s$r0_post <- par_base$r0_pre
  
  tag <- paste0("fixed_phi1_", sprintf("%.3f", phi1_pick),
                "_lambda_", nm, "_", sprintf("%.2f", par_s$lambda_pre))
  
  df_s <- simulate_regime(par_s, t_end = 900, dt = 0.05)
  
  write_csv_quiet(df_s, fs::path(OUT_DIR, paste0("scenario_sim_", tag, ".csv")))
  plot_scenario_timeseries(df_s, tag, out_dir = OUT_DIR)
  plot_scenario_phase_end(df_s, tag, out_dir = OUT_DIR)
  
  scenario_tbl[[tag]] <- summarize_tail(df_s, name = tag, frac_tail = 0.25) %>%
    mutate(phi1 = phi1_pick,
           lambda_pre = par_s$lambda_pre,
           lambda_post = par_s$lambda_post,
           lambda_star = lambda_star,
           lambda_rel = par_s$lambda_pre / lambda_star,
           scenario_type = "fixed_lambda",
           sim_file = fs::path(OUT_DIR, paste0("scenario_sim_", tag, ".csv")),
           sim_exists = TRUE)
}

# ---- G2) moving Hopf patriarch: lambda drifts from pre->post across t_crash
par_m <- par_base
par_m$lambda_pre  <- pmin(pmax(0.85 * lambda_star, ranges$lambda[1]), ranges$lambda[2])
par_m$lambda_post <- pmin(pmax(1.20 * lambda_star, ranges$lambda[1]), ranges$lambda[2])

tag_m <- paste0("movingHopf_phi1_", sprintf("%.3f", phi1_pick),
                "_lamPre_", sprintf("%.2f", par_m$lambda_pre),
                "_lamPost_", sprintf("%.2f", par_m$lambda_post),
                "_dbar_", sprintf("%.2f", par_m$d_bar))

df_m <- simulate_regime(par_m, t_end = 900, dt = 0.05)

write_csv_quiet(df_m, fs::path(OUT_DIR, paste0("scenario_sim_", tag_m, ".csv")))
plot_scenario_timeseries(df_m, tag_m, out_dir = OUT_DIR)
plot_scenario_phase_end(df_m, tag_m, out_dir = OUT_DIR)

scenario_tbl[[tag_m]] <- summarize_tail(df_m, name = tag_m, frac_tail = 0.25) %>%
  mutate(phi1 = phi1_pick,
         lambda_pre = par_m$lambda_pre,
         lambda_post = par_m$lambda_post,
         lambda_star = lambda_star,
         lambda_rel = par_m$lambda_pre / lambda_star,
         scenario_type = "moving_hopf",
         sim_file = fs::path(OUT_DIR, paste0("scenario_sim_", tag_m, ".csv")),
         sim_exists = TRUE)

scenario_tbl <- bind_rows(scenario_tbl)

write_csv_quiet(scenario_tbl, fs::path(OUT_DIR, "TABLE_scenarios_tail_diagnostics.csv"))

# scenario tidy (simple flags)
scenario_tidy <- scenario_tbl %>%
  mutate(
    cycle_flag = is.finite(period_e) & period_e > 0 & is.finite(e_amp) & e_amp > 1e-3
  ) %>%
  select(run, time_end,
         e_mean, e_sd, e_amp,
         omega_mean, omega_sd, omega_amp,
         d_mean, d_sd, d_amp,
         u_mean, u_sd, u_amp,
         pA_mean, pA_sd, pA_amp,
         xi_mean, xi_sd, xi_amp,
         gY_mean, period_e,
         phi1, lambda_pre, lambda_post, lambda_star, lambda_rel,
         cycle_flag, scenario_type)

write_csv_quiet(scenario_tidy, fs::path(OUT_DIR, "TABLE_scenarios_tidy.csv"))

# scenario audit (reclassification from full-path + amplitude thresholds)
scen_audit <- scenario_tidy %>%
  rowwise() %>%
  mutate(
    sim_file = fs::path(OUT_DIR, paste0("scenario_sim_", run, ".csv")),
    sim_exists = fs::file_exists(sim_file),
    regime2 = if (sim_exists) {
      df <- readr::read_csv(sim_file, show_col_types = FALSE)
      classify_scenario(df)
    } else NA_character_
  ) %>%
  ungroup()

write_csv_quiet(scen_audit, fs::path(OUT_DIR, "TABLE_scenarios_audit_reclassified.csv"))

cat("[E] saved: TABLE_scenarios_tail_diagnostics.csv + TABLE_scenarios_tidy.csv + TABLE_scenarios_audit_reclassified.csv\n")

## ============================================================
## 15) DRIVER H: Scenario summary plots (amplitude + period)
## ============================================================
cat("[F] Exporting scenario summary plots...\n")

scen_long <- scenario_tidy %>%
  select(run, lambda_rel, scenario_type, cycle_flag,
         e_amp, omega_amp, d_amp, u_amp, pA_amp, xi_amp, period_e) %>%
  pivot_longer(cols = c(e_amp, omega_amp, d_amp, u_amp, pA_amp, xi_amp),
               names_to = "metric", values_to = "value")

p_amp <- ggplot(scen_long, aes(x = lambda_rel, y = value, label = run, shape = scenario_type)) +
  geom_point(size = 3) +
  geom_text(check_overlap = TRUE, nudge_y = 0.002, size = 3) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  labs(
    title = "Tail amplitudes vs λ_pre / λ* (scenario comparison)",
    x = "lambda_rel = λ_pre / λ*",
    y = "tail amplitude (q95 - q05)",
    shape = "scenario type"
  )

ggsave(fs::path(OUT_DIR, "PLOT_scenarios_amplitudes_vs_lambdaRel.png"),
       p_amp, width = 8.2, height = 10.0, dpi = 160)

p_per <- ggplot(scenario_tidy, aes(x = lambda_rel, y = period_e, label = run, shape = scenario_type)) +
  geom_point(size = 3) +
  geom_text(check_overlap = TRUE, nudge_y = 0.05, size = 3) +
  labs(
    title = "Period proxy: peak-to-peak period of e(t) in tail window",
    x = "lambda_rel = λ_pre / λ*",
    y = "period_e (median peak-to-peak distance)",
    shape = "scenario type"
  )

ggsave(fs::path(OUT_DIR, "PLOT_scenarios_period_vs_lambdaRel.png"),
       p_per, width = 8.2, height = 5.8, dpi = 160)

cat("[F] saved: PLOT_scenarios_amplitudes_vs_lambdaRel.png + PLOT_scenarios_period_vs_lambdaRel.png\n")

## ============================================================
## 16) DRIVER I: Manifest
## ============================================================
manifest <- make_outputs_manifest(OUT_DIR)
write_csv_quiet(manifest, fs::path(OUT_DIR, "TABLE_outputs_review_manifest.csv"))

cat("\nDONE. Exports in:\n", OUT_DIR, "\n\nKey files:\n",
    "- hopf_plane_scores.csv\n",
    "- TABLE_plane_ranking_tidy.csv\n",
    "- PLOT_plane_selection_scatter.png\n",
    "- hopf_plane_phi1_lambda.csv + FIG_A_hopf_plane_signH_phi1_lambda.png\n",
    "- hopf_plane_r0_lambda.csv  + FIG_B_hopf_plane_signH_r0_lambda.png\n",
    "- hopf_plane_phi1_i.csv     + FIG_C_hopf_plane_signH_phi1_i.png\n",
    "- TABLE_plane_H_audit.csv + PLOT_plane_H_quantiles.png\n",
    "- TABLE_scenarios_tail_diagnostics.csv\n",
    "- TABLE_scenarios_tidy.csv\n",
    "- TABLE_scenarios_audit_reclassified.csv\n",
    "- scenario_sim_*.csv + scenario_timeseries_*.png + phase_end_*.png\n",
    "- PLOT_scenarios_amplitudes_vs_lambdaRel.png\n",
    "- PLOT_scenarios_period_vs_lambdaRel.png\n",
    "- TABLE_outputs_review_manifest.csv\n\n")
