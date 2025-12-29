############################################################
## wealth_goodwin_hopf_map_residual_u_RUNALL.R
##
## Purpose (single script, all outputs):
##  - Residual follower u: u adjusts to u_d; no feedback into (e, omega, d).
##  - Enforce u_d = u_star at SS via calibrating varsigma0.
##  - Reduced system Hopf mapping (SS -> J* -> RH -> H) across planes.
##  - Rank planes (ROBUST scoring), export FIG A/FIG B Hopf maps with sign(H).
##  - Find lambda* and simulate below/near/above scenarios.
##  - Export scenario sims, timeseries panels, PHASE plots with start/end labels,
##    tail diagnostics table, and scenario amplitude/period summary plots.
##  - Export review manifest + plane H-audit + scenario reclassification audit.
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
## GLOBAL OUTPUT DIR (ALL EXPORTS GO HERE)
## ============================================================
OUT_DIR <- fs::path("outputs", "wealth_goodwin", "grid_search")
fs::dir_create(OUT_DIR)

## ============================================================
## 0) Helpers
## ============================================================
inv_logit <- function(x) 1 / (1 + exp(-x))

require_params <- function(par, keys) {
  miss <- setdiff(keys, names(par))
  if (length(miss) > 0) stop("Missing parameters: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

`%||%` <- function(a, b) if (is.null(a) || !is.finite(a)) b else a

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

## ============================================================
## Phillips + profit
## ============================================================
Phi_fun    <- function(e, par) par$phi0 + par$phi1 * e
Phi_e_fun  <- function(e, par) par$phi1
pi_fun     <- function(omega) 1 - omega

## Reduced net profit rate (u fixed at u_star)
r_net_reduced <- function(omega, d, par) {
  (par$u_star / par$sigma) * (pi_fun(omega) - par$i * d)
}

## Logistic investment share
kappa_fun <- function(r, par) {
  par$kappa_min + (par$kappa_max - par$kappa_min) /
    (1 + exp(-par$lambda * (r - par$r0)))
}
kappa_r_fun <- function(r, par) {
  s <- inv_logit(par$lambda * (r - par$r0))
  (par$kappa_max - par$kappa_min) * par$lambda * s * (1 - s)
}

## Owners saving
sA_fun <- function(chiA, par) {
  if (isTRUE(par$sA_const)) return(par$sA0)
  s <- inv_logit(par$sA_lambda * (chiA - par$sA_mid))
  par$sA_min + (par$sA_max - par$sA_min) * s
}

## Leakage psi(d)
psi_fun <- function(d, par) {
  if (isTRUE(par$psi_zero)) return(0)
  inv_logit(par$psi_lambda * (d - par$psi_mid)) * par$psi_max
}

## Non-owners induced consumption "share"
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
## 1) Reduced steady state (u = u_star)
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
  
  ## d* identity
  d_star <- (k_target - (1 / A) * r_star) / g_n
  if (d_star <= eps) return(list(ok = FALSE, reason = "d* <= 0 (no interior debt steady state)"))
  
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
## 2) Enforce u_d = u_star at SS by calibrating varsigma0
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
## 3) Reduced Jacobian (e, omega, d) at SS + RH + H
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
  eig <- eigen(Jp$J)$values
  
  list(
    ok = TRUE,
    par = par,
    steady_state = ss,
    jacobian = Jp$J,
    eigenvalues = eig,
    rh = rh,
    jac_parts = list(kappa_r = Jp$kappa_r, Phi_e = Jp$Phi_e)
  )
}

## ============================================================
## Hopf gap + robust root finder (scan -> bracket -> uniroot)
## ============================================================
hopf_gap_reduced <- function(par, calibrate_varsigma0 = TRUE) {
  out <- analyze_reduced(par, calibrate_varsigma0 = calibrate_varsigma0)
  if (!isTRUE(out$ok)) return(NA_real_)
  out$rh$H
}

find_hopf_mu <- function(par, mu, lo, hi,
                         n_grid = 300,
                         tol = 1e-8, maxit = 200,
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
## 4) Hopf mapping utilities on a 2D parameter plane
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

## Robust plane scoring (fixes the “phi1-i degeneracy gets #1” problem)
score_hopf_plane_v2 <- function(df, H_band = 1e-3, tol = 1e-6) {
  
  n_ok <- sum(df$ok, na.rm = TRUE)
  
  dfF <- df %>% filter(ok, feasible_RH, is.finite(H))
  if (nrow(dfF) == 0) {
    return(tibble(
      n_ok = n_ok,
      n_feasible = 0,
      feasible_share = 0,
      hopf_band_share = NA_real_,
      share_H_pos = NA_real_,
      share_H_neg = NA_real_,
      H_min = NA_real_,
      H_max = NA_real_,
      H_med = NA_real_,
      H_q05 = NA_real_,
      H_q95 = NA_real_,
      sd_H = NA_real_,
      sign_change = FALSE,
      has_robust_sign_change = FALSE,
      degeneracy_flag = TRUE
    ))
  }
  
  Hmin <- min(dfF$H)
  Hmax <- max(dfF$H)
  
  share_pos <- mean(dfF$H >  tol)
  share_neg <- mean(dfF$H < -tol)
  
  has_robust_sign_change <- (Hmin < -tol) && (Hmax > tol)
  sign_change <- has_robust_sign_change
  
  sdH <- stats::sd(dfF$H)
  degeneracy_flag <- (sdH < 1e-8) || (!has_robust_sign_change)
  
  tibble(
    n_ok = n_ok,
    n_feasible = nrow(dfF),
    feasible_share = nrow(dfF) / max(1, n_ok),
    hopf_band_share = mean(abs(dfF$H) < H_band),
    share_H_pos = share_pos,
    share_H_neg = share_neg,
    H_min = Hmin,
    H_max = Hmax,
    H_med = median(dfF$H),
    H_q05 = as.numeric(quantile(dfF$H, 0.05)),
    H_q95 = as.numeric(quantile(dfF$H, 0.95)),
    sd_H = sdH,
    sign_change = sign_change,
    has_robust_sign_change = has_robust_sign_change,
    degeneracy_flag = degeneracy_flag
  )
}

compare_planes <- function(par_base, planes, ranges, n = 80,
                           H_band = 1e-3, tol = 1e-6,
                           calibrate_varsigma0 = TRUE) {
  
  out <- vector("list", length(planes))
  
  for (j in seq_along(planes)) {
    p1 <- planes[[j]][1]
    p2 <- planes[[j]][2]
    
    p1_rng <- ranges[[p1]]
    p2_rng <- ranges[[p2]]
    stopifnot(length(p1_rng) == 2, length(p2_rng) == 2)
    
    df <- hopf_plane(
      par_base,
      p1 = p1, p2 = p2,
      p1_vals = grid_values(p1_rng[1], p1_rng[2], n),
      p2_vals = grid_values(p2_rng[1], p2_rng[2], n),
      calibrate_varsigma0 = calibrate_varsigma0,
      H_band = H_band
    )
    
    sc <- score_hopf_plane_v2(df, H_band = H_band, tol = tol) %>%
      mutate(p1 = p1, p2 = p2)
    
    out[[j]] <- sc
  }
  
  bind_rows(out) %>%
    arrange(
      desc(has_robust_sign_change),
      degeneracy_flag,
      desc(feasible_share),
      desc(sd_H))
  
}

plot_hopf_plane_signH <- function(df, p1, p2, title = NULL) {
  
  df2 <- df %>%
    mutate(
      status = case_when(
        !ok ~ "SS infeasible",
        ok & !feasible_RH ~ "RH prereq fail (a1/a2/a3)",
        ok & feasible_RH & is.finite(H) & H > 0 ~ "H>0 (stable focus)",
        ok & feasible_RH & is.finite(H) & H < 0 ~ "H<0 (unstable)",
        ok & feasible_RH & is.finite(H) & abs(H) <= 1e-10 ~ "H≈0",
        TRUE ~ "Other"
      ),
      status = factor(status, levels = c(
        "SS infeasible",
        "RH prereq fail (a1/a2/a3)",
        "H>0 (stable focus)",
        "H<0 (unstable)",
        "H≈0",
        "Other"
      ))
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

## ============================================================
## 5) Full simulation with RESIDUAL u (no feedback into e,omega,d)
## ============================================================
rhs_full_residual_u <- function(t, state, par) {
  e <- state[["e"]]
  w <- state[["omega"]]
  d <- state[["d"]]
  u <- state[["u"]]
  
  ## guardrails
  e <- min(max(e, 1e-9), 1 - 1e-9)
  w <- min(max(w, 1e-9), 1 - 1e-9)
  d <- max(d, 1e-12)
  u <- max(u, 1e-9)
  
  ## core at u_star
  r   <- r_net_reduced(w, d, par)
  kap <- kappa_fun(r, par)
  
  gK <- (par$u_star * kap) / par$sigma - par$delta
  gY <- gK
  
  g_n <- par$alpha + par$beta
  
  Phi <- Phi_fun(e, par)
  de  <- (gY - g_n) * e
  dw  <- (Phi - par$alpha) * w
  dd  <- kap - (1 - w) + par$i * d - d * gY
  
  ## residual follower u
  u_d <- u_desired_uStar(t, w, d, par)
  if (!is.finite(u_d)) {
    du <- 0
    gU <- 0
  } else {
    du <- par$lambda_u * (u_d - u)
    gU <- du / u
  }
  
  ## diagnostics for inequality (uses realized u)
  rho  <- par$sigma / u
  chiA <- par$p * rho + d
  chiN <- -d
  xi   <- chiA / abs(chiN)
  
  Delta <- Delta_fun_uStar(t, w, d, par)
  
  list(
    c(de, dw, dd, du),
    c(r = r, kappa = kap, gK = gK, gY = gY,
      Delta = Delta, u_d = u_d, gU = gU,
      rho = rho, chiA = chiA, chiN = chiN, xi = xi)
  )
}

simulate_full_residual_u <- function(par, state0 = NULL,
                                     t_end = 700, dt = 0.05,
                                     perturb = c(e = 0.01, omega = -0.01, d = 0.05, u = -0.05),
                                     calibrate_varsigma0 = TRUE) {
  
  if (calibrate_varsigma0) {
    cal <- calibrate_varsigma0_at_ss(par)
    if (!isTRUE(cal$ok)) stop("Calibration failed: ", cal$reason)
    par <- cal$par
    ss  <- cal$ss
  } else {
    ss <- steady_state_reduced(par)
    if (!isTRUE(ss$ok)) stop("No admissible reduced SS: ", ss$reason)
  }
  
  if (is.null(state0)) {
    state0 <- c(
      e     = ss$e + perturb[["e"]],
      omega = ss$omega + perturb[["omega"]],
      d     = ss$d + perturb[["d"]],
      u     = par$u_star + perturb[["u"]]
    )
  }
  
  times <- seq(0, t_end, by = dt)
  
  sol <- ode(
    y = state0,
    times = times,
    func = rhs_full_residual_u,
    parms = par,
    method = "lsoda"
  )
  
  as.data.frame(sol) |> as_tibble()
}

## ============================================================
## 6) Tail diagnostics
## ============================================================
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

period_peaks <- function(t, x, min_prom = NULL) {
  x <- as.numeric(x); t <- as.numeric(t)
  ok <- is.finite(x) & is.finite(t)
  x <- x[ok]; t <- t[ok]
  if (length(x) < 200) return(NA_real_)
  
  dx1 <- diff(x)
  peaks <- which(dx1[-1] < 0 & dx1[-length(dx1)] > 0) + 1
  if (length(peaks) < 3) return(NA_real_)
  
  if (!is.null(min_prom)) {
    keep <- x[peaks] > (median(x) + min_prom)
    peaks <- peaks[keep]
    if (length(peaks) < 3) return(NA_real_)
  }
  
  median(diff(t[peaks]))
}

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
    
    xi_mean = mean(dft$xi, na.rm = TRUE),
    xi_sd   = sd(dft$xi, na.rm = TRUE),
    xi_amp  = amp(log(pmax(dft$xi, 1e-12))),
    
    gY_mean = mean(dft$gY, na.rm = TRUE),
    
    period_e = period_peaks(dft$time, dft$e)
  )
}

## ============================================================
## 7) Plotting: timeseries panels + phase plots with START/END labels
## ============================================================
plot_scenario_timeseries <- function(df, tag, out_dir = OUT_DIR) {
  p <- df %>%
    select(time, e, omega, d, u, xi) %>%
    pivot_longer(-time, names_to = "var", values_to = "value") %>%
    ggplot(aes(x = time, y = value)) +
    geom_line() +
    facet_wrap(~var, scales = "free_y", ncol = 1) +
    labs(title = paste0("Scenario time paths (", tag, ")"), x = "t", y = NULL)
  
  ggsave(fs::path(out_dir, paste0("scenario_timeseries_", tag, ".png")),
         p, width = 7.2, height = 9.0, dpi = 160)
  p
}

phase_end_plot <- function(dat, x, y, tag,
                           out_dir = OUT_DIR,
                           tail_frac = 0.25,
                           arrow_every = 40) {
  
  dat <- dplyr::arrange(dat, time)
  
  n <- nrow(dat)
  i0 <- max(1, floor((1 - tail_frac) * n))
  dat <- dat %>%
    dplyr::mutate(idx = dplyr::row_number(),
                  is_tail = idx >= i0)
  
  p_start <- dat[1, ]
  p_end   <- dat[n, ]
  
  ## arrows along tail
  tail_dat <- dat %>% dplyr::filter(is_tail) %>%
    dplyr::mutate(
      x0 = .data[[x]], y0 = .data[[y]],
      x1 = dplyr::lead(.data[[x]]), y1 = dplyr::lead(.data[[y]])
    )
  
  arrows <- tail_dat %>% dplyr::filter(!is.na(x1), idx %% arrow_every == 0)
  
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
  p1 <- phase_end_plot(df, x = "e", y = "omega", tag = tag, out_dir = out_dir)
  p2 <- phase_end_plot(df, x = "d", y = "xi",    tag = tag, out_dir = out_dir)
  invisible(list(p_omega_e = p1, p_xi_d = p2))
}

## ============================================================
## 8) Post-processing: tidy tables + selection plots
## ============================================================
make_plane_decision_outputs <- function(scores, out_dir = OUT_DIR) {
  
  # safe numeric cleaner for columns (NA/Inf -> fallback)
  clean_num <- function(x, fallback = 0) {
    x <- tidyr::replace_na(x, fallback)
    x[!is.finite(x)] <- fallback
    x
  }
  
  plane_tbl <- scores %>%
    mutate(
      plane = paste0("(", p1, ", ", p2, ")"),
      feasible_share   = clean_num(feasible_share, 0),
      hopf_band_share  = clean_num(hopf_band_share, 0),
      sd_H             = clean_num(sd_H, 0),
      H_min            = clean_num(H_min, NA_real_),
      H_max            = clean_num(H_max, NA_real_),
      
      # Simple composite score: reward feasible area + near-H density,
      # reward robust sign-change, penalize degeneracy.
      score = 100 * feasible_share +
        25  * hopf_band_share +
        20  * as.integer(has_robust_sign_change) -
        30  * as.integer(degeneracy_flag)
    ) %>%
    arrange(desc(has_robust_sign_change), degeneracy_flag, desc(score)) %>%
    select(
      plane, p1, p2,
      has_robust_sign_change, degeneracy_flag, sign_change,
      feasible_share, hopf_band_share, sd_H,
      H_min, H_max,
      n_feasible, n_ok, score
    )
  
  readr::write_csv(plane_tbl, fs::path(out_dir, "TABLE_plane_ranking_tidy.csv"))
  
  p_sel <- plane_tbl %>%
    ggplot(aes(x = feasible_share, y = hopf_band_share, label = plane)) +
    geom_point(aes(shape = has_robust_sign_change), size = 3) +
    geom_text(check_overlap = TRUE, nudge_y = 0.01, size = 3) +
    labs(
      title = "Plane selection: feasibility vs near-Hopf thickness",
      subtitle = "Shape = robust sign-change (Hopf candidate). Higher is better.",
      x = "feasible_share (RH feasible among SS-ok)",
      y = "hopf_band_share (|H| < band within feasible RH)"
    )
  
  ggsave(fs::path(out_dir, "PLOT_plane_selection_scatter.png"),
         p_sel, width = 9.2, height = 6.0, dpi = 160)
  
  list(plane_tbl = plane_tbl, plot = p_sel)
}


make_scenario_outputs <- function(scen_tbl, out_dir = OUT_DIR) {
  scen_tbl2 <- scen_tbl %>%
    mutate(
      lambda_rel = lambda / lambda_star,
      cycle_flag = (is.finite(period_e) & period_e > 0 & is.finite(e_amp) & e_amp > 1e-3)
    )
  
  write_csv(scen_tbl2, fs::path(out_dir, "TABLE_scenarios_tidy.csv"))
  
  scen_long <- scen_tbl2 %>%
    select(run, lambda_rel, e_amp, d_amp, xi_amp, period_e, cycle_flag) %>%
    pivot_longer(cols = c(e_amp, d_amp, xi_amp), names_to = "metric", values_to = "value")
  
  p_amp <- ggplot(scen_long, aes(x = lambda_rel, y = value, label = run, shape = cycle_flag)) +
    geom_point(size = 3) +
    geom_text(check_overlap = TRUE, nudge_y = 0.002, size = 3) +
    facet_wrap(~metric, scales = "free_y", ncol = 1) +
    labs(
      title = "Cycle emergence diagnostic: tail amplitudes vs λ/λ*",
      x = "lambda_rel = λ / λ*",
      y = "tail amplitude (q95 - q05)",
      shape = "cycle_flag"
    )
  
  ggsave(fs::path(out_dir, "PLOT_scenarios_amplitudes_vs_lambdaRel.png"),
         p_amp, width = 7.6, height = 9.0, dpi = 160)
  
  p_per <- ggplot(scen_tbl2, aes(x = lambda_rel, y = period_e, label = run, shape = cycle_flag)) +
    geom_point(size = 3) +
    geom_text(check_overlap = TRUE, nudge_y = 0.02, size = 3) +
    labs(
      title = "Period proxy: peak-to-peak period of e(t) in tail window",
      x = "lambda_rel = λ / λ*",
      y = "period_e (median peak-to-peak distance)",
      shape = "cycle_flag"
    )
  
  ggsave(fs::path(out_dir, "PLOT_scenarios_period_vs_lambdaRel.png"),
         p_per, width = 8.2, height = 5.5, dpi = 160)
  
  list(scen_tbl = scen_tbl2, p_amp = p_amp, p_per = p_per)
}

## ============================================================
## 9) Audit exports: outputs manifest + plane H audit + scenario reclass audit
## ============================================================
make_outputs_manifest <- function(out_dir = OUT_DIR,
                                  reviewed_figA = FALSE,
                                  reviewed_figB = FALSE,
                                  reviewed_Haudit = FALSE) {
  
  expected <- tribble(
    ~file, ~purpose, ~reviewed,
    "hopf_plane_scores.csv", "Raw plane ranking metrics (ROBUST fields)", FALSE,
    "TABLE_plane_ranking_tidy.csv", "Tidy plane ranking table", FALSE,
    "PLOT_plane_selection_scatter.png", "Plane selection scatter plot", FALSE,
    "hopf_plane_phi1_lambda.csv", "Plane grid points for FIG A", FALSE,
    "hopf_plane_phi1_i.csv", "Plane grid points for FIG B (if chosen)", FALSE,
    "FIG_A_hopf_plane_signH_phi1_lambda.png", "Hopf map A (phi1, lambda)", reviewed_figA,
    "FIG_B_hopf_plane_signH_<...>.png", "Hopf map B (chosen plane)", reviewed_figB,
    "TABLE_scenarios_tail_diagnostics.csv", "Tail-window diagnostics per scenario", FALSE,
    "TABLE_scenarios_tidy.csv", "Scenario table (lambda_rel, cycle_flag)", FALSE,
    "PLOT_scenarios_amplitudes_vs_lambdaRel.png", "Amplitudes vs lambda_rel", FALSE,
    "PLOT_scenarios_period_vs_lambdaRel.png", "Period vs lambda_rel", FALSE,
    "phase_end_*", "Phase plots with start/end labels", FALSE,
    "TABLE_plane_H_audit.csv", "Plane audit: H quantiles (feasible RH only)", FALSE,
    "PLOT_plane_H_quantiles.png", "Plane audit plot: H median + 5–95%", reviewed_Haudit,
    "TABLE_scenarios_audit_reclassified.csv", "Scenario audit: whole-path domain + amp thresholds", FALSE
  )
  
  # check existence with a basic wildcard support
  expected %>%
    mutate(
      exists = case_when(
        grepl("\\*$", file) ~ TRUE, # wildcard: assume generated if you asked for it
        grepl("<\\.\\.\\.>", file) ~ TRUE, # placeholder: treated as “expected”
        TRUE ~ fs::file_exists(fs::path(out_dir, file))
      ),
      status = case_when(
        reviewed ~ "Reviewed",
        exists & !reviewed ~ "Ready to review",
        !exists ~ "Missing",
        TRUE ~ "Unknown"
      )
    ) %>%
    select(file, exists, reviewed, status, purpose)
}

audit_plane_files <- function(out_dir = OUT_DIR) {
  csv_files <- fs::dir_ls(out_dir, regexp = "\\.csv$")
  bn <- basename(csv_files)
  plane_files <- csv_files[grepl("^hopf_plane_", bn) & !grepl("scores", bn)]
  
  if (length(plane_files) == 0) {
    write_csv(tibble(), fs::path(out_dir, "TABLE_plane_H_audit.csv"))
    return(invisible(NULL))
  }
  
  audit_one <- function(df, plane_name) {
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
  
  audits <- lapply(plane_files, function(f) {
    df <- read_csv_safe(f)
    if (is.null(df)) return(NULL)
    nm <- tools::file_path_sans_ext(basename(f))
    audit_one(df, nm)
  }) %>% bind_rows()
  
  write_csv(audits, fs::path(out_dir, "TABLE_plane_H_audit.csv"))
  
  if (nrow(audits) > 0) {
    pH <- audits %>%
      mutate(plane = factor(plane, levels = plane)) %>%
      ggplot(aes(x = plane, y = H_med)) +
      geom_point() +
      geom_errorbar(aes(ymin = H_q05, ymax = H_q95), width = 0.15) +
      coord_flip() +
      labs(
        title = "Plane audit: H median with 5–95% interval (feasible RH only)",
        x = NULL, y = "H"
      )
    ggsave(fs::path(out_dir, "PLOT_plane_H_quantiles.png"), pH, width = 8.5, height = 5.2, dpi = 160)
  }
  
  invisible(audits)
}

scenario_domain_flags <- function(df) {
  tibble(
    ever_bad_e     = any(!is.finite(df$e)     | df$e <= 0     | df$e >= 1),
    ever_bad_omega = any(!is.finite(df$omega) | df$omega <= 0 | df$omega >= 1),
    ever_bad_d     = any(!is.finite(df$d)     | df$d <= 0),
    ever_bad_u     = any(!is.finite(df$u)     | df$u <= 0)
  ) %>%
    mutate(domain_ok_all = !(ever_bad_e | ever_bad_omega | ever_bad_d | ever_bad_u))
}

classify_scenario_v2 <- function(df, frac_tail = 0.25,
                                 amp_tol_e = 0.01,
                                 amp_tol_omega = 0.005) {
  
  dom <- scenario_domain_flags(df)
  dft <- tail_window(df, frac_tail)
  eA <- amp(dft$e)
  wA <- amp(dft$omega)
  
  if (!dom$domain_ok_all) return("Blow-up / out-of-domain")
  
  if (is.finite(eA) && is.finite(wA) && (eA >= amp_tol_e || wA >= amp_tol_omega)) {
    return("Cycle / sustained (nontrivial amp)")
  }
  
  "Convergent / damped (tiny amp)"
}

scenario_audit_reclass <- function(out_dir = OUT_DIR) {
  scen_tidy <- read_csv_safe(fs::path(out_dir, "TABLE_scenarios_tidy.csv"))
  if (is.null(scen_tidy)) return(invisible(NULL))
  
  scen_audit <- scen_tidy %>%
    rowwise() %>%
    mutate(
      sim_file = fs::path(out_dir, paste0("scenario_sim_", run, ".csv")),
      sim_exists = fs::file_exists(sim_file),
      regime2 = if (sim_exists) {
        df <- read_csv_safe(sim_file)
        classify_scenario_v2(df)
      } else NA_character_
    ) %>%
    ungroup()
  
  write_csv(scen_audit, fs::path(out_dir, "TABLE_scenarios_audit_reclassified.csv"))
  invisible(scen_audit)
}

## ============================================================
## 10) DRIVER: baseline + planes + figures + scenarios + summaries
## ============================================================

## -----------------------------
## Baseline calibration
## -----------------------------
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
  
  ## utilization adjustment speed (residual follower only)
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
  
  ## optional leakage psi(d)
  psi_zero  = TRUE,
  psi_max   = 0.00, psi_lambda = 5.0, psi_mid = 1.0
)

## -----------------------------
## Candidate planes + ranges
## -----------------------------
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

## ============================================================
## Deliverable 1: rank planes + decision plot (ROBUST)
## ============================================================
scores <- compare_planes(
  par_base,
  planes = planes,
  ranges = ranges,
  n = 80,
  H_band = 1e-3,
  tol = 1e-6,
  calibrate_varsigma0 = TRUE
)

write_csv(scores, fs::path(OUT_DIR, "hopf_plane_scores.csv"))
cat("\n[DRIVER] saved: hopf_plane_scores.csv\n")
print(scores)

plane_outputs <- make_plane_decision_outputs(scores, out_dir = OUT_DIR)

## ============================================================
## Deliverable 2: FIG A (phi1, lambda) sign(H)
## ============================================================
df_A <- hopf_plane(
  par_base,
  p1 = "phi1", p2 = "lambda",
  p1_vals = grid_values(ranges$phi1[1], ranges$phi1[2], 140),
  p2_vals = grid_values(ranges$lambda[1], ranges$lambda[2], 140),
  calibrate_varsigma0 = TRUE,
  H_band = 1e-3
)

write_csv(df_A, fs::path(OUT_DIR, "hopf_plane_phi1_lambda.csv"))

pA <- plot_hopf_plane_signH(df_A, "phi1", "lambda",
                            title = "FIG A: Hopf map (lambda vs phi1) with sign(H)")
ggsave(fs::path(OUT_DIR, "FIG_A_hopf_plane_signH_phi1_lambda.png"),
       pA, width = 7.2, height = 6.0, dpi = 160)

## ============================================================
## Deliverable 3: FIG B (use top robust plane, avoid duplicating A)
## ============================================================
top <- scores %>% slice(1)
p1B <- top$p1
p2B <- top$p2

if (identical(c(p1B, p2B), c("phi1","lambda"))) {
  p1B <- "r0"; p2B <- "lambda"
}

df_B <- hopf_plane(
  par_base,
  p1 = p1B, p2 = p2B,
  p1_vals = grid_values(ranges[[p1B]][1], ranges[[p1B]][2], 140),
  p2_vals = grid_values(ranges[[p2B]][1], ranges[[p2B]][2], 140),
  calibrate_varsigma0 = TRUE,
  H_band = 1e-3
)

csv_B_name <- paste0("hopf_plane_", p1B, "_", p2B, ".csv")
png_B_name <- paste0("FIG_B_hopf_plane_signH_", p1B, "_", p2B, ".png")

write_csv(df_B, fs::path(OUT_DIR, csv_B_name))

pB <- plot_hopf_plane_signH(df_B, p1B, p2B,
                            title = paste0("FIG B: Hopf map (", p2B, " vs ", p1B, ") with sign(H)"))
ggsave(fs::path(OUT_DIR, png_B_name),
       pB, width = 7.2, height = 6.0, dpi = 160)

## ============================================================
## Deliverable 4: scenarios below / near / above Hopf in lambda
## ============================================================
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

cat(sprintf("\n[DRIVER] Hopf approx: lambda* = %.6f at phi1 = %.4f\n", lambda_star, phi1_pick))

lambda_vals <- c(
  below = 0.70 * lambda_star,
  near  = 0.98 * lambda_star,
  above = 1.20 * lambda_star
)

lambda_vals <- pmin(pmax(lambda_vals, ranges$lambda[1]), ranges$lambda[2])

scenario_tbl <- list()

for (nm in names(lambda_vals)) {
  
  par_s <- par_base
  par_s$phi1 <- phi1_pick
  par_s$lambda <- as.numeric(lambda_vals[[nm]])
  
  tag <- paste0("phi1_", sprintf("%.3f", phi1_pick),
                "_lambda_", nm, "_", sprintf("%.2f", par_s$lambda))
  
  df_s <- simulate_full_residual_u(
    par_s,
    t_end = 700, dt = 0.05,
    calibrate_varsigma0 = TRUE
  )
  
  ## export sim CSV
  write_csv(df_s, fs::path(OUT_DIR, paste0("scenario_sim_", tag, ".csv")))
  
  ## export plots
  plot_scenario_timeseries(df_s, tag, out_dir = OUT_DIR)
  plot_scenario_phase_end(df_s, tag, out_dir = OUT_DIR)
  
  ## tail summary
  scenario_tbl[[nm]] <- summarize_tail(df_s, name = tag, frac_tail = 0.25) %>%
    mutate(phi1 = phi1_pick, lambda = par_s$lambda, lambda_star = lambda_star)
}

scenario_tbl <- bind_rows(scenario_tbl)
write_csv(scenario_tbl, fs::path(OUT_DIR, "TABLE_scenarios_tail_diagnostics.csv"))
cat("\n[DRIVER] saved: TABLE_scenarios_tail_diagnostics.csv\n")
print(scenario_tbl)

## ============================================================
## Deliverable 5: scenario tidy + amplitude/period plots
## ============================================================
scen_outputs <- make_scenario_outputs(scenario_tbl, out_dir = OUT_DIR)

## ============================================================
## Deliverable 6: audits (plane H-audit + scenario reclass) + manifest
## ============================================================
audit_plane_files(OUT_DIR)
scenario_audit_reclass(OUT_DIR)

manifest <- make_outputs_manifest(
  OUT_DIR,
  reviewed_figA = TRUE,
  reviewed_figB = TRUE,
  reviewed_Haudit = TRUE
)

write_csv(manifest, fs::path(OUT_DIR, "TABLE_outputs_review_manifest.csv"))

cat("\nDONE. Exports in:\n", OUT_DIR, "\n\nKey files:\n",
    "- hopf_plane_scores.csv\n",
    "- TABLE_plane_ranking_tidy.csv\n",
    "- PLOT_plane_selection_scatter.png\n",
    "- hopf_plane_phi1_lambda.csv\n",
    "- ", csv_B_name, "\n",
    "- FIG_A_hopf_plane_signH_phi1_lambda.png\n",
    "- ", png_B_name, "\n",
    "- TABLE_scenarios_tail_diagnostics.csv\n",
    "- TABLE_scenarios_tidy.csv\n",
    "- PLOT_scenarios_amplitudes_vs_lambdaRel.png\n",
    "- PLOT_scenarios_period_vs_lambdaRel.png\n",
    "- phase_end_* (start/end labeled)\n",
    "- TABLE_plane_H_audit.csv + PLOT_plane_H_quantiles.png\n",
    "- TABLE_scenarios_audit_reclassified.csv\n",
    "- TABLE_outputs_review_manifest.csv\n",
    sep = "")
