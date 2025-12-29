  ############################################################
## wealth_goodwin_hopf_cycles.R  (CONSOLIDATED)
## Goodwin–Minsky core (e, omega, d) + utilization adjustment u(t)
## + wealth inequality diagnostics (chi_A, chi_N, xi)
##
## Includes:
##  - Reduced interior steady state (u = u_star, sigma fixed)
##  - Reduced analytic Jacobian, eigenvalues
##  - RH + Hopf functional H
##  - Robust Hopf search (scan -> bracket -> uniroot, ignores NA)
##  - Finite-diff transversality dH/dmu at mu = mu*
##  - Full simulation with u(t) adjustment and inequality diagnostics
##  - Plots saved in "wealth_goodwin/"
############################################################

## -----------------------------
## 0) Packages
## -----------------------------
suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(fs)
  library(tibble)
})

## -----------------------------
## 1) Helpers
## -----------------------------
inv_logit <- function(x) 1 / (1 + exp(-x))

require_params <- function(par, keys) {
  miss <- setdiff(keys, names(par))
  if (length(miss) > 0) stop("Missing parameters: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

## Phillips (linear)
Phi_fun <- function(e, par) par$phi0 + par$phi1 * e
Phi_e_fun <- function(e, par) par$phi1

pi_fun <- function(omega) 1 - omega

## r_net = (u/sigma)(1 - omega - i d)
r_net_fun <- function(omega, d, u, par) {
  (u / par$sigma) * (pi_fun(omega) - par$i * d)
}

## Investment share kappa(r) (logistic)
kappa_fun <- function(r, par) {
  par$kappa_min + (par$kappa_max - par$kappa_min) /
    (1 + exp(-par$lambda * (r - par$r0)))
}
kappa_r_fun <- function(r, par) {
  s <- inv_logit(par$lambda * (r - par$r0))
  (par$kappa_max - par$kappa_min) * par$lambda * s * (1 - s)
}

## Owners saving propensity s_A(chi_A)
sA_fun <- function(chiA, par) {
  if (isTRUE(par$sA_const)) return(par$sA0)
  s <- inv_logit(par$sA_lambda * (chiA - par$sA_mid))
  par$sA_min + (par$sA_max - par$sA_min) * s
}

## Credit/other leakage psi(d) (optional)
psi_fun <- function(d, par) {
  if (isTRUE(par$psi_zero)) return(0)
  inv_logit(par$psi_lambda * (d - par$psi_mid)) * par$psi_max
}

## Non-owners induced consumption share out of net disposable income share yN = omega - i d
cN_ind_fun <- function(yN, par) {
  if (yN <= 0) return(par$cN_minus)
  par$cN_minus + (par$cN_plus - par$cN_minus) * (yN / (yN + par$cN_eta))
}

## Autonomous baseline consumption share out of capacity output Y^p (trend)
cN_bar_fun <- function(t, par) par$varsigma0 * exp(par$varsigma1 * t)

## Demand wedge Delta (shares of Y): Δ = 1 - leakages - induced comps - owners cons - investment
Delta_fun <- function(t, omega, d, u, par) {
  r <- r_net_fun(omega, d, u, par)
  kap <- kappa_fun(r, par)
  
  rho <- par$sigma / u
  chiA <- par$p * rho + d
  sA <- sA_fun(chiA, par)
  
  yN <- omega - par$i * d
  cN_ind <- cN_ind_fun(yN, par)
  
  psi <- psi_fun(d, par)
  cA <- (1 - sA) * (1 - omega)
  
  1 - psi - cN_ind - cA - kap
}

## Desired utilization: u_d = cbar/Delta
u_desired_fun <- function(t, omega, d, u, par) {
  Delta <- Delta_fun(t, omega, d, u, par)
  if (!is.finite(Delta) || Delta <= 1e-10) return(NA_real_)
  cbar <- cN_bar_fun(t, par)
  cbar / Delta
}

## -----------------------------
## 2) Reduced steady state (u = u_star)
## -----------------------------
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
  ## gY = A*kappa - delta = g_n => kappa* = (g_n + delta)/A
  k_target <- (g_n + par$delta) / A
  
  if (k_target <= par$kappa_min + eps || k_target >= par$kappa_max - eps) {
    return(list(ok = FALSE, reason = "kappa_target outside (kappa_min, kappa_max)"))
  }
  
  s <- (k_target - par$kappa_min) / (par$kappa_max - par$kappa_min)
  if (s <= eps || s >= 1 - eps) return(list(ok = FALSE, reason = "logit inversion s not in (0,1)"))
  
  r_star <- par$r0 + (1 / par$lambda) * log(s / (1 - s))
  
  ## d* = (k_target - (1/A)*r*) / g_n
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

## -----------------------------
## 3) Reduced Jacobian at SS (state: e, omega, d)
## -----------------------------
jacobian_reduced_at_ss <- function(ss, par) {
  e <- ss$e; omega <- ss$omega; d <- ss$d
  g_n <- ss$g_n
  A <- ss$A
  r <- ss$r
  
  kap_r <- kappa_r_fun(r, par)
  Phi_e <- Phi_e_fun(e, par)
  
  J <- matrix(0, 3, 3)
  J[1,2] <- - e * (A^2) * kap_r
  J[1,3] <- - e * (A^2) * par$i * kap_r
  
  J[2,1] <- omega * Phi_e
  
  J[3,2] <- 1 - A * kap_r + d * (A^2) * kap_r
  J[3,3] <- (par$i - g_n) - A * par$i * kap_r + d * (A^2) * par$i * kap_r
  
  list(J = J, kappa_r = kap_r, Phi_e = Phi_e)
}

## -----------------------------
## 4) RH + Hopf H for reduced Jacobian
## -----------------------------
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
  stable_RH <- (a1 > 0) && (a2 > 0) && (a3 > 0) && (a1 * a2 > a3)
  
  list(a1 = a1, a2 = a2, a3 = a3, H = H, stable_RH = stable_RH)
}

analyze_reduced <- function(par) {
  ss <- steady_state_reduced(par)
  if (!isTRUE(ss$ok)) return(ss)
  Jp <- jacobian_reduced_at_ss(ss, par)
  rh <- rh_hopf_reduced(Jp$J)
  eig <- eigen(Jp$J)$values
  
  list(
    ok = TRUE,
    steady_state = ss,
    jacobian = Jp$J,
    eigenvalues = eig,
    rh = rh,
    jac_parts = list(kappa_r = Jp$kappa_r, Phi_e = Jp$Phi_e)
  )
}

hopf_gap_reduced <- function(par) {
  out <- analyze_reduced(par)
  if (!isTRUE(out$ok)) return(NA_real_)
  out$rh$H
}

dH_dmu_fd_reduced <- function(par, mu, h = 1e-4) {
  if (!mu %in% names(par)) stop("Parameter '", mu, "' not in par list.")
  par_p <- par; par_m <- par
  par_p[[mu]] <- par[[mu]] + h
  par_m[[mu]] <- par[[mu]] - h
  Hp <- hopf_gap_reduced(par_p)
  Hm <- hopf_gap_reduced(par_m)
  if (!is.finite(Hp) || !is.finite(Hm)) return(NA_real_)
  (Hp - Hm) / (2*h)
}

## -----------------------------
## 5) Robust Hopf finder: scan -> bracket -> uniroot
## -----------------------------
find_hopf_mu <- function(par, mu, lo, hi,
                         n_grid = 250,
                         tol = 1e-8, maxit = 200,
                         verbose = TRUE) {
  
  f <- function(x) {
    par2 <- par
    par2[[mu]] <- x
    hopf_gap_reduced(par2)
  }
  
  xs <- seq(lo, hi, length.out = n_grid)
  Hs <- vapply(xs, f, numeric(1))
  
  ok <- is.finite(Hs)
  if (!any(ok)) stop("No admissible steady state anywhere on [lo,hi].")
  
  xs_ok <- xs[ok]
  Hs_ok <- Hs[ok]
  
  if (verbose) {
    cat(sprintf("\n[scan] admissible points: %d / %d\n", length(xs_ok), length(xs)))
    cat(sprintf("[scan] H range on admissible set: [%.4g, %.4g]\n", min(Hs_ok), max(Hs_ok)))
  }
  
  sgn <- sign(Hs_ok)
  idx <- which(sgn[-1] * sgn[-length(sgn)] < 0)
  
  if (length(idx) == 0) {
    stop("No sign change in H on admissible subset of [lo,hi]. Try different bracket or parameter.")
  }
  
  a <- xs_ok[idx[1]]
  b <- xs_ok[idx[1] + 1]
  
  if (verbose) cat(sprintf("[scan] bracketing Hopf on [%g, %g]\n", a, b))
  
  uniroot(f, interval = c(a, b), tol = tol, maxiter = maxit)$root
}

## -----------------------------
## 6) Full simulation with utilization adjustment (state: e, omega, d, u)
## -----------------------------
rhs_full <- function(t, state, par) {
  e <- state[["e"]]; w <- state[["omega"]]; d <- state[["d"]]; u <- state[["u"]]
  
  ## guardrails
  e <- min(max(e, 1e-9), 1 - 1e-9)
  w <- min(max(w, 1e-9), 1 - 1e-9)
  d <- max(d, 1e-12)
  u <- max(u, 1e-9)
  
  r <- r_net_fun(w, d, u, par)
  kap <- kappa_fun(r, par)
  
  ## gK = I/K - delta; with K/Y = sigma/u => I/K = (kap Y)/K = kap/(K/Y) = kap/(sigma/u) = u*kap/sigma
  gK <- (u * kap) / par$sigma - par$delta
  
  u_d <- u_desired_fun(t, w, d, u, par)
  if (!is.finite(u_d)) {
    du <- 0
    gU <- 0
  } else {
    du <- par$lambda_u * (u_d - u)
    gU <- du / u
  }
  
  gY <- gK + gU
  g_n <- par$alpha + par$beta
  
  Phi <- Phi_fun(e, par)
  de <- (gY - g_n) * e
  dw <- (Phi - par$alpha) * w
  
  ## debt ratio (financing-gap backbone)
  dd <- kap - (1 - w) + par$i * d - d * gY
  
  ## inequality diagnostics
  rho <- par$sigma / u
  chiA <- par$p * rho + d
  chiN <- -d
  xi <- chiA / abs(chiN)
  
  Delta <- Delta_fun(t, w, d, u, par)
  
  list(
    c(de, dw, dd, du),
    c(r = r, kappa = kap, gK = gK, gU = gU, gY = gY,
      Delta = Delta, u_d = u_d,
      rho = rho, chiA = chiA, chiN = chiN, xi = xi)
  )
}

simulate_full <- function(par, state0 = NULL,
                          t_end = 700, dt = 0.05,
                          perturb = c(e = 0.01, omega = -0.01, d = 0.05, u = -0.05)) {
  
  if (is.null(state0)) {
    ss <- steady_state_reduced(par)
    if (!isTRUE(ss$ok)) stop("No admissible reduced steady state: ", ss$reason)
    
    state0 <- c(
      e = ss$e + perturb[["e"]],
      omega = ss$omega + perturb[["omega"]],
      d = ss$d + perturb[["d"]],
      u = par$u_star + perturb[["u"]]
    )
  }
  
  times <- seq(0, t_end, by = dt)
  
  sol <- ode(
    y = state0,
    times = times,
    func = rhs_full,
    parms = par,
    method = "lsoda"
  )
  
  as.data.frame(sol) |> as_tibble()
}

## -----------------------------
## 7) Plot helpers (save to wealth_goodwin/)
## -----------------------------
plot_states_u <- function(df, tag, out_dir) {
  p <- df |>
    select(time, e, omega, d, u) |>
    pivot_longer(-time, names_to = "var", values_to = "value") |>
    ggplot(aes(x = time, y = value)) +
    geom_line() +
    facet_wrap(~var, scales = "free_y", ncol = 1) +
    labs(title = paste0("States + utilization (", tag, ")"), x = "t", y = NULL)
  ggsave(file.path(out_dir, paste0("states_u_", tag, ".png")), p, width = 7, height = 8, dpi = 160)
  p
}

plot_finance <- function(df, tag, out_dir) {
  p <- df |>
    select(time, r, kappa, gK, gU, gY, Delta, u_d) |>
    pivot_longer(-time, names_to = "var", values_to = "value") |>
    ggplot(aes(x = time, y = value)) +
    geom_line() +
    facet_wrap(~var, scales = "free_y", ncol = 1) +
    labs(title = paste0("Finance + demand wedge (", tag, ")"), x = "t", y = NULL)
  ggsave(file.path(out_dir, paste0("finance_", tag, ".png")), p, width = 7, height = 9, dpi = 160)
  p
}

plot_inequality <- function(df, tag, out_dir) {
  p <- df |>
    select(time, rho, chiA, chiN, xi) |>
    pivot_longer(-time, names_to = "var", values_to = "value") |>
    ggplot(aes(x = time, y = value)) +
    geom_line() +
    facet_wrap(~var, scales = "free_y", ncol = 1) +
    labs(title = paste0("Wealth + inequality diagnostics (", tag, ")"), x = "t", y = NULL)
  ggsave(file.path(out_dir, paste0("inequality_", tag, ".png")), p, width = 7, height = 9, dpi = 160)
  p
}

## -----------------------------
## 8) Reports + near-Hopf runner
## -----------------------------
run_report_reduced <- function(par) {
  out <- analyze_reduced(par)
  if (!isTRUE(out$ok)) {
    cat("FAILED:", out$reason, "\n")
    return(invisible(out))
  }
  ss <- out$steady_state
  rh <- out$rh
  
  cat("\n=== REDUCED STEADY STATE (u = u_star) ===\n")
  cat(sprintf("e* = %.4f | omega* = %.4f | d* = %.4f | u* = %.4f\n",
              ss$e, ss$omega, ss$d, par$u_star))
  cat(sprintf("r* = %.6f | kappa* = %.6f | g_n = %.4f | A=u*/sigma = %.6f\n",
              ss$r, ss$kappa, ss$g_n, ss$A))
  
  cat("\n=== REDUCED EIGENVALUES ===\n")
  print(out$eigenvalues)
  
  cat("\n=== REDUCED RH + HOPF ===\n")
  cat(sprintf("a1=%.6g a2=%.6g a3=%.6g | H=%.6g | stable_RH=%s\n",
              rh$a1, rh$a2, rh$a3, rh$H, rh$stable_RH))
  
  invisible(out)
}

run_near_hopf <- function(par, mu = "sigma", lo = 1.5, hi = 6.0,
                          eps = 0.03, t_end = 700, dt = 0.05,
                          out_dir = "wealth_goodwin",
                          n_grid = 250) {
  
  dir_create(out_dir)
  
  mu_star <- find_hopf_mu(par, mu, lo, hi, n_grid = n_grid, verbose = TRUE)
  cat(sprintf("\nHopf approx at %s* = %.8f\n", mu, mu_star))
  
  ## transversality at Hopf (finite diff at mu*)
  par_star <- par
  par_star[[mu]] <- mu_star
  dH <- dH_dmu_fd_reduced(par_star, mu, h = abs(mu_star) * 1e-4 + 1e-8)
  cat(sprintf("Finite-diff transversality: dH/d%s |_{%s*} ≈ %.6g\n", mu, mu, dH))
  
  ## below/above
  par_lo <- par_star; par_hi <- par_star
  par_lo[[mu]] <- mu_star * (1 - eps)
  par_hi[[mu]] <- mu_star * (1 + eps)
  
  cat("\n--- Reduced report: below ---\n")
  run_report_reduced(par_lo)
  cat("\n--- Reduced report: above ---\n")
  run_report_reduced(par_hi)
  
  ## full sims
  tag_lo <- paste0(mu, "_belowHopf")
  tag_hi <- paste0(mu, "_aboveHopf")
  
  df_lo <- simulate_full(par_lo, t_end = t_end, dt = dt)
  df_hi <- simulate_full(par_hi, t_end = t_end, dt = dt)
  
  plot_states_u(df_lo, tag_lo, out_dir)
  plot_finance(df_lo, tag_lo, out_dir)
  plot_inequality(df_lo, tag_lo, out_dir)
  
  plot_states_u(df_hi, tag_hi, out_dir)
  plot_finance(df_hi, tag_hi, out_dir)
  plot_inequality(df_hi, tag_hi, out_dir)
  
  invisible(list(mu_star = mu_star, dH = dH, below = df_lo, above = df_hi))
}

############################################################
## 9) Example run
############################################################
par <- list(
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
  
  ## utilization adjustment speed
  lambda_u = 2.0,
  
  ## price (wealth accounting)
  p = 1.0,
  
  ## autonomous non-owners consumption share out of capacity output
  varsigma0 = 0.55,
  varsigma1 = 0.00,  ## keep 0 for Hopf neighborhood experiments
  
  ## induced non-owner consumption
  cN_minus = 0.02,
  cN_plus  = 0.20,
  cN_eta   = 0.10,
  
  ## owners saving
  sA_const  = TRUE,
  sA0       = 0.40,
  sA_min    = 0.20, sA_max = 0.80, sA_lambda = 2.0, sA_mid = 2.0,
  
  ## optional leakage / credit term psi(d)
  psi_zero  = TRUE,
  psi_max   = 0.00, psi_lambda = 5.0, psi_mid = 1.0
)

out_dir <- "wealth_goodwin"
dir_create(out_dir)

## baseline reduced diagnostics
run_report_reduced(par)

## near-Hopf: vary sigma (ensure bracket feasible)
res <- run_near_hopf(
  par,
  mu = "lambda",
  lo = 10, hi = 150,
  eps = 0.03,
  t_end = 700, dt = 0.05,
  out_dir = out_dir,
  n_grid = 400
)


## ============================================================
## Regime classification utilities for res$below / res$above
## ============================================================

tail_window <- function(df, frac = 0.35) {
  n <- nrow(df)
  i0 <- max(1, floor((1 - frac) * n))
  df[i0:n, , drop = FALSE]
}

lin_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(NA_real_)
  coef(lm(y[ok] ~ x[ok]))[2]
}

amp <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NA_real_)
  quantile(x, 0.95) - quantile(x, 0.05)
}

cv <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NA_real_)
  sd(x) / (abs(mean(x)) + 1e-12)
}

## crude cycle detector:
## - amplitude in tail > small threshold
## - AND amplitude is not exploding (tail amp not >> earlier amp)
cycle_signature <- function(df, var = "e", frac_tail = 0.35, explode_factor = 2.0) {
  df_tail <- tail_window(df, frac_tail)
  df_early <- df[1:floor(nrow(df)*(1-frac_tail)), , drop=FALSE]
  
  a_tail <- amp(df_tail[[var]])
  a_early <- amp(df_early[[var]])
  
  if (!is.finite(a_tail) || !is.finite(a_early)) return(list(cycle = NA, explode = NA, a_tail = a_tail))
  
  explode <- (a_tail > explode_factor * a_early) && (a_tail > 1e-3)
  cycle   <- (a_tail > 1e-3) && !explode
  
  list(cycle = cycle, explode = explode, a_tail = a_tail, a_early = a_early)
}

summarize_run <- function(df, name = "run", frac_tail = 0.35) {
  df_tail <- tail_window(df, frac_tail)
  
  ## slopes (tail)
  s_logxi <- lin_slope(df_tail$time, log(pmax(df_tail$xi, 1e-12)))
  s_d     <- lin_slope(df_tail$time, df_tail$d)
  s_u     <- lin_slope(df_tail$time, df_tail$u)
  s_e     <- lin_slope(df_tail$time, df_tail$e)
  
  ## amplitudes (tail)
  a_e <- amp(df_tail$e)
  a_w <- amp(df_tail$omega)
  a_d <- amp(df_tail$d)
  a_u <- amp(df_tail$u)
  a_x <- amp(log(pmax(df_tail$xi, 1e-12)))
  
  ## growth diagnostics (tail means)
  gY_bar <- mean(df_tail$gY, na.rm = TRUE)
  gK_bar <- mean(df_tail$gK, na.rm = TRUE)
  gU_bar <- mean(df_tail$gU, na.rm = TRUE)
  
  ## cycle/explosion tests on a couple of key vars
  cyc_e <- cycle_signature(df, "e", frac_tail)
  cyc_d <- cycle_signature(df, "d", frac_tail)
  
  ## boundary hits (proxy for “crisis corridor”)
  hit_bound <- any(df_tail$e < 0.02 | df_tail$e > 0.98 |
                     df_tail$omega < 0.02 | df_tail$omega > 0.98 |
                     df_tail$d > quantile(df$d, 0.999, na.rm = TRUE) * 0.999, na.rm = TRUE)
  
  tibble::tibble(
    run = name,
    time_end = max(df$time, na.rm=TRUE),
    
    ## tail slopes
    slope_log_xi = s_logxi,
    slope_d      = s_d,
    slope_u      = s_u,
    slope_e      = s_e,
    
    ## tail amplitudes
    amp_e = a_e, amp_omega = a_w, amp_d = a_d, amp_u = a_u, amp_log_xi = a_x,
    
    ## tail means
    gY_bar = gY_bar, gK_bar = gK_bar, gU_bar = gU_bar,
    
    ## cycle/explosion flags (heuristic)
    cycle_e = cyc_e$cycle, explode_e = cyc_e$explode,
    cycle_d = cyc_d$cycle, explode_d = cyc_d$explode,
    
    ## boundary
    hit_bound = hit_bound
  )
}

classify_macro <- function(row, eps_amp = 5e-3, eps_slope = 1e-4) {
  ## Crisis corridor if exploding or boundary hits
  if (isTRUE(row$hit_bound) || isTRUE(row$explode_e) || isTRUE(row$explode_d)) {
    return("Crisis corridor / runaway")
  }
  
  ## Stationary equilibrium if amplitudes tiny AND slopes tiny
  if (is.finite(row$amp_e) && row$amp_e < eps_amp &&
      is.finite(row$amp_d) && row$amp_d < eps_amp &&
      abs(row$slope_e %||% 0) < eps_slope &&
      abs(row$slope_d %||% 0) < eps_slope) {
    return("Stable equilibrium (damped)")
  }
  
  ## Cyclical if sustained amplitude without explosion
  if (isTRUE(row$cycle_e) || isTRUE(row$cycle_d)) {
    return("Endogenous cycle (limit-cycle candidate)")
  }
  
  ## Otherwise: slow drift / quasi-cycle
  "Quasi-cyclical / drifting"
}

classify_inequality <- function(row, eps = 1e-5) {
  ## Rising if tail drift of log xi is positive
  if (!is.finite(row$slope_log_xi)) return("Unknown")
  if (row$slope_log_xi > eps) "Rising inequality" else "Contained inequality"
}

`%||%` <- function(a, b) if (is.null(a) || !is.finite(a)) b else a

regime_map <- function(sumdf) {
  sumdf %>%
    rowwise() %>%
    mutate(
      macro_regime = classify_macro(cur_data()),
      ineq_regime  = classify_inequality(cur_data()),
      regime_2x2   = paste(macro_regime, "×", ineq_regime)
    ) %>%
    ungroup()
}

## ------------------------------------------------------------
## Phase plots (saved)
## ------------------------------------------------------------
phase_plot <- function(df, x, y, tag, out_dir) {
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_path(alpha = 0.8) +
    labs(title = paste0("Phase: ", y, " vs ", x, " (", tag, ")"), x = x, y = y)
  ggsave(file.path(out_dir, paste0("phase_", y, "_vs_", x, "_", tag, ".png")),
         p, width = 6.5, height = 5.5, dpi = 160)
  p
}

## ============================================================
## APPLY TO YOUR RESULTS
## ============================================================

out_dir <- "wealth_goodwin"
dir_create(out_dir)

## assumes you already have: res <- run_near_hopf(... mu="lambda" ...)
s_below <- summarize_run(res$below, name = "lambda_belowHopf")
s_above <- summarize_run(res$above, name = "lambda_aboveHopf")

summary_tbl <- bind_rows(s_below, s_above) %>%
  regime_map()

print(summary_tbl)

write.csv(summary_tbl, file.path(out_dir, "regime_summary.csv"), row.names = FALSE)

## Phase plots (tail only is often cleaner, but keep full first)
phase_plot(res$below, "e", "omega", "lambda_belowHopf", out_dir)
phase_plot(res$below, "omega", "d", "lambda_belowHopf", out_dir)
phase_plot(res$below, "d", "u", "lambda_belowHopf", out_dir)

phase_plot(res$above, "e", "omega", "lambda_aboveHopf", out_dir)
phase_plot(res$above, "omega", "d", "lambda_aboveHopf", out_dir)
phase_plot(res$above, "d", "u", "lambda_aboveHopf", out_dir)




## ============================================================
## Plotting upgrades: show the END / FINAL trajectory clearly
## ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

tail_window <- function(df, frac = 0.25) {
  n <- nrow(df)
  i0 <- max(1, floor((1 - frac) * n))
  df[i0:n, , drop = FALSE]
}

add_tail_flag <- function(df, frac = 0.25) {
  n <- nrow(df)
  i0 <- max(1, floor((1 - frac) * n))
  df %>%
    mutate(
      idx = row_number(),
      is_tail = idx >= i0,
      t_tail_start = time[i0]
    )
}

## ---------- 1) Time-series: full path + tail emphasized + end marker ----------
plot_timeseries_end <- function(df, vars = c("e","omega","d","u","xi"),
                                tag, out_dir = "wealth_goodwin",
                                tail_frac = 0.25) {
  
  df2 <- add_tail_flag(df, tail_frac)
  
  df_long <- df2 %>%
    select(time, all_of(vars), is_tail, t_tail_start) %>%
    pivot_longer(cols = all_of(vars), names_to = "var", values_to = "value")
  
  # end points (per variable)
  df_end <- df_long %>%
    group_by(var) %>%
    slice_tail(n = 1) %>%
    ungroup()
  
  p <- ggplot(df_long, aes(x = time, y = value, group = var)) +
    geom_line(aes(alpha = is_tail)) +
    geom_vline(xintercept = unique(df_long$t_tail_start), linetype = "dashed") +
    geom_point(data = df_end, size = 2) +
    geom_text(
      data = df_end,
      aes(label = sprintf("end=%.3g", value)),
      hjust = -0.05, vjust = 0.5,
      check_overlap = TRUE
    ) +
    facet_wrap(~var, scales = "free_y", ncol = 1) +
    scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1.0), guide = "none") +
    labs(
      title = paste0("Time paths with tail + end point (", tag, ")"),
      subtitle = paste0("Tail window = last ", round(100*tail_frac), "%; dashed line = tail start"),
      x = "t", y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(5.5, 35, 5.5, 5.5))
  
  ggsave(file.path(out_dir, paste0("timeseries_end_", tag, ".png")),
         p, width = 7.5, height = 9, dpi = 160)
  
  p
}

## ---------- 2) Tail-only zoom time-series (to see the attractor clearly) ----------
plot_tail_zoom <- function(df, vars = c("e","omega","d","u","xi"),
                           tag, out_dir = "wealth_goodwin",
                           tail_frac = 0.25) {
  
  df_tail <- tail_window(df, tail_frac) %>%
    select(time, all_of(vars)) %>%
    pivot_longer(-time, names_to = "var", values_to = "value")
  
  p <- ggplot(df_tail, aes(x = time, y = value)) +
    geom_line() +
    facet_wrap(~var, scales = "free_y", ncol = 1) +
    labs(
      title = paste0("Tail zoom (final dynamics) (", tag, ")"),
      subtitle = paste0("Only last ", round(100*tail_frac), "% of time"),
      x = "t", y = NULL
    )
  
  ggsave(file.path(out_dir, paste0("tail_zoom_", tag, ".png")),
         p, width = 7, height = 9, dpi = 160)
  
  p
}

## ---------- 3) Phase plot: color by time + start/end + direction arrows ----------
df <- df %>% arrange(time)
stopifnot(all(diff(df$time) > 0))
phase_plot_end <- function(df, x, y, tag, out_dir = "wealth_goodwin",
                           tail_frac = 0.25, arrow_every = 150) {
  
  df2 <- add_tail_flag(df, tail_frac)
  
  # start/end markers
  df_start <- df2 %>% slice_head(n = 1)
  df_end   <- df2 %>% slice_tail(n = 1)
  
  # arrows along the *tail* (so “final trajectory” is explicit)
  df_tail <- df2 %>% filter(is_tail)
  n_tail <- nrow(df_tail)
  idx <- seq(1, n_tail - 1, by = arrow_every)
  segs <- tibble(
    x0 = df_tail[[x]][idx],
    y0 = df_tail[[y]][idx],
    x1 = df_tail[[x]][idx + 1],
    y1 = df_tail[[y]][idx + 1]
  )
  
  p <- ggplot(df2, aes(x = .data[[x]], y = .data[[y]])) +
    geom_path(aes(alpha = is_tail), linewidth = 0.7) +
    geom_segment(
      data = segs,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      arrow = arrow(type = "closed", length = unit(0.10, "inches")),
      inherit.aes = FALSE
    ) +
    geom_point(data = df_start, size = 2, inherit.aes = FALSE,
               aes(x = .data[[x]], y = .data[[y]])) +
    geom_text(data = df_start, inherit.aes = FALSE,
              aes(x = .data[[x]], y = .data[[y]], label = "start"),
              hjust = -0.1, vjust = 1.1, check_overlap = TRUE) +
    geom_point(data = df_end, size = 2, inherit.aes = FALSE,
               aes(x = .data[[x]], y = .data[[y]])) +
    geom_text(data = df_end, inherit.aes = FALSE,
              aes(x = .data[[x]], y = .data[[y]], label = "end"),
              hjust = -0.1, vjust = -0.6, check_overlap = TRUE) +
    scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1.0), guide = "none") +
    labs(
      title = paste0("Phase plot with final-trajectory arrows (", tag, ")"),
      subtitle = paste0("Arrows drawn on tail window (last ", round(100*tail_frac), "%)"),
      x = x, y = y
    )
  
  ggsave(file.path(out_dir, paste0("phase_end_", y, "_vs_", x, "_", tag, ".png")),
         p, width = 6.5, height = 5.5, dpi = 160)
  
  p
}

## ---------- 4) Convenience wrapper: make the “final trajectory pack” ----------
final_trajectory_pack <- function(df, tag, out_dir = "wealth_goodwin", tail_frac = 0.25) {
  dir.create(out_dir, showWarnings = FALSE)
  
  plot_timeseries_end(df, vars = c("e","omega","d","u","xi"), tag = tag, out_dir = out_dir, tail_frac = tail_frac)
  plot_tail_zoom(df, vars = c("e","omega","d","u","xi"), tag = tag, out_dir = out_dir, tail_frac = tail_frac)
  
  phase_plot_end(df, "e", "omega", tag = tag, out_dir = out_dir, tail_frac = tail_frac)
  phase_plot_end(df, "omega", "d", tag = tag, out_dir = out_dir, tail_frac = tail_frac)
  phase_plot_end(df, "d", "u", tag = tag, out_dir = out_dir, tail_frac = tail_frac)
  
  # inequality-focused phase (optional, but usually useful)
  if ("xi" %in% names(df)) {
    phase_plot_end(df, "d", "xi", tag = tag, out_dir = out_dir, tail_frac = tail_frac)
    phase_plot_end(df, "u", "xi", tag = tag, out_dir = out_dir, tail_frac = tail_frac)
  }
  
  invisible(TRUE)
}


library(ggplot2)
library(viridis)   # for scale_color_viridis_c()

ggplot(df, aes(x = e, y = omega, color = time)) +
  geom_path() +
  scale_color_viridis_c()

## ============================================================
## Plotting upgrades: show the END / FINAL trajectory clearly
##   - Time series: tail emphasized + end labels
##   - Tail zoom: last X% only
##   - Phase plots: tail arrows + start/end markers
##   - Optional: phase plot colored by time (viridis)
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
  library(fs)
  # optional, only used by phase_plot_time_color()
  # install.packages("viridis") if needed
  # library(viridis)
})

## ---------- utilities ----------
tail_window <- function(df, frac = 0.25) {
  df <- dplyr::arrange(df, time)
  n <- nrow(df)
  i0 <- max(1, floor((1 - frac) * n))
  df[i0:n, , drop = FALSE]
}

add_tail_flag <- function(df, frac = 0.25) {
  df <- dplyr::arrange(df, time)
  stopifnot(all(diff(df$time) > 0))
  
  n <- nrow(df)
  i0 <- max(1, floor((1 - frac) * n))
  
  df %>%
    dplyr::mutate(
      idx = dplyr::row_number(),
      is_tail = idx >= i0,
      t_tail_start = time[i0]
    )
}

## ---------- 1) Time-series: full path + tail emphasized + end marker ----------
plot_timeseries_end <- function(df,
                                vars = c("e","omega","d","u","xi"),
                                tag,
                                out_dir = "wealth_goodwin",
                                tail_frac = 0.25) {
  
  df2 <- add_tail_flag(df, tail_frac)
  
  vars_ok <- intersect(vars, names(df2))
  stopifnot(length(vars_ok) > 0)
  
  df_long <- df2 %>%
    dplyr::select(time, dplyr::all_of(vars_ok), is_tail, t_tail_start) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(vars_ok), names_to = "var", values_to = "value")
  
  df_end <- df_long %>%
    dplyr::group_by(var) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup()
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = time, y = value)) +
    ggplot2::geom_line(ggplot2::aes(alpha = is_tail)) +
    ggplot2::geom_vline(xintercept = unique(df_long$t_tail_start), linetype = "dashed") +
    ggplot2::geom_point(data = df_end, size = 2) +
    ggplot2::geom_text(
      data = df_end,
      ggplot2::aes(label = sprintf("end=%.3g", value)),
      hjust = -0.05, vjust = 0.5,
      check_overlap = TRUE
    ) +
    ggplot2::facet_wrap(~var, scales = "free_y", ncol = 1) +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1.0), guide = "none") +
    ggplot2::labs(
      title = paste0("Time paths with tail + end point (", tag, ")"),
      subtitle = paste0("Tail window = last ", round(100*tail_frac), "%; dashed line = tail start"),
      x = "t", y = NULL
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(5.5, 35, 5.5, 5.5))
  
  fs::dir_create(out_dir)
  ggplot2::ggsave(file.path(out_dir, paste0("timeseries_end_", tag, ".png")),
                  p, width = 7.5, height = 9, dpi = 160)
  
  p
}

## ---------- 2) Tail-only zoom time-series ----------
plot_tail_zoom <- function(df,
                           vars = c("e","omega","d","u","xi"),
                           tag,
                           out_dir = "wealth_goodwin",
                           tail_frac = 0.25) {
  
  df_tail <- tail_window(df, tail_frac)
  
  vars_ok <- intersect(vars, names(df_tail))
  stopifnot(length(vars_ok) > 0)
  
  df_long <- df_tail %>%
    dplyr::select(time, dplyr::all_of(vars_ok)) %>%
    tidyr::pivot_longer(-time, names_to = "var", values_to = "value")
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = time, y = value)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~var, scales = "free_y", ncol = 1) +
    ggplot2::labs(
      title = paste0("Tail zoom (final dynamics) (", tag, ")"),
      subtitle = paste0("Only last ", round(100*tail_frac), "% of time"),
      x = "t", y = NULL
    )
  
  fs::dir_create(out_dir)
  ggplot2::ggsave(file.path(out_dir, paste0("tail_zoom_", tag, ".png")),
                  p, width = 7, height = 9, dpi = 160)
  
  p
}

## ---------- 3) Phase plot: tail arrows + start/end ----------
phase_end_plot <- function(dat, x, y, tag,
                           out_dir = "wealth_goodwin",
                           tail_frac = 0.25,
                           arrow_every = 40) {
  
  dat <- dplyr::arrange(dat, time)
  
  n <- nrow(dat)
  i0 <- max(1, floor((1 - tail_frac) * n))
  dat <- dat %>%
    dplyr::mutate(idx = dplyr::row_number(),
                  is_tail = idx >= i0)
  
  # start/end points
  p_start <- dat[1, ]
  p_end   <- dat[n, ]
  
  # arrows along tail (subsample)
  tail_dat <- dat %>% dplyr::filter(is_tail)
  tail_dat <- tail_dat %>% dplyr::mutate(
    x0 = .data[[x]],
    y0 = .data[[y]],
    x1 = dplyr::lead(.data[[x]]),
    y1 = dplyr::lead(.data[[y]])
  )
  arrows <- tail_dat %>%
    dplyr::filter(!is.na(x1), idx %% arrow_every == 0)
  
  p <- ggplot(dat, aes(x = .data[[x]], y = .data[[y]])) +
    geom_path(aes(alpha = is_tail)) +
    scale_alpha_manual(values = c(`FALSE` = 0.25, `TRUE` = 1.0), guide = "none") +
    geom_point(data = p_start, size = 2) +
    geom_point(data = p_end, size = 2) +
    geom_segment(
      data = arrows,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      arrow = arrow(length = unit(0.12, "inches")),
      inherit.aes = FALSE
    ) +
    labs(
      title = paste0("Phase (tail emphasized + direction): ", y, " vs ", x, " (", tag, ")"),
      subtitle = paste0("Tail = last ", round(100*tail_frac), "% | points mark start/end"),
      x = x, y = y
    )
  
  ggsave(file.path(out_dir, paste0("phase_end_", y, "_vs_", x, "_", tag, ".png")),
         p, width = 6.5, height = 5.5, dpi = 160)
  
  p
}

## ---------- Optional: phase plot colored by time (needs viridis) ----------
phase_plot_time_color <- function(df, x, y, tag,
                                  out_dir = "wealth_goodwin") {
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' not installed. Install it or skip phase_plot_time_color().")
  }
  
  df <- dplyr::arrange(df, time)
  stopifnot(all(diff(df$time) > 0))
  stopifnot(x %in% names(df), y %in% names(df))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = time)) +
    ggplot2::geom_path(linewidth = 0.7) +
    viridis::scale_color_viridis_c() +
    ggplot2::labs(
      title = paste0("Phase plot colored by time (", tag, ")"),
      x = x, y = y, color = "time"
    )
  
  fs::dir_create(out_dir)
  ggplot2::ggsave(file.path(out_dir, paste0("phase_timecolor_", y, "_vs_", x, "_", tag, ".png")),
                  p, width = 6.5, height = 5.5, dpi = 160)
  
  p
}

## ---------- 4) Convenience wrapper: “final trajectory pack” ----------
final_trajectory_pack <- function(df, tag,
                                  out_dir = "wealth_goodwin",
                                  tail_frac = 0.25,
                                  arrow_every = 150) {
  
  fs::dir_create(out_dir)
  
  plot_timeseries_end(df, vars = c("e","omega","d","u","xi"),
                      tag = tag, out_dir = out_dir, tail_frac = tail_frac)
  
  plot_tail_zoom(df, vars = c("e","omega","d","u","xi"),
                 tag = tag, out_dir = out_dir, tail_frac = tail_frac)
  
  phase_plot_end(df, "e",     "omega", tag = tag, out_dir = out_dir,
                 tail_frac = tail_frac, arrow_every = arrow_every)
  
  phase_plot_end(df, "omega", "d",     tag = tag, out_dir = out_dir,
                 tail_frac = tail_frac, arrow_every = arrow_every)
  
  if ("u" %in% names(df)) {
    phase_plot_end(df, "d", "u", tag = tag, out_dir = out_dir,
                   tail_frac = tail_frac, arrow_every = arrow_every)
  }
  
  if ("xi" %in% names(df)) {
    phase_plot_end(df, "d", "xi", tag = tag, out_dir = out_dir,
                   tail_frac = tail_frac, arrow_every = arrow_every)
    if ("u" %in% names(df)) {
      phase_plot_end(df, "u", "xi", tag = tag, out_dir = out_dir,
                     tail_frac = tail_frac, arrow_every = arrow_every)
    }
  }
  
  invisible(TRUE)
}


final_trajectory_pack(res$below, tag = "lambda_belowHopf", out_dir = "wealth_goodwin", tail_frac = 0.25)
final_trajectory_pack(res$above, tag = "lambda_aboveHopf", out_dir = "wealth_goodwin", tail_frac = 0.25)




phase_end_plot(res$below, "e", "omega", "lambda_belowHopf", out_dir)
phase_end_plot(res$above, "e", "omega", "lambda_aboveHopf", out_dir)

phase_end_plot(res$below, "omega", "d", "lambda_belowHopf", out_dir)
phase_end_plot(res$above, "omega", "d", "lambda_aboveHopf", out_dir)

phase_end_plot(res$below, "d", "u", "lambda_belowHopf", out_dir)
phase_end_plot(res$above, "d", "u", "lambda_aboveHopf", out_dir)

