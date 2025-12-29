############################################################
## keen_goodwin_minsky_3d_slowfinance.R
## Keen-orthodox Goodwin–Minsky 3D with Buffett discipline
## + Slow-finance option (time-scale separation via tau_d)
##
## Key idea:
##   d-dot = (1/tau_d) * [kappa - pi + i d - d gY]
## with tau_d >> 1 => debt moves slowly (long wave),
## while (e, omega) keep distributive cycles (short wave).
############################################################

suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scatterplot3d)
  library(fs)
  library(tibble)
})

## -----------------------------
## 0) Helpers + core blocks
## -----------------------------
inv_logit <- function(x) 1 / (1 + exp(-x))

require_params <- function(par, keys) {
  miss <- setdiff(keys, names(par))
  if (length(miss) > 0) stop("Missing parameters: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

## Buffett switch
S_fun <- function(d, par) inv_logit(par$psi * (d / par$dbar - 1))
S_d_fun <- function(d, par) {
  S <- S_fun(d, par)
  (par$psi / par$dbar) * S * (1 - S)
}

## Shares + profit rate
pi_fun <- function(omega) 1 - omega
r_fun  <- function(omega, d, par) (pi_fun(omega) - par$i * d) / par$sigma

## Investment + derivative
kappa_fun <- function(r, par) {
  par$kappa_min + (par$kappa_max - par$kappa_min) /
    (1 + exp(-par$lambda * (r - par$r0)))
}
kappa_r_fun <- function(r, par) {
  s <- inv_logit(par$lambda * (r - par$r0))
  (par$kappa_max - par$kappa_min) * par$lambda * s * (1 - s)
}

## Output growth (capacity branch closure)
gY_fun <- function(r, par) kappa_fun(r, par) / par$sigma - par$delta

## Phillips (disciplined)
Phi_fun <- function(e, d, par) {
  par$phi0 + par$phi1 * e - par$eta * S_fun(d, par)
}

## -----------------------------
## 1) Interior steady state (closed form)
## (Note: steady state does NOT depend on tau_d: time scale only)
## -----------------------------
steady_state <- function(par, strict = FALSE, e_band = c(0.75, 0.95), eps = 1e-12) {
  require_params(par, c(
    "alpha","beta","delta","sigma","i",
    "kappa_min","kappa_max","lambda","r0",
    "phi0","phi1","eta","dbar","psi"
  ))
  
  g_n <- par$alpha + par$beta
  if (g_n <= 0) return(list(ok = FALSE, reason = "g_n <= 0 (alpha+beta must be > 0)"))
  
  k_target <- par$sigma * (g_n + par$delta)
  
  if (k_target <= par$kappa_min + eps || k_target >= par$kappa_max - eps)
    return(list(ok = FALSE, reason = "kappa_target outside (kappa_min, kappa_max)"))
  
  s <- (k_target - par$kappa_min) / (par$kappa_max - par$kappa_min)
  if (s <= eps || s >= 1 - eps) return(list(ok = FALSE, reason = "s not in (0,1)"))
  
  r_star <- par$r0 + (1 / par$lambda) * log(s / (1 - s))
  d_star <- (k_target - par$sigma * r_star) / g_n
  if (d_star <= eps) return(list(ok = FALSE, reason = "d* <= 0 (no interior debt steady state)"))
  
  pi_star <- k_target + d_star * (par$i - g_n)
  omega_star <- 1 - pi_star
  if (omega_star <= eps || omega_star >= 1 - eps)
    return(list(ok = FALSE, reason = "omega* not in (0,1)"))
  
  if (par$phi1 <= 0) return(list(ok = FALSE, reason = "phi1 must be > 0"))
  e_star <- (par$alpha - par$phi0 + par$eta * S_fun(d_star, par)) / par$phi1
  if (e_star <= eps || e_star >= 1 - eps)
    return(list(ok = FALSE, reason = "e* not in (0,1)"))
  
  if (strict && (e_star < e_band[1] || e_star > e_band[2]))
    return(list(ok = FALSE, reason = "e* outside plausibility band"))
  
  list(ok = TRUE, e = e_star, omega = omega_star, d = d_star, r = r_star,
       kappa = k_target, g_n = g_n)
}

## -----------------------------
## 2) Analytic Jacobian at steady state
## -----------------------------
jacobian_at_ss <- function(ss, par) {
  e <- ss$e; omega <- ss$omega; d <- ss$d
  g_n <- ss$g_n; r <- ss$r
  
  kap_r <- kappa_r_fun(r, par)
  Sd <- S_d_fun(d, par)
  
  A <- - e * kap_r / (par$sigma^2)
  B <- - e * par$i * kap_r / (par$sigma^2)
  C <- omega * par$phi1
  D <- - omega * par$eta * Sd
  E <- 1 - (1 - d / par$sigma) * (kap_r / par$sigma)
  F <- (par$i - g_n) - (1 - d / par$sigma) * (par$i * kap_r / par$sigma)
  
  J <- matrix(c(
    0, A, B,
    C, 0, D,
    0, E, F
  ), nrow = 3, byrow = TRUE)
  
  list(J = J, A = A, B = B, C = C, D = D, E = E, F = F,
       kappa_r = kap_r, S_d = Sd)
}

rh_hopf <- function(Jparts) {
  A <- Jparts$A; B <- Jparts$B; C <- Jparts$C
  D <- Jparts$D; E <- Jparts$E; F <- Jparts$F
  
  a1 <- -F
  a2 <- -A * C - E * D
  a3 <- C * (A * F - B * E)
  
  H <- a1 * a2 - a3
  stable_RH <- (a1 > 0) && (a2 > 0) && (a3 > 0) && (a1 * a2 > a3)
  
  list(a1 = a1, a2 = a2, a3 = a3, H = H, stable_RH = stable_RH)
}

analyze_model <- function(par, strict = FALSE, e_band = c(0.75, 0.95)) {
  ss <- steady_state(par, strict = strict, e_band = e_band)
  if (!isTRUE(ss$ok)) return(ss)
  
  Jp <- jacobian_at_ss(ss, par)
  rh <- rh_hopf(Jp)
  eig <- eigen(Jp$J)$values
  
  list(ok = TRUE, steady_state = ss, jacobian = Jp$J,
       jac_parts = Jp[c("A","B","C","D","E","F","kappa_r","S_d")],
       eigenvalues = eig, rh = rh)
}

hopf_gap <- function(par) {
  out <- analyze_model(par, strict = FALSE)
  if (!isTRUE(out$ok)) return(NA_real_)
  out$rh$H
}

dH_dmu_fd <- function(par, mu, h = 1e-3) {
  if (!mu %in% names(par)) stop("mu not found in par: ", mu)
  par_p <- par; par_m <- par
  par_p[[mu]] <- par[[mu]] + h
  par_m[[mu]] <- par[[mu]] - h
  Hp <- hopf_gap(par_p)
  Hm <- hopf_gap(par_m)
  if (is.na(Hp) || is.na(Hm)) return(NA_real_)
  (Hp - Hm) / (2*h)
}

run_report <- function(par) {
  out <- analyze_model(par, strict = FALSE)
  if (!isTRUE(out$ok)) {
    cat("FAILED:", out$reason, "\n")
    return(invisible(out))
  }
  ss <- out$steady_state
  rh <- out$rh
  jp <- out$jac_parts
  
  cat("\n=== STEADY STATE ===\n")
  cat(sprintf("e* = %.4f | omega* = %.4f | d* = %.4f | r* = %.5f\n",
              ss$e, ss$omega, ss$d, ss$r))
  cat(sprintf("kappa_target = %.4f | g_n = %.4f\n", ss$kappa, ss$g_n))
  
  cat("\n=== JACOBIAN PARTS ===\n")
  cat(sprintf("A=%.6f B=%.6f C=%.6f D=%.6f E=%.6f F=%.6f\n",
              jp$A, jp$B, jp$C, jp$D, jp$E, jp$F))
  cat(sprintf("kappa_r=%.6f | S_d=%.6f\n", jp$kappa_r, jp$S_d))
  
  cat("\n=== EIGENVALUES ===\n")
  print(out$eigenvalues)
  
  cat("\n=== ROUTH–HURWITZ ===\n")
  cat(sprintf("a1=%.6g a2=%.6g a3=%.6g | H=%.6g | stable_RH=%s\n",
              rh$a1, rh$a2, rh$a3, rh$H, rh$stable_RH))
  
  invisible(out)
}

## -----------------------------
## 3) ODE system (with slow-finance option)
## -----------------------------
rhs_keengm_slowfinance <- function(t, state, par) {
  e <- state[["e"]]
  w <- state[["omega"]]
  d <- state[["d"]]
  
  ## numerical guards (not “resets”)
  d_floor <- if (!is.null(par$d_floor)) par$d_floor else 1e-9
  e_safe  <- min(max(e, 1e-9), 1 - 1e-9)
  w_safe  <- min(max(w, 1e-9), 1 - 1e-9)
  d_safe  <- max(d, d_floor)
  
  ## auxiliaries
  pi <- pi_fun(w_safe)
  r  <- (pi - par$i * d_safe) / par$sigma
  kap <- kappa_fun(r, par)
  gY  <- kap / par$sigma - par$delta
  g_n <- par$alpha + par$beta
  Phi <- Phi_fun(e_safe, d_safe, par)
  
  ## core dynamics
  de <- (gY - g_n) * e_safe
  dw <- (Phi - par$alpha) * w_safe
  
  ## debt dynamics (choose mode)
  mode <- if (!is.null(par$debt_mode)) par$debt_mode else "dynamic"  # "dynamic" or "frozen"
  if (mode == "frozen") {
    dd <- 0
  } else {
    tau_d <- if (!is.null(par$tau_d)) par$tau_d else 1
    dd_raw <- kap - pi + par$i * d_safe - d_safe * gY
    dd <- dd_raw / tau_d
    if (d <= d_floor && dd < 0) dd <- 0
  }
  
  ## per-capita output in levels (index)
  y_pc <- exp(par$alpha * t) * e_safe
  
  list(
    c(de, dw, dd),
    c(pi = pi, r = r, kappa = kap, gY = gY, Phi = Phi, y_pc = y_pc, g_n = g_n)
  )
}

simulate_model <- function(par,
                           state0 = NULL,
                           t_end = 800,
                           dt = 0.1,
                           perturb = c(e = 0.01, omega = -0.01, d = 0.10)) {
  
  if (is.null(state0)) {
    ss <- steady_state(par, strict = FALSE)
    if (!isTRUE(ss$ok)) stop("No admissible steady state: ", ss$reason)
    state0 <- c(
      e = ss$e + perturb[["e"]],
      omega = ss$omega + perturb[["omega"]],
      d = ss$d + perturb[["d"]]
    )
  }
  
  times <- seq(0, t_end, by = dt)
  
  sol <- ode(
    y = state0,
    times = times,
    func = rhs_keengm_slowfinance,
    parms = par,
    method = "lsoda"
  )
  
  as.data.frame(sol) |>
    as_tibble() |>
    mutate(
      y_pc_0 = first(y_pc),
      y_pc_norm = y_pc / y_pc_0
    )
}

## -----------------------------
## 4) Plots
## -----------------------------
plot_states <- function(df, tag, out_dir = "outputs") {
  p <- df |>
    select(time, e, omega, d) |>
    pivot_longer(-time, names_to = "var", values_to = "value") |>
    ggplot(aes(x = time, y = value)) +
    geom_line() +
    facet_wrap(~var, scales = "free_y", ncol = 1) +
    labs(title = paste0("States over time (", tag, ")"), x = "t", y = NULL)
  
  ggsave(file.path(out_dir, paste0("states_", tag, ".png")), p, width = 7, height = 7, dpi = 160)
  p
}

plot_finance <- function(df, tag, out_dir = "outputs") {
  p <- df |>
    select(time, r, pi, kappa, gY) |>
    pivot_longer(-time, names_to = "var", values_to = "value") |>
    ggplot(aes(x = time, y = value)) +
    geom_line() +
    facet_wrap(~var, scales = "free_y", ncol = 1) +
    labs(title = paste0("Finance block over time (", tag, ")"), x = "t", y = NULL)
  
  ggsave(file.path(out_dir, paste0("finance_", tag, ".png")), p, width = 7, height = 7, dpi = 160)
  p
}

plot_y_pc <- function(df, tag, out_dir = "outputs") {
  p <- ggplot(df, aes(x = time, y = y_pc_norm)) +
    geom_line() +
    labs(title = paste0("Per-capita output (levels, normalized) (", tag, ")"),
         x = "t", y = "y_pc (index, t=0=1)")
  
  ggsave(file.path(out_dir, paste0("y_pc_", tag, ".png")), p, width = 7, height = 4, dpi = 160)
  p
}

plot_3d_traj <- function(df, tag, out_dir = "outputs", thin = 5) {
  df2 <- df |> slice(seq(1, n(), by = thin))
  png(file.path(out_dir, paste0("traj3d_", tag, ".png")), width = 900, height = 700)
  scatterplot3d(df2$e, df2$omega, df2$d,
                type = "l",
                xlab = "e", ylab = "omega", zlab = "d",
                main = paste0("3D trajectory (e, omega, d) (", tag, ")"))
  dev.off()
  invisible(TRUE)
}

run_plots <- function(par, tag, t_end = 800, dt = 0.1, state0 = NULL, out_dir = "outputs") {
  dir_create(out_dir)
  df <- simulate_model(par, state0 = state0, t_end = t_end, dt = dt)
  plot_states(df, tag, out_dir)
  plot_finance(df, tag, out_dir)
  plot_y_pc(df, tag, out_dir)
  plot_3d_traj(df, tag, out_dir)
  invisible(df)
}

## -----------------------------
## 5) Example runs: “Goodwin fast + Finance slow”
## -----------------------------
par <- list(
  alpha = 0.02, beta = 0.01,
  delta = 0.02, sigma = 3,
  i = 0.04,
  kappa_min = 0.00, kappa_max = 0.30,
  lambda = 20, r0 = 0.04,
  phi0 = -0.06, phi1 = 0.10,
  eta = 0.02, dbar = 1.0, psi = 45,   # psi > 36 = stable side locally
  d_floor = 1e-9,
  tau_d = 30,                          # <<< slow finance knob
  debt_mode = "dynamic"                # "dynamic" or "frozen"
)

## Local diagnostics
out <- run_report(par)
cat(sprintf("\nFinite-diff dH/dpsi ≈ %.6g\n", dH_dmu_fd(par, "psi", h = 1e-2)))

## (A) Slow-finance dynamics: distributive cycles + long debt wave
df_slow <- run_plots(par, tag = "slowFinance_tau30_psi45", t_end = 1200, dt = 0.2)

## (B) Sanity check: freeze debt -> pure Goodwin-ish core (given this closure)
par$debt_mode <- "frozen"
df_frozen <- run_plots(par, tag = "debtFrozen_GoodwinCore", t_end = 400, dt = 0.05)

## (C) Fast finance comparison (optional): can get messy fast
par$debt_mode <- "dynamic"
par$tau_d <- 1
par$psi <- 10
df_fast <- run_plots(par, tag = "fastFinance_tau1_psi10", t_end = 300, dt = 0.05)
