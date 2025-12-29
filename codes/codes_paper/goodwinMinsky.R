############################################################
## keen_goodwin_minsky_3d_full.R
## Keen-orthodox Goodwin–Minsky 3D core + Buffett discipline
##
## Includes:
##  - Closed-form interior steady state (admissibility gates)
##  - Analytic Jacobian at steady state
##  - Eigenvalues
##  - Routh–Hurwitz coefficients + Hopf functional H
##  - Hopf targeting for psi (when d* ~ dbar)
##  - Finite-diff dH/dmu with steady-state re-solve
##  - Simulation + plots (states, finance, y_pc, 3D traj)
############################################################

## -----------------------------
## 0) Packages (used only for simulation/plots)
## -----------------------------
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
## 1) Helpers
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

## Profit share and net profit rate
pi_fun <- function(omega) 1 - omega

r_fun <- function(omega, d, par) {
  (pi_fun(omega) - par$i * d) / par$sigma
}

## Investment (logistic Keen)
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
## 2) Interior steady state (closed form)
## -----------------------------
steady_state <- function(par,
                         strict = FALSE,
                         e_band = c(0.75, 0.95),
                         eps = 1e-12) {
  
  require_params(par, c(
    "alpha","beta","delta","sigma","i",
    "kappa_min","kappa_max","lambda","r0",
    "phi0","phi1",
    "eta","dbar","psi"
  ))
  
  g_n <- par$alpha + par$beta
  if (g_n <= 0) return(list(ok = FALSE, reason = "g_n <= 0 (alpha+beta must be >0)"))
  
  k_target <- par$sigma * (g_n + par$delta)
  
  ## Gate 1: k_target within logistic bounds
  if (k_target <= par$kappa_min + eps || k_target >= par$kappa_max - eps) {
    return(list(ok = FALSE, reason = "kappa_target outside (kappa_min, kappa_max)"))
  }
  
  s <- (k_target - par$kappa_min) / (par$kappa_max - par$kappa_min)
  if (s <= eps || s >= 1 - eps) {
    return(list(ok = FALSE, reason = "s not in (0,1)"))
  }
  
  ## r* from logistic inversion
  r_star <- par$r0 + (1 / par$lambda) * log(s / (1 - s))
  
  ## d* from identity: r* = (k_target - d*g_n)/sigma
  d_star <- (k_target - par$sigma * r_star) / g_n
  
  ## Gate 2: d*>0
  if (d_star <= eps) {
    return(list(ok = FALSE, reason = "d* <= 0 (no interior debt steady state)"))
  }
  
  ## profit share and wage share
  pi_star <- k_target + d_star * (par$i - g_n)
  omega_star <- 1 - pi_star
  
  ## Gate 3: omega in (0,1)
  if (omega_star <= eps || omega_star >= 1 - eps) {
    return(list(ok = FALSE, reason = "omega* not in (0,1)"))
  }
  
  ## e* from Phillips block
  if (par$phi1 <= 0) return(list(ok = FALSE, reason = "phi1 must be > 0"))
  e_star <- (par$alpha - par$phi0 + par$eta * S_fun(d_star, par)) / par$phi1
  
  ## Gate 4: e in (0,1)
  if (e_star <= eps || e_star >= 1 - eps) {
    return(list(ok = FALSE, reason = "e* not in (0,1)"))
  }
  
  ## Optional plausibility band
  if (strict && (e_star < e_band[1] || e_star > e_band[2])) {
    return(list(ok = FALSE, reason = sprintf(
      "e* outside plausibility band [%.2f, %.2f]",
      e_band[1], e_band[2]
    )))
  }
  
  list(
    ok = TRUE,
    e = e_star,
    omega = omega_star,
    d = d_star,
    r = r_star,
    kappa = k_target,
    g_n = g_n
  )
}

## -----------------------------
## 3) Analytic Jacobian at steady state
## -----------------------------
jacobian_at_ss <- function(ss, par) {
  e <- ss$e; omega <- ss$omega; d <- ss$d
  g_n <- ss$g_n
  r <- ss$r
  
  kap_r <- kappa_r_fun(r, par)
  Sd <- S_d_fun(d, par)
  
  ## Sparse form J = [0 A B; C 0 D; 0 E F]
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
  
  list(
    J = J,
    A = A, B = B, C = C, D = D, E = E, F = F,
    kappa_r = kap_r, S_d = Sd
  )
}

## -----------------------------
## 4) RH coefficients + Hopf functional
## -----------------------------
rh_hopf <- function(Jparts) {
  A <- Jparts$A; B <- Jparts$B; C <- Jparts$C
  D <- Jparts$D; E <- Jparts$E; F <- Jparts$F
  
  ## det(sI - J) = s^3 + a1 s^2 + a2 s + a3
  a1 <- -F
  a2 <- -A * C - E * D
  a3 <- C * (A * F - B * E)
  
  H <- a1 * a2 - a3
  stable_RH <- (a1 > 0) && (a2 > 0) && (a3 > 0) && (a1 * a2 > a3)
  
  ## Sparse identity check: H = E*(F*D + C*B)
  H_sparse <- E * (F * D + C * B)
  
  list(a1 = a1, a2 = a2, a3 = a3, H = H, H_sparse = H_sparse, stable_RH = stable_RH)
}

## -----------------------------
## 5) One-shot analysis
## -----------------------------
analyze_model <- function(par, strict = FALSE, e_band = c(0.75, 0.95)) {
  ss <- steady_state(par, strict = strict, e_band = e_band)
  if (!isTRUE(ss$ok)) return(ss)
  
  Jp <- jacobian_at_ss(ss, par)
  rh <- rh_hopf(Jp)
  eig <- eigen(Jp$J)$values
  
  list(
    ok = TRUE,
    steady_state = ss,
    jacobian = Jp$J,
    jac_parts = Jp[c("A","B","C","D","E","F","kappa_r","S_d")],
    eigenvalues = eig,
    rh = rh
  )
}

## -----------------------------
## 6) Hopf gap and dH/dmu (finite diff, re-solving SS)
## -----------------------------
hopf_gap <- function(par, strict = FALSE, e_band = c(0.75, 0.95)) {
  out <- analyze_model(par, strict = strict, e_band = e_band)
  if (!isTRUE(out$ok)) return(NA_real_)
  out$rh$H
}

dH_dmu_fd <- function(par, mu, h = 1e-3, strict = FALSE, e_band = c(0.75, 0.95)) {
  if (!mu %in% names(par)) stop("Parameter '", mu, "' not in par list.")
  if (!is.finite(par[[mu]])) stop("par[['", mu, "']] must be finite.")
  
  par_p <- par; par_m <- par
  par_p[[mu]] <- par[[mu]] + h
  par_m[[mu]] <- par[[mu]] - h
  
  Hp <- hopf_gap(par_p, strict = strict, e_band = e_band)
  Hm <- hopf_gap(par_m, strict = strict, e_band = e_band)
  
  if (is.na(Hp) || is.na(Hm)) return(NA_real_)
  (Hp - Hm) / (2*h)
}

## -----------------------------
## 7) Hopf-targeting for psi (best when d* ~ dbar)
## -----------------------------
psi_hopf_guess <- function(out, par, tol = 1e-6) {
  if (!isTRUE(out$ok)) stop("Provide a valid analyze_model() output.")
  
  ss <- out$steady_state
  Jp <- out$jac_parts
  
  omega <- ss$omega
  d <- ss$d
  
  C <- Jp$C; B <- Jp$B; F <- Jp$F
  if (abs(F) < 1e-12) return(NA_real_)
  
  ## Target D on Hopf boundary: FD + CB = 0 => D_H = -(C*B)/F
  D_H <- - (C * B) / F
  
  S <- S_fun(d, par)
  Sfac <- S * (1 - S)
  if (Sfac < tol) return(NA_real_)
  
  ## D = -omega*eta*(psi/dbar)*S(1-S) => psi = (-D_H)*dbar/(omega*eta*S(1-S))
  psi_H <- (-D_H) * par$dbar / (omega * par$eta * Sfac)
  psi_H
}

## -----------------------------
## 8) Quick report printer
## -----------------------------
run_report <- function(par, strict = FALSE, e_band = c(0.75, 0.95)) {
  out <- analyze_model(par, strict = strict, e_band = e_band)
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
  cat(sprintf("a1=%.6g a2=%.6g a3=%.6g\n", rh$a1, rh$a2, rh$a3))
  cat(sprintf("H=%.6g (sparse check %.6g) | stable_RH=%s\n",
              rh$H, rh$H_sparse, rh$stable_RH))
  
  psi_guess <- psi_hopf_guess(out, par)
  if (is.finite(psi_guess)) {
    cat(sprintf("\nHopf psi-guess (from FD+CB=0): psi_H ≈ %.4f\n", psi_guess))
  } else {
    cat("\nHopf psi-guess: NA (switch saturated or F ~ 0)\n")
  }
  
  invisible(out)
}

############################################################
## 9) Simulation + plots
############################################################

## ODE RHS (Keen orthodox, no resets)
rhs_keengm <- function(t, state, par) {
  e <- state[["e"]]
  w <- state[["omega"]]
  d <- state[["d"]]
  
  ## guardrails (numerical only)
  d_floor <- if (!is.null(par$d_floor)) par$d_floor else 1e-9
  d_safe  <- max(d, d_floor)
  e_safe  <- min(max(e, 1e-9), 1 - 1e-9)
  w_safe  <- min(max(w, 1e-9), 1 - 1e-9)
  
  ## auxiliaries
  pi <- pi_fun(w_safe)
  r  <- (pi - par$i * d_safe) / par$sigma
  kap <- kappa_fun(r, par)
  gY  <- kap / par$sigma - par$delta
  g_n <- par$alpha + par$beta
  Phi <- Phi_fun(e_safe, d_safe, par)
  
  ## dynamics
  de <- (gY - g_n) * e_safe
  dw <- (Phi - par$alpha) * w_safe
  dd <- kap - pi + par$i * d_safe - d_safe * gY
  
  ## prevent negative drift at the floor
  if (d <= d_floor && dd < 0) dd <- 0
  
  ## per-capita output level (index)
  y_pc <- exp(par$alpha * t) * e_safe
  
  list(
    c(de, dw, dd),
    c(pi = pi, r = r, kappa = kap, gY = gY, Phi = Phi, y_pc = y_pc, g_n = g_n)
  )
}

simulate_model <- function(par, state0 = NULL,
                           t_end = 250, dt = 0.05,
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
    func = rhs_keengm,
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

plot_states <- function(df, tag, out_dir = "outputs/minsky/") {
  p <- df |>
    select(time, e, omega, d) |>
    pivot_longer(-time, names_to = "var", values_to = "value") |>
    ggplot(aes(x = time, y = value)) +
    geom_line() +
    facet_wrap(~var, scales = "free_y", ncol = 1) +
    labs(title = paste0("States over time (", tag, ")"), x = "t", y = NULL)
  
  ggsave(file.path(out_dir, paste0("states_", tag, ".png")),
         p, width = 7, height = 7, dpi = 160)
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
  
  ggsave(file.path(out_dir, paste0("finance_", tag, ".png")),
         p, width = 7, height = 7, dpi = 160)
  p
}

plot_y_pc <- function(df, tag, out_dir = "outputs") {
  p <- ggplot(df, aes(x = time, y = y_pc_norm)) +
    geom_line() +
    labs(title = paste0("Per-capita output level (normalized) (", tag, ")"),
         x = "t", y = "y_pc (index, t=0 = 1)")
  
  ggsave(file.path(out_dir, paste0("y_pc_", tag, ".png")),
         p, width = 7, height = 4, dpi = 160)
  p
}

plot_3d_traj <- function(df, tag, out_dir = "outputs", thin = 5) {
  df2 <- df |> slice(seq(1, n(), by = thin))
  
  png(file.path(out_dir, paste0("traj3d_", tag, ".png")), width = 900, height = 700)
  scatterplot3d(
    df2$e, df2$omega, df2$d,
    type = "l",
    xlab = "e", ylab = "omega", zlab = "d",
    main = paste0("3D trajectory (e, omega, d) (", tag, ")")
  )
  dev.off()
  
  invisible(TRUE)
}

run_plots <- function(par, tag, t_end = 250, dt = 0.05,
                      state0 = NULL, out_dir = "outputs") {
  
  dir_create(out_dir)
  
  df <- simulate_model(par, state0 = state0, t_end = t_end, dt = dt)
  
  plot_states(df, tag, out_dir)
  plot_finance(df, tag, out_dir)
  plot_y_pc(df, tag, out_dir)
  plot_3d_traj(df, tag, out_dir)
  
  invisible(df)
}

############################################################
## 10) Example run (your working set)
############################################################
par <- list(
  alpha = 0.02, beta = 0.01,
  delta = 0.02, sigma = 3,
  i = 0.04,
  kappa_min = 0.00, kappa_max = 0.30,
  lambda = 20, r0 = 0.04,
  phi0 = -0.06, phi1 = 0.10,
  eta = 0.02, dbar = 1.0, psi = 10,
  d_floor = 1e-9
)

## Diagnostics
out <- run_report(par, strict = FALSE)
dHpsi <- dH_dmu_fd(par, "psi", h = 1e-2)
cat(sprintf("\nFinite-diff dH/dpsi ≈ %.6g\n", dHpsi))

## Plots: unstable (psi below Hopf)
par$psi <- 10
df_B <- run_plots(par, tag = "psi10_unstable", t_end = 250, dt = 0.05)

## Plots: stable (psi above Hopf)
par$psi <- 45
df_A <- run_plots(par, tag = "psi45_stable", t_end = 250, dt = 0.05)






