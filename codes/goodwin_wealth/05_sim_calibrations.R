# ============================================================
# simulate_9_calibrations.R
# 9 simulations around the Hopf boundary with sigma fixed.
# States: (e, omega, d). Finance auxiliary: f_bar(r) (Option 1 reduction)
#
# Exports to:
# outputs/wealth_goodwin/grid_search/best_results/sims/
#   - scenario_summary.csv
#   - trajectories/scenario_XX.csv
#   - plots/time_*.png, phase_*.png
# ============================================================

suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

# ----------------------------
# Paths
# ----------------------------
base_dir <- "outputs/wealth_goodwin/grid_search/best_results"
sim_dir  <- file.path(base_dir, "sims")
traj_dir <- file.path(sim_dir, "trajectories")
plot_dir <- file.path(sim_dir, "plots")
dir.create(sim_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(traj_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helper utilities
# ----------------------------
clamp01 <- function(x) pmin(1, pmax(0, x))
rf_tilde_from_rf <- function(rF, psi) -clamp01(psi) * rF

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
# Base parameters (match your calibration code)
# ----------------------------
par_base <- list(
  sigma = 2.0,
  delta = 0.05,
  alpha = 0.02,
  beta  = 0.01,
  i     = 0.03,
  
  phi0 = 0.00,
  phi1 = 0.60,  # will be computed endogenously
  phi2 = 0.70,
  phi3 = 5.0,
  phi4 = 1.0,
  
  kappa_min = 0.01,
  kappa_max = 0.35,
  kappa0    = 0.05,
  kappa1    = 40.0,
  
  lam0 = 0.00,
  lam1 = 40.0,
  
  rF = 0.08,
  rF_tilde = 0.00,
  psi = 0.20
)

# ----------------------------
# Model primitives
# ----------------------------
kappa_fun <- function(r, p) logistic_bounded(r, p$kappa_min, p$kappa_max, p$kappa0, p$kappa1)
kappa_r   <- function(r, p) logistic_bounded_d(r, p$kappa_min, p$kappa_max, p$kappa0, p$kappa1)

g_fun <- function(r, p) kappa_fun(r, p) / p$sigma - p$delta

lambda_fun <- function(r, p) {
  rF_eff <- p$rF + p$rF_tilde
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
# Backbone steady state (pinned by sigma, gn, delta via kappa*)
# ----------------------------
steady_state_backbone <- function(p) {
  gn <- p$alpha + p$beta
  if (!is.finite(gn) || gn <= 0) return(list(ok = FALSE, reason = "gn_nonpositive"))
  
  kappa_star <- p$sigma * (gn + p$delta)
  if (!(kappa_star > p$kappa_min && kappa_star < p$kappa_max)) {
    return(list(ok = FALSE, reason = "kappa_star_out_of_bounds"))
  }
  
  r_star <- logistic_bounded_inv(kappa_star, p$kappa_min, p$kappa_max, p$kappa0, p$kappa1)
  d_star <- (kappa_star - p$sigma * r_star) / gn
  omega_star <- 1 - p$i * d_star - p$sigma * r_star
  
  ok <- all(is.finite(c(r_star, d_star, omega_star))) && d_star >= 0 && omega_star > 0 && omega_star < 1
  list(
    ok = ok,
    reason = if (ok) "OK" else "backbone_invalid",
    gn = gn,
    kappa_star = kappa_star,
    r_star = r_star,
    d_star = d_star,
    omega_star = omega_star
  )
}

# ----------------------------
# Reduction: f_bar(r) (Option 1 logic)
# f_bar = iotaF/gn, iotaF = lam * sF, sF = r/(1-lam)
# ----------------------------
f_bar <- function(r, p, gn) {
  lam <- lambda_fun(r, p)
  if (!is.finite(lam) || lam >= 0.999) return(NA_real_)
  sF <- r / (1 - lam)
  iotaF <- lam * sF
  iotaF / gn
}

# ----------------------------
# Full SS + phi1 endogenous to hit e_target
# ----------------------------
steady_state_full <- function(p, bb, phi2, e_target) {
  gn <- bb$gn
  r_star <- bb$r_star
  d_star <- bb$d_star
  omega_star <- bb$omega_star
  
  # finance at r_star
  f_star <- f_bar(r_star, p, gn)
  if (!is.finite(f_star) || f_star < 0) return(list(ok = FALSE, reason = "f_star_invalid"))
  
  p$phi2 <- phi2
  Z_star <- Z_fun(d_star, f_star, p)
  
  numer <- (p$alpha - p$phi0 + phi2 * Z_star)
  if (!is.finite(numer) || numer <= 0) return(list(ok = FALSE, reason = "phi1_numer_nonpositive"))
  phi1 <- numer / e_target
  e_star <- numer / phi1
  
  ok <- is.finite(phi1) && is.finite(e_star) && phi1 > 0 && e_star > 0
  list(
    ok = ok,
    reason = if (ok) "OK" else "ss_invalid",
    r_star = r_star,
    d_star = d_star,
    omega_star = omega_star,
    f_star = f_star,
    Z_star = Z_star,
    phi1 = phi1,
    e_star = e_star
  )
}

# ----------------------------
# Jacobian at SS (for local classification)
# ----------------------------
jacobian_at_ss <- function(ss, bb, p, phi2) {
  gn <- bb$gn
  r  <- ss$r_star
  d  <- ss$d_star
  e  <- ss$e_star
  om <- ss$omega_star
  phi1 <- ss$phi1
  
  r_om <- -1 / p$sigma
  r_d  <- -p$i / p$sigma
  
  kap_rv <- kappa_r(r, p)
  g_r   <- kap_rv / p$sigma
  
  # f_r under reduction f_bar(r)
  lam <- lambda_fun(r, p)
  lam_rv <- lambda_r(r, p)
  if (lam >= 0.999) lam <- 0.999
  # f = [lam*r/(1-lam)]/gn
  # derivative wrt r:
  # f_r = (1/gn)* d/dr [lam*r/(1-lam)]
  # Let q = lam/(1-lam). Then f = (r*q)/gn.
  # q_r = lam_r/(1-lam)^2
  q <- lam / (1 - lam)
  q_r <- lam_rv / ((1 - lam)^2)
  f_r <- (q + r * q_r) / gn
  
  f <- ss$f_star
  Zd <- Z_d_fun(d, f, p)
  Zf <- Z_f_fun(d, f, p)
  
  Z_om <- Zf * f_r * r_om
  Z_dv <- Zd + Zf * f_r * r_d
  
  # \dot e = (g(r) - gn) e
  J11 <- 0
  J12 <- e * g_r * r_om
  J13 <- e * g_r * r_d
  
  # \dot omega = omega*(phi0 + phi1 e - alpha - phi2 Z)
  J21 <- om * phi1
  J22 <- -om * phi2 * Z_om
  J23 <- -om * phi2 * Z_dv
  
  # \dot d = kappa(r) - (1-omega) + i d - d*g(r)
  # partials:
  # d/domega: +1 + kappa_r*r_om - d*g_r*r_om
  # d/dd: i - gn? (no: g depends on r(d), and r_d enters) and -g - d*g_r*r_d
  kap_r_term_om <- kap_rv * r_om
  kap_r_term_d  <- kap_rv * r_d
  
  # g(r) level
  g <- g_fun(r, p)
  
  J31 <- 0
  J32 <- 1 + kap_r_term_om - d * g_r * r_om
  J33 <- kap_r_term_d + p$i - g - d * g_r * r_d
  
  matrix(c(J11,J12,J13,
           J21,J22,J23,
           J31,J32,J33), nrow = 3, byrow = TRUE)
}

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
  tibble(a1=a1, a2=a2, a3=a3, H=H, tr=tr, det=detJ)
}

# ----------------------------
# ODE system (3D with f_bar(r) as auxiliary)
# ----------------------------
ode_3d <- function(t, state, pars) {
  e <- state[["e"]]
  om <- state[["omega"]]
  d <- state[["d"]]
  
  p <- pars$p
  gn <- pars$gn
  
  # profit rate proxy
  r <- (1 - om - p$i * d) / p$sigma
  
  kap <- kappa_fun(r, p)
  g <- kap / p$sigma - p$delta
  
  f <- f_bar(r, p, gn)
  if (!is.finite(f)) f <- 1e6  # force blow-up if invalid
  
  Z <- Z_fun(d, f, p)
  
  de <- (g - gn) * e
  domega <- om * (p$phi0 + p$phi1 * e - p$alpha - p$phi2 * Z)
  dd <- kap - (1 - om) + p$i * d - d * g
  
  list(c(de, domega, dd),
       c(r=r, kappa=kap, g=g, f=f, Z=Z))
}

# ============================================================
# 9 calibrations (the “cross” around Hopf)
# Choose psi narratively; use rF_eff as control knob
# ============================================================
psi_choice <- 0.20

rF_eff_vals <- c(0.038, 0.040, 0.043)
phi2_vals   <- c(0.60, 0.70, 0.80)

scenarios <- expand_grid(rF_eff = rF_eff_vals, phi2 = phi2_vals) %>%
  arrange(phi2, rF_eff) %>%
  mutate(
    scenario_id = sprintf("S%02d", row_number()),
    psi = psi_choice,
    rF = rF_eff / (1 - psi),
    rF_tilde = rf_tilde_from_rf(rF, psi)
  )

write_csv(scenarios, file.path(sim_dir, "scenarios.csv"))

# ============================================================
# Simulation settings
# ============================================================
times <- seq(0, 400, by = 0.1)  # adjust horizon if you want longer
pert  <- 0.01                  # 1% perturbation around SS

# ============================================================
# Run all scenarios
# ============================================================
summary_rows <- list()

for (ii in seq_len(nrow(scenarios))) {
  
  sc <- scenarios[ii, ]
  
  # build parameter list for scenario
  p <- par_base
  p$psi <- sc$psi
  p$rF  <- sc$rF
  p$rF_tilde <- sc$rF_tilde
  p$phi2 <- sc$phi2
  
  # backbone
  bb <- steady_state_backbone(p)
  if (!isTRUE(bb$ok)) {
    summary_rows[[ii]] <- tibble(
      scenario_id=sc$scenario_id, ok=FALSE, fail="backbone",
      rF_eff=sc$rF_eff, psi=sc$psi, rF=sc$rF, phi2=sc$phi2
    )
    next
  }
  
  # SS + phi1 to hit target e*
  ss <- steady_state_full(p, bb, phi2=sc$phi2, e_target=0.94)
  if (!isTRUE(ss$ok)) {
    summary_rows[[ii]] <- tibble(
      scenario_id=sc$scenario_id, ok=FALSE, fail=ss$reason,
      rF_eff=sc$rF_eff, psi=sc$psi, rF=sc$rF, phi2=sc$phi2
    )
    next
  }
  p$phi1 <- ss$phi1  # IMPORTANT: use endogenous phi1
  
  # local stability (Jacobian)
  J <- jacobian_at_ss(ss, bb, p, phi2=sc$phi2)
  ev <- eigen(J, only.values = TRUE)$values
  maxRe <- max(Re(ev))
  maxIm <- max(abs(Im(ev)))
  rh <- rh_coeffs(J)
  local_regime <- ifelse(maxRe < 0, "stable", "unstable")
  local_osc <- ifelse(maxIm > 1e-6, TRUE, FALSE)
  
  # initial condition near SS
  y0 <- c(
    e     = ss$e_star * (1 + pert),
    omega = ss$omega_star * (1 - pert),
    d     = ss$d_star * (1 + pert)
  )
  
  out <- ode(
    y = y0, times = times, func = ode_3d,
    parms = list(p=p, gn=bb$gn),
    method = "lsoda"
  )
  
  df <- as.data.frame(out) %>%
    as_tibble() %>%
    mutate(
      scenario_id = sc$scenario_id,
      rF_eff = sc$rF_eff, psi=sc$psi, rF=sc$rF,
      phi2 = sc$phi2, phi1 = p$phi1,
      e_star = ss$e_star, omega_star = ss$omega_star, d_star = ss$d_star
    )
  
  # quick boundedness flags
  df <- df %>%
    mutate(
      bad = (e <= 0) | (omega <= 0) | (omega >= 1) | (d < 0) | is.na(r) | is.na(f) | (abs(f) > 1e5)
    )
  
  diverged <- any(df$bad)
  
  # tail diagnostics (last 30%)
  n <- nrow(df)
  i0 <- max(1, floor(0.70 * n))
  tail <- df[i0:n, ]
  
  # amplitude & period proxy
  amp_e  <- max(tail$e, na.rm=TRUE) - min(tail$e, na.rm=TRUE)
  amp_om <- max(tail$omega, na.rm=TRUE) - min(tail$omega, na.rm=TRUE)
  
  # crude cycle detector: std in tail
  sd_e  <- sd(tail$e, na.rm=TRUE)
  sd_om <- sd(tail$omega, na.rm=TRUE)
  
  # export trajectory
  traj_file <- file.path(traj_dir, paste0("scenario_", sc$scenario_id, ".csv"))
  write_csv(df, traj_file)
  
  # plots: time series
  p_time <- df %>%
    select(time, e, omega, d) %>%
    pivot_longer(-time) %>%
    ggplot(aes(x=time, y=value)) +
    geom_line() +
    facet_wrap(~name, scales="free_y", ncol=1) +
    labs(
      title = paste0(sc$scenario_id,
                     " | rF_eff=", sprintf("%.3f", sc$rF_eff),
                     " psi=", sprintf("%.2f", sc$psi),
                     " phi2=", sprintf("%.2f", sc$phi2),
                     " (local ", local_regime, ")"),
      x="time", y=""
    ) +
    theme_minimal(base_size = 12)
  
  ggsave(file.path(plot_dir, paste0("time_", sc$scenario_id, ".png")),
         p_time, width=7, height=7, dpi=220)
  
  # plots: phase (omega vs e)
  p_phase <- df %>%
    ggplot(aes(x=e, y=omega)) +
    geom_path(alpha=0.8) +
    geom_point(data = df %>% slice(1), aes(x=e, y=omega), size=2) +
    labs(
      title = paste0("Phase: omega vs e | ", sc$scenario_id,
                     " | rF_eff=", sprintf("%.3f", sc$rF_eff),
                     " phi2=", sprintf("%.2f", sc$phi2)),
      x="e", y="omega"
    ) +
    theme_minimal(base_size = 12)
  
  ggsave(file.path(plot_dir, paste0("phase_omega_e_", sc$scenario_id, ".png")),
         p_phase, width=7, height=5, dpi=220)
  
  # record summary
  summary_rows[[ii]] <- tibble(
    scenario_id = sc$scenario_id,
    ok = TRUE,
    diverged = diverged,
    rF_eff = sc$rF_eff, psi = sc$psi, rF = sc$rF,
    phi2 = sc$phi2, phi1 = p$phi1,
    e_star = ss$e_star, omega_star = ss$omega_star, d_star = ss$d_star,
    local_maxReEig = maxRe,
    local_maxImEig = maxIm,
    local_regime = paste0(local_regime, ifelse(local_osc, "_osc", "_real")),
    a1 = rh$a1, a2=rh$a2, a3=rh$a3, H=rh$H,
    tail_amp_e = amp_e,
    tail_amp_omega = amp_om,
    tail_sd_e = sd_e,
    tail_sd_omega = sd_om
  )
}

summary_df <- bind_rows(summary_rows)

write_csv(summary_df, file.path(sim_dir, "scenario_summary.csv"))

# A quick “shortlist” filter: not diverged, decent oscillation in tail
shortlist <- summary_df %>%
  filter(ok, !diverged) %>%
  arrange(desc(tail_sd_omega), desc(tail_sd_e), abs(local_maxReEig)) %>%
  mutate(candidate = (tail_sd_omega > 0.005 | tail_sd_e > 0.005))

write_csv(shortlist, file.path(sim_dir, "scenario_shortlist.csv"))

writeLines(capture.output(sessionInfo()), file.path(sim_dir, "sessionInfo.txt"))

message("Done. Outputs in: ", sim_dir)
message("Key files: scenario_summary.csv, scenario_shortlist.csv, trajectories/*.csv, plots/*.png")
