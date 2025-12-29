############################################################
## 02_goodwinUnionAffiliation.R
## Minimal Goodwin (2D) + Union Affiliation as bifurcation:
##   mu(s) = zeta * (s - s0)
##
## System (in original Goodwin variables, deviations auxiliary):
##   de/dt     = -a * (omega - omega*)
##   domega/dt =  b * (e - e*) + mu(s)*(omega - omega*) - c*(omega - omega*)^3
##
## Features:
## - Interior fixed point is (e*, omega*) for any fixed s
## - Scan in s for Hopf boundary (trace=mu crosses 0, det=ab>0)
## - Simulate for {below, near, above} s0
## - Uses logit states to keep e, omega in (0,1)
## - Saves figures to ./outputs
############################################################

suppressPackageStartupMessages({
  library(deSolve)
  library(numDeriv)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(fs)
})

## ---------------------------------------------------------
## 0) Utility
## ---------------------------------------------------------
inv_logit <- function(x) 1 / (1 + exp(-x))
logit     <- function(p) log(p / (1 - p))
clamp01   <- function(u, eps = 1e-10) pmin(pmax(u, eps), 1 - eps)

dir_create("outputs")

## ---------------------------------------------------------
## 1) Model blocks (Goodwin + institutional flip)
## ---------------------------------------------------------
mu_fun <- function(s, par) par$zeta * (s - par$s0)

F_vec <- function(x, par) {
  e     <- x[1]
  omega <- x[2]
  s     <- par$s_fixed
  
  mu <- mu_fun(s, par)
  
  ## deviations (auxiliary)
  u  <- e     - par$e_star
  v  <- omega - par$omega_star
  r2 <- u^2 + v^2
  
  ## Normal-form style radial saturation:
  ## - de keeps Goodwin core (-a*v) but adds radial damping on u
  ## - domega keeps Phillips core (b*u) + institutional flip (mu*v) + radial damping on v
  de     <- -par$a * v - par$c * u * r2
  domega <-  par$b * u + mu * v - par$c * v * r2
  
  c(de, domega)
}

## ---------------------------------------------------------
## 2) Fixed point (trivial here)
## ---------------------------------------------------------
fixed_point_solver <- function(par) {
  list(
    e_star     = par$e_star,
    omega_star = par$omega_star,
    s_fixed    = par$s_fixed,
    mu         = mu_fun(par$s_fixed, par)
  )
}

## ---------------------------------------------------------
## 3) Local diagnostics (Jacobian + eigenvalues + Hopf flags)
## ---------------------------------------------------------
local_diagnostics <- function(par) {
  fp <- fixed_point_solver(par)
  x_star <- c(fp$e_star, fp$omega_star)
  
  J <- jacobian(func = function(x) F_vec(x, par), x = x_star)
  ev <- eigen(J, only.values = TRUE)$values
  
  trJ  <- sum(diag(J))
  detJ <- det(J)
  
  ## Hopf boundary in 2D: tr=0 and det>0 (transversality holds here)
  list(
    fp=fp, J=J, eigenvalues=ev,
    tr=trJ, det=detJ,
    hopf_boundary = (abs(trJ) < 1e-8 && detJ > 0)
  )
}

## ---------------------------------------------------------
## 4) ODE in transformed states (logit(e), logit(omega))
## ---------------------------------------------------------
goodwin_union_logstates <- function(t, state, par) {
  x <- state["x"]
  y <- state["y"]
  
  e     <- inv_logit(x)
  omega <- inv_logit(y)
  
  vec <- F_vec(c(e, omega), par)
  de <- vec[1]; domega <- vec[2]
  
  ## logit derivatives
  dx <- de     / (e * (1 - e))
  dy <- domega / (omega * (1 - omega))
  
  list(c(dx, dy))
}

## ---------------------------------------------------------
## 5) Calibration
## ---------------------------------------------------------
## Choose the interior reference point (Goodwin "center")
par0 <- list(
  ## Fixed point location (pick targets)
  e_star     = 0.90,
  omega_star = 0.65,
  
  ## Goodwin feedback strengths
  a = 0.20,   # omega -> e
  b = 0.25,   # e -> omega
  
  ## Saturation (bounds the cycle when mu>0)
  c = 5.0,
  
  ## Institutional flip
  s0   = 0.20,  # threshold unionization (Hopf boundary at s = s0)
  zeta = 1.50,  # flip strength
  
  ## s is treated as a fixed parameter for each run (set later)
  s_fixed = NA_real_
)

cat("\n--- Baseline diagnostics at s = s0 (should be Hopf boundary) ---\n")
par_tmp <- par0; par_tmp$s_fixed <- par0$s0
base <- local_diagnostics(par_tmp)
print(base$fp)
cat("Trace:", base$tr, "Det:", base$det, "\n")
cat("Eigenvalues:\n"); print(base$eigenvalues)

## ---------------------------------------------------------
## 6) Scan in s for stability + Hopf boundary
## ---------------------------------------------------------
s_grid <- seq(0.01, 0.99, by = 0.005)

scan <- lapply(s_grid, function(sval) {
  par <- par0
  par$s_fixed <- sval
  
  out <- local_diagnostics(par)
  tibble(
    s = sval,
    mu = out$fp$mu,
    tr = out$tr,
    det = out$det,
    maxRe = max(Re(out$eigenvalues)),
    complex = any(abs(Im(out$eigenvalues)) > 1e-10)
  )
}) %>% bind_rows()

p_mu <- ggplot(scan, aes(s, mu)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.4) +
  geom_line(linewidth=0.7) +
  labs(title="Institutional flip term μ(s) = ζ(s - s0)",
       subtitle="Hopf boundary at μ=0 ⇔ s=s0",
       x="Unionization s", y="μ(s)") +
  theme_minimal()

ggsave("outputs/union/fig_mu_scan_s.png", p_mu, width=9, height=4, dpi=300)

p_stab <- ggplot(scan, aes(s, maxRe)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.4) +
  geom_line(linewidth=0.7) +
  labs(title="Local stability proxy: max Re(eigenvalue) at (e*, ω*)",
       subtitle="maxRe<0 stable, maxRe>0 unstable; crossing at s=s0",
       x="Unionization s", y="max Re(λ)") +
  theme_minimal()

ggsave("outputs/union/fig_stability_scan_s.png", p_stab, width=9, height=4, dpi=300)

## ---------------------------------------------------------
## 7) Choose {below, near, above} around s0 and simulate
## ---------------------------------------------------------
s_set <- c(
  below = max(0.01, par0$s0 - 0.05),
  near  = par0$s0,
  above = min(0.99, par0$s0 + 0.05)
)

cat("\n--- s values to simulate ---\n")
print(s_set)

## Initial conditions: small perturbation near the fixed point
e0     <- clamp01(par0$e_star + 0.01)
omega0 <- clamp01(par0$omega_star - 0.01)
state0 <- c(x = logit(e0), y = logit(omega0))

times <- seq(0, 800, by = 0.05)

run_one <- function(sval, tag) {
  par <- par0
  par$s_fixed <- sval
  
  sol <- ode(
    y = state0,
    times = times,
    func = goodwin_union_logstates,
    parms = par,
    method = "lsoda",
    rtol = 1e-9,
    atol = 1e-11
  )
  
  df <- as.data.frame(sol) %>%
    mutate(
      e     = inv_logit(x),
      omega = inv_logit(y),
      s     = sval,
      mu    = mu_fun(sval, par),
      tag   = tag
    ) %>%
    select(time, tag, s, mu, e, omega)
  
  df
}

df_all <- bind_rows(
  run_one(s_set["below"], "below"),
  run_one(s_set["near"],  "near"),
  run_one(s_set["above"], "above")
)

## ---------------------------------------------------------
## 8) Plots
## ---------------------------------------------------------
df_long <- df_all %>%
  pivot_longer(cols=c(e, omega), names_to="var", values_to="value") %>%
  mutate(var = factor(var, levels=c("e","omega"),
                      labels=c("Employment rate (e)", "Wage share (ω)")))

p_time <- ggplot(df_long, aes(time, value)) +
  geom_line(linewidth=0.6) +
  facet_grid(var ~ tag, scales="free_y") +
  labs(title="Goodwin with Union Affiliation: time paths around s0",
       subtitle="below: μ<0 damped; near: μ≈0 marginal; above: μ>0 bounded cycle (via cubic saturation)",
       x="time", y="") +
  theme_minimal() +
  theme(legend.position="none")

ggsave("outputs/union/fig_timepaths_e_omega.png", p_time, width=12, height=6, dpi=300)

## Phase portrait
p_phase <- ggplot(df_all, aes(e, omega)) +
  geom_path(linewidth=0.5) +
  facet_wrap(~tag, scales="free") +
  labs(title="Phase portrait: (e, ω) trajectories by institutional regime",
       x="e", y="ω") +
  theme_minimal()

ggsave("outputs/union/fig_phase_e_omega.png", p_phase, width=10, height=4, dpi=300)

## Quick ranges
summ <- df_all %>%
  group_by(tag) %>%
  summarise(
    s = first(s),
    mu = first(mu),
    e_min = min(e), e_max = max(e),
    omega_min = min(omega), omega_max = max(omega),
    .groups="drop"
  )

cat("\n--- Range summary by regime tag ---\n")
print(summ)

cat("\nSaved outputs in ./outputs:\n",
    "- fig_mu_scan_s.png\n",
    "- fig_stability_scan_s.png\n",
    "- fig_timepaths_e_omega.png\n",
    "- fig_phase_e_omega.png\n", sep = "")

############################################################
## End 02_goodwinUnionAffiliation.R
############################################################
