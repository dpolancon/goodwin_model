############################################################
## 02_goodwinUnionAffiliation_v3.R
## Goodwin (2D) + Union Affiliation as bifurcation parameter
##
## Fixes vs v2:
##  (1) Uses the CORRECT weighted radius for Hopf normal form:
##        R^2 = b*u^2 + a*v^2   (not u^2 + v^2)
##      where u=e-e*, v=ω-ω*.
##  (2) Phase portrait uses LAST window of time (not burn-in cutoff),
##      so you see a thin asymptotic orbit even when convergence is slow.
##  (3) Cleaner regime separation: "below" and "above" farther from s0,
##      "near" slightly above s0 but not microscopic.
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
## 1) Model blocks
## ---------------------------------------------------------
mu_fun <- function(s, par) par$zeta * (s - par$s0)

F_vec <- function(x, par) {
  e     <- x[1]
  omega <- x[2]
  s     <- par$s_fixed
  
  mu <- mu_fun(s, par)
  
  ## deviations around fixed point
  u <- e     - par$e_star
  v <- omega - par$omega_star
  
  ## IMPORTANT: weighted radius consistent with anisotropic linear part
  ## In scaled coords: U = sqrt(b)u, V = sqrt(a)v, so R^2 = U^2+V^2 = b u^2 + a v^2
  R2 <- par$b * u^2 + par$a * v^2
  
  ## Hopf normal-form style saturation (supercritical if c>0)
  de     <- -par$a * v - par$c * u * R2
  domega <-  par$b * u + mu * v - par$c * v * R2
  
  c(de, domega)
}

fixed_point_solver <- function(par) {
  list(
    e_star     = par$e_star,
    omega_star = par$omega_star,
    s_fixed    = par$s_fixed,
    mu         = mu_fun(par$s_fixed, par)
  )
}

local_diagnostics <- function(par) {
  fp <- fixed_point_solver(par)
  x_star <- c(fp$e_star, fp$omega_star)
  
  J <- jacobian(func = function(x) F_vec(x, par), x = x_star)
  ev <- eigen(J, only.values = TRUE)$values
  
  list(
    fp = fp,
    J  = J,
    eigenvalues = ev,
    tr  = sum(diag(J)),
    det = det(J),
    maxRe = max(Re(ev))
  )
}

## ---------------------------------------------------------
## 2) ODE in transformed states (logit(e), logit(omega))
## ---------------------------------------------------------
goodwin_union_logstates <- function(t, state, par) {
  x <- state["x"]
  y <- state["y"]
  
  e     <- inv_logit(x)
  omega <- inv_logit(y)
  
  vec <- F_vec(c(e, omega), par)
  de <- vec[1]; domega <- vec[2]
  
  dx <- de     / (e * (1 - e))
  dy <- domega / (omega * (1 - omega))
  
  list(c(dx, dy))
}

## ---------------------------------------------------------
## 3) Calibration
## ---------------------------------------------------------
par0 <- list(
  e_star     = 0.90,
  omega_star = 0.65,
  
  a = 0.20,
  b = 0.25,
  
  ## saturation strength (bigger => smaller cycle)
  c = 25.0,
  
  s0   = 0.15,
  zeta = 1.50,
  
  s_fixed = NA_real_
)

cat("\n--- Diagnostics at s=s0 ---\n")
par_tmp <- par0; par_tmp$s_fixed <- par0$s0
diag0 <- local_diagnostics(par_tmp)
print(diag0$fp)
cat("Trace:", diag0$tr, "Det:", diag0$det, "maxRe:", diag0$maxRe, "\n")
cat("Eigenvalues:\n"); print(diag0$eigenvalues)

## ---------------------------------------------------------
## 4) Scan in s (stability + mu)
## ---------------------------------------------------------
s_grid <- seq(0.01, 0.99, by = 0.005)

scan <- lapply(s_grid, function(sval) {
  par <- par0; par$s_fixed <- sval
  out <- local_diagnostics(par)
  tibble(
    s = sval,
    mu = out$fp$mu,
    tr = out$tr,
    det = out$det,
    maxRe = out$maxRe,
    complex = any(abs(Im(out$eigenvalues)) > 1e-10)
  )
}) %>% bind_rows()

p_mu <- ggplot(scan, aes(s, mu)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  labs(title="Institutional flip term μ(s) = ζ(s − s0)",
       subtitle="Hopf boundary (linear): μ=0 ⇔ s=s0",
       x="Unionization s", y="μ(s)") +
  theme_minimal()

ggsave("outputs/v3_mu_scan_s.png", p_mu, width=9, height=4, dpi=300)

p_stab <- ggplot(scan, aes(s, maxRe)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  labs(title="Local stability proxy: max Re(eigenvalue) at (e*, ω*)",
       subtitle="Crossing at s=s0; kink where μ^2 = 4ab (focus→node)",
       x="Unionization s", y="max Re(λ)") +
  theme_minimal()

ggsave("outputs/v3_stability_scan_s.png", p_stab, width=9, height=4, dpi=300)

## ---------------------------------------------------------
## 5) Simulate {below, near, above}
## ---------------------------------------------------------
ds_far   <- 0.20     # farther from s0 so regimes are visually clean
eps_near <- 0.02     # not microscopic; avoids super-slow convergence

s_set <- c(
  below = max(0.01, par0$s0 - ds_far),
  near  = min(0.99, par0$s0 + eps_near),
  above = min(0.99, par0$s0 + ds_far)
)

cat("\n--- s values to simulate ---\n")
print(s_set)
cat("Corresponding mu:\n")
print(sapply(s_set, function(ss) mu_fun(ss, par0)))

## Initial conditions
e0     <- clamp01(par0$e_star + 0.01)
omega0 <- clamp01(par0$omega_star - 0.01)
state0 <- c(x = logit(e0), y = logit(omega0))

times <- seq(0, 1200, by = 0.05)

## Phase window: plot only the last T_last units to get thin orbit
T_last <- 150
t_end  <- max(times)
t_min_phase <- t_end - T_last

run_one <- function(sval, tag) {
  par <- par0
  par$s_fixed <- sval
  
  sol <- ode(
    y = state0,
    times = times,
    func = goodwin_union_logstates,
    parms = par,
    method = "lsoda",
    rtol = 1e-10,
    atol = 1e-12
  )
  
  as.data.frame(sol) %>%
    mutate(
      e     = inv_logit(x),
      omega = inv_logit(y),
      s     = sval,
      mu    = mu_fun(sval, par),
      tag   = tag
    ) %>%
    select(time, tag, s, mu, e, omega)
}

df_all <- bind_rows(
  run_one(s_set["below"], "below"),
  run_one(s_set["near"],  "near"),
  run_one(s_set["above"], "above")
)

## ---------------------------------------------------------
## 6) Plots
## ---------------------------------------------------------
## Time paths
df_long <- df_all %>%
  pivot_longer(cols=c(e, omega), names_to="var", values_to="value") %>%
  mutate(var = factor(var, levels=c("e","omega"),
                      labels=c("Employment rate (e)", "Wage share (ω)")))

p_time <- ggplot(df_long, aes(time, value)) +
  geom_line(linewidth=0.6) +
  facet_grid(var ~ tag, scales="free_y") +
  labs(title="Goodwin with Union Affiliation: time paths around s0 (v3)",
       subtitle=paste0("Weighted radius saturation. Phase plot uses last ", T_last, " time units."),
       x="time", y="") +
  theme_minimal() +
  theme(legend.position="none")

ggsave("outputs/v3_timepaths_e_omega.png", p_time, width=12, height=6, dpi=300)

## Phase portrait: last window only
df_phase <- df_all %>% filter(time >= t_min_phase)

p_phase <- ggplot(df_phase, aes(e, omega)) +
  geom_path(linewidth=0.6) +
  facet_wrap(~tag, scales="free") +
  labs(title=paste0("Phase portrait (last ", T_last, " time units): asymptotic orbit"),
       x="e", y="ω") +
  theme_minimal()

ggsave("outputs/v3_phase_e_omega_lastwindow.png", p_phase, width=10, height=4, dpi=300)

## Summary ranges
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
    "- v3_mu_scan_s.png\n",
    "- v3_stability_scan_s.png\n",
    "- v3_timepaths_e_omega.png\n",
    "- v3_phase_e_omega_lastwindow.png\n", sep="")

############################################################
## End 02_goodwinUnionAffiliation_v3.R
############################################################
