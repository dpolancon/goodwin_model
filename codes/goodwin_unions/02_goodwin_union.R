############################################################
## 02_goodwinUnionAffiliation_v2.R
## Goodwin (2D) + Union Affiliation as bifurcation parameter
##
## Key modifications vs v1:
## (i) Radial saturation (normal-form style) to guarantee
##      bounded limit cycle when mu>0.
## (ii) "near" set to s0 + eps (not exactly s0) to avoid
##      conservative center + numerical annulus.
## (iii) Phase portrait plotted AFTER burn-in to show thin
##      asymptotic orbit (limit cycle) instead of transient.
##
## System (original variables; deviations auxiliary):
##   u = e - e*
##   v = ω - ω*
##   r^2 = u^2 + v^2
##   de/dt     = -a*v - c*u*r^2
##   dω/dt     =  b*u + μ(s)*v - c*v*r^2
##   μ(s) = ζ(s - s0)
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
  
  ## deviations
  u  <- e     - par$e_star
  v  <- omega - par$omega_star
  r2 <- u^2 + v^2
  
  ## Normal-form style radial saturation (bounded cycles for mu>0)
  de     <- -par$a * v - par$c * u * r2
  domega <-  par$b * u + mu * v - par$c * v * r2
  
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
  
  ## logit derivatives
  dx <- de     / (e * (1 - e))
  dy <- domega / (omega * (1 - omega))
  
  list(c(dx, dy))
}

## ---------------------------------------------------------
## 3) Calibration
## ---------------------------------------------------------
par0 <- list(
  ## Fixed point location
  e_star     = 0.90,
  omega_star = 0.65,
  
  ## Goodwin feedback
  a = 0.20,
  b = 0.25,
  
  ## Radial saturation strength
  c = 5.0,
  
  ## Institutional flip
  s0   = 0.20,
  zeta = 1.50,
  
  ## will be set per run
  s_fixed = NA_real_
)

cat("\n--- Diagnostics at s=s0 (Hopf boundary in linearization: tr=0, det=ab>0) ---\n")
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

ggsave("outputs/v2_mu_scan_s.png", p_mu, width=9, height=4, dpi=300)

p_stab <- ggplot(scan, aes(s, maxRe)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  labs(title="Local stability proxy: max Re(eigenvalue) at (e*, ω*)",
       subtitle="maxRe<0 stable, maxRe>0 unstable; crossing at s=s0",
       x="Unionization s", y="max Re(λ)") +
  theme_minimal()

ggsave("outputs/v2_stability_scan_s.png", p_stab, width=9, height=4, dpi=300)

## ---------------------------------------------------------
## 5) Simulate {below, near, above} with improved 'near'
## ---------------------------------------------------------
eps_near <- 0.005  # small >0 so near is slightly supercritical
ds_far   <- 0.05

s_set <- c(
  below = max(0.01, par0$s0 - ds_far),
  near  = min(0.99, par0$s0 + eps_near),
  above = min(0.99, par0$s0 + ds_far)
)

cat("\n--- s values to simulate ---\n")
print(s_set)

## Initial conditions (small perturbation)
e0     <- clamp01(par0$e_star + 0.01)
omega0 <- clamp01(par0$omega_star - 0.01)
state0 <- c(x = logit(e0), y = logit(omega0))

times <- seq(0, 800, by = 0.05)
burnin <- 200  # for phase portrait thinning (plot asymptotic orbit only)

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
  labs(title="Goodwin with Union Affiliation: time paths around s0 (v2)",
       subtitle=paste0("below: μ<0 damped; near: μ>0 small cycle; above: μ>0 larger cycle.  burn-in=", burnin, " for phase plot"),
       x="time", y="") +
  theme_minimal() +
  theme(legend.position="none")

ggsave("outputs/v2_timepaths_e_omega.png", p_time, width=12, height=6, dpi=300)

## Phase portrait AFTER burn-in (thin asymptotic orbit)
df_phase <- df_all %>% filter(time >= burnin)

p_phase <- ggplot(df_phase, aes(e, omega)) +
  geom_path(linewidth=0.5) +
  facet_wrap(~tag, scales="free") +
  labs(title=paste0("Phase portrait after burn-in (t ≥ ", burnin, "): asymptotic orbit"),
       x="e", y="ω") +
  theme_minimal()

ggsave("outputs/v2_phase_e_omega_burnin.png", p_phase, width=10, height=4, dpi=300)

## Summary ranges (full sample)
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
    "- v2_mu_scan_s.png\n",
    "- v2_stability_scan_s.png\n",
    "- v2_timepaths_e_omega.png\n",
    "- v2_phase_e_omega_burnin.png\n", sep="")

#####################################################
## End 02_goodwinUnionAffiliation_v2.R
#####################################################
