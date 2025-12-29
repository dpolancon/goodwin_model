############################################################
## 02_goodwinMinsky.R
## Keen-orthodox Goodwin–Minsky (3D) with Option 1 fragility:
##   i(d) = i0 + i1*d  (endogenous risk premium / stress)
##
## Features:
## - High-employment calibration via Phillips intercept (targets e*)
## - Robust scan in lambda for Hopf functional H = a1*a2 - a3
## - Simulation for {below, near, above} Hopf candidate
## - Low initial debt (history: rising leverage cycles)
## - Adds per-capita output in levels: y_pc(t)
## - Saves figures to ./outputs
############################################################

suppressPackageStartupMessages({
  library(deSolve)
  library(numDeriv)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scatterplot3d)
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
## 1) Model blocks (Keen-orthodox + Option 1 fragility)
## ---------------------------------------------------------
Phi <- function(e, par) par$phi0 + par$phi1 * e

kappa <- function(r, par) {
  par$kmin + (par$kmax - par$kmin) / (1 + exp(-par$lambda * (r - par$r0)))
}

pi_fun <- function(omega) 1 - omega

## Option 1: nonlinear effective interest i(d)
i_eff <- function(d, par) par$i0 + par$i1 * d

r_fun <- function(omega, d, par) {
  (pi_fun(omega) - i_eff(d, par) * d) / par$sigma
}

gY_fun <- function(r, par) kappa(r, par) / par$sigma - par$delta

F_vec <- function(x, par) {
  e     <- x[1]
  omega <- x[2]
  d     <- x[3]
  
  pi  <- pi_fun(omega)
  r   <- r_fun(omega, d, par)
  kap <- kappa(r, par)
  gY  <- kap / par$sigma - par$delta
  g_n <- par$alpha + par$beta
  
  de     <- (gY - g_n) * e
  domega <- (Phi(e, par) - par$alpha) * omega
  
  ## Debt ratio dynamics with nonlinear interest:
  dd <- kap - pi + i_eff(d, par) * d - d * gY
  
  c(de, domega, dd)
}

## ---------------------------------------------------------
## 2) Fixed point solver (interior)
##     - e* from Phillips: Phi(e*) = alpha
##     - then solve for r* from kappa(r*) = sigma*(g_n + delta)
##     - then enforce d* stationarity at gY=g_n and omega* admissible
##
## NOTE:
## With nonlinear i(d), omega* depends on d* through i(d*)d*.
## We solve for d* by rootfinding stationarity in dd at the fixed point.
## ---------------------------------------------------------
e_star_from_phillips <- function(par) (par$alpha - par$phi0) / par$phi1

fixed_point_solver <- function(par, tol = 1e-10) {
  g_n <- par$alpha + par$beta
  
  ## 1) e*
  e_star <- e_star_from_phillips(par)
  if (!is.finite(e_star) || e_star <= 0 || e_star >= 1) {
    stop(sprintf("Inadmissible e*: %.6f (needs 0<e*<1).", e_star))
  }
  
  ## 2) r* from growth condition
  k_target <- par$sigma * (g_n + par$delta)
  kappa_minus <- function(r) kappa(r, par) - k_target
  
  r_lo <- -2; r_hi <- 2
  for (iter in 1:60) {
    if (kappa_minus(r_lo) <= 0 && kappa_minus(r_hi) >= 0) break
    r_lo <- r_lo * 1.5
    r_hi <- r_hi * 1.5
  }
  if (!(kappa_minus(r_lo) <= 0 && kappa_minus(r_hi) >= 0)) {
    stop("Could not bracket r*: kappa(r) never hits target.")
  }
  r_star <- uniroot(kappa_minus, lower = r_lo, upper = r_hi, tol = tol)$root
  
  ## 3) solve for d* using dd=0 at fixed point
  ## At fixed point: gY = g_n and kappa = k_target.
  ## Need omega consistent with r equation:
  ## r_star = (pi - i_eff(d)*d)/sigma  and pi = 1 - omega.
  ## Rearranged: omega(d) = 1 - [sigma*r_star + i_eff(d)*d]
  omega_of_d <- function(d) 1 - (par$sigma * r_star + i_eff(d, par) * d)
  
  dd_at_fp <- function(d) {
    omega <- omega_of_d(d)
    if (!is.finite(omega) || omega <= 0 || omega >= 1) return(NA_real_)
    pi <- 1 - omega
    kap <- k_target
    gY  <- g_n
    ## dd = kap - pi + i(d)*d - d*gY
    kap - pi + i_eff(d, par) * d - d * gY
  }
  
  ## bracket d*
  d_grid <- seq(1e-8, 50, length.out = 2000)
  f_grid <- sapply(d_grid, dd_at_fp)
  ok_idx <- which(is.finite(f_grid))
  if (length(ok_idx) < 10) stop("No feasible omega(d) region for d*; calibration too tight.")
  
  d_grid2 <- d_grid[ok_idx]
  f_grid2 <- f_grid[ok_idx]
  
  sgn <- sign(f_grid2)
  chg <- which(sgn[-1] * sgn[-length(sgn)] < 0)
  
  if (length(chg) == 0) {
    ## No sign change: pick closest-to-zero as "pseudo fp" and flag
    j <- which.min(abs(f_grid2))
    d_star <- d_grid2[j]
    omega_star <- omega_of_d(d_star)
    ## still require positive debt for interior fp
    if (!is.finite(d_star) || d_star <= 0) stop("Inadmissible d* (no root, and closest is nonpositive).")
    ## residuals
    x_star <- c(e_star, omega_star, d_star)
    res <- F_vec(x_star, par)
    return(list(
      e_star=e_star, omega_star=omega_star, d_star=d_star,
      r_star=r_star, kappa_star=kappa(r_star, par),
      gY_star=gY_fun(r_star, par),
      residuals=res,
      note="No dd-root: used closest-to-zero dd(d) on admissible omega(d)."
    ))
  }
  
  ## Use first sign change interval
  k <- chg[1]
  lo <- d_grid2[k]; hi <- d_grid2[k+1]
  d_star <- uniroot(function(d) dd_at_fp(d), lower = lo, upper = hi, tol = tol)$root
  omega_star <- omega_of_d(d_star)
  
  if (!(is.finite(d_star) && d_star > 0)) stop(sprintf("Inadmissible d*: %.6f", d_star))
  if (!(omega_star > 0 && omega_star < 1)) stop(sprintf("Inadmissible omega*: %.6f", omega_star))
  
  x_star <- c(e_star, omega_star, d_star)
  res <- F_vec(x_star, par)
  
  list(
    e_star = e_star,
    omega_star = omega_star,
    d_star = d_star,
    r_star = r_star,
    kappa_star = kappa(r_star, par),
    gY_star = gY_fun(r_star, par),
    residuals = res,
    note = NA_character_
  )
}

## ---------------------------------------------------------
## 3) Local Hopf diagnostic: H = a1*a2 - a3
## ---------------------------------------------------------
routh_hurwitz_from_J <- function(J) {
  trJ <- sum(diag(J))
  a1  <- -trJ
  
  M11 <- det(J[2:3, 2:3])
  M22 <- det(J[c(1,3), c(1,3)])
  M33 <- det(J[1:2, 1:2])
  a2  <- M11 + M22 + M33
  
  a3  <- -det(J)
  
  H <- a1 * a2 - a3
  list(a1=a1, a2=a2, a3=a3, HopfGap=H)
}

local_diagnostics <- function(par) {
  fp <- fixed_point_solver(par)
  x_star <- c(fp$e_star, fp$omega_star, fp$d_star)
  
  J <- jacobian(func = function(x) F_vec(x, par), x = x_star)
  ev <- eigen(J, only.values = TRUE)$values
  rh <- routh_hurwitz_from_J(J)
  
  list(fp=fp, J=J, eigenvalues=ev, rh=rh)
}

HopfGap_safe <- function(par) {
  tmp <- try(local_diagnostics(par), silent = TRUE)
  if (inherits(tmp, "try-error")) return(list(ok=FALSE, H=NA_real_, msg=as.character(tmp)))
  list(ok=TRUE, H=tmp$rh$HopfGap, fp=tmp$fp, ev=tmp$eigenvalues, rh=tmp$rh)
}

## ---------------------------------------------------------
## 4) ODE in transformed states (logit(e), logit(omega), log(d))
## ---------------------------------------------------------
goodwin_minsky_logstates <- function(t, state, par) {
  x <- state["x"]
  y <- state["y"]
  z <- state["z"]
  
  e     <- inv_logit(x)
  omega <- inv_logit(y)
  d     <- exp(z)
  
  vec <- F_vec(c(e, omega, d), par)
  de <- vec[1]; domega <- vec[2]; dd <- vec[3]
  
  dx <- de / (e * (1 - e))
  dy <- domega / (omega * (1 - omega))
  dz <- dd / d
  
  list(c(dx, dy, dz))
}

## ---------------------------------------------------------
## 5) Baseline calibration (high employment target)
## ---------------------------------------------------------
## Target fixed point employment in [0.75, 0.95]
e_star_target <- 0.90

## Phillips slope: keep moderate to avoid labor-market apocalypse
phi1 <- 0.25
phi0 <- 0.02 - phi1 * e_star_target  # ensures e*=(alpha-phi0)/phi1 = e_star_target

par0 <- list(
  alpha = 0.02,
  beta  = 0.015,
  
  sigma = 3,
  delta = 0.015,
  
  ## Option 1 interest schedule
  i0 = 0.04,     # baseline interest
  i1 = 0.0025,    # fragility (tune this)
  
  ## Phillips
  phi0 = phi0,
  phi1 = phi1,
  
  ## Investment (logistic)
  kmin   = 0.00,
  kmax   = 0.40,
  lambda = 40,   # will be scanned
  r0     = 0.05
)

cat("\n--- Baseline fixed point + eigen diagnostics ---\n")
base <- HopfGap_safe(par0)
if (!base$ok) stop(base$msg)
print(base$fp)
print(base$rh)
cat("Eigenvalues:\n"); print(base$ev)

## ---------------------------------------------------------
## 6) Robust scan in lambda (admissibility + H(lambda))
## ---------------------------------------------------------
lam_grid <- seq(10, 160, by = 1)

scan <- lapply(lam_grid, function(lam) {
  par <- par0
  par$lambda <- lam
  out <- HopfGap_safe(par)
  tibble(
    lambda = lam,
    ok = out$ok,
    H = out$H,
    maxRe = if (out$ok) max(Re(out$ev)) else NA_real_,
    msg = ifelse(out$ok, NA_character_, out$msg),
    e_star = if (out$ok) out$fp$e_star else NA_real_,
    omega_star = if (out$ok) out$fp$omega_star else NA_real_,
    d_star = if (out$ok) out$fp$d_star else NA_real_
  )
}) %>% bind_rows()

cat("\n--- Scan summary ---\n")
print(scan %>% count(ok))

p_ok <- ggplot(scan, aes(x=lambda, y=as.numeric(ok))) +
  geom_line(linewidth=0.6) +
  scale_y_continuous(breaks=c(0,1), labels=c("inadmissible","admissible")) +
  labs(title="Admissible interior fixed point across λ (with i(d)=i0+i1 d)",
       x="Investment aggressiveness (λ)", y="") +
  theme_minimal()

ggsave("outputs/fig02_ok_lambda.png", p_ok, width=9, height=3.5, dpi=300)

p_H <- ggplot(scan %>% filter(ok), aes(lambda, H)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.4) +
  geom_line(linewidth=0.7) +
  labs(title="Hopf functional H(λ) on admissible region",
       subtitle="H = a1*a2 - a3 (local Hopf boundary if it crosses 0)",
       x="λ", y="H") +
  theme_minimal()

ggsave("outputs/fig02_H_lambda.png", p_H, width=9, height=4.5, dpi=300)

## ---------------------------------------------------------
## 7) Locate Hopf bracket in lambda (within admissible set)
## ---------------------------------------------------------
scan_ok <- scan %>% filter(ok) %>% arrange(lambda)

idx <- which(sign(scan_ok$H[-1]) * sign(scan_ok$H[-nrow(scan_ok)]) < 0)

if (length(idx) == 0) {
  cat("\nNo sign change in H(λ) over admissible region.\n")
  cat("You can still simulate around a 'near-Hopf' point by minimizing |H|.\n")
  j <- which.min(abs(scan_ok$H))
  lambda0 <- scan_ok$lambda[j]
  lambda_set <- c(max(1, lambda0 - 5), lambda0, lambda0 + 5)
  names(lambda_set) <- c("below","near","above")
} else {
  k <- idx[1]
  lo <- scan_ok$lambda[k]
  hi <- scan_ok$lambda[k+1]
  
  H_of_lambda <- function(lam) {
    par <- par0
    par$lambda <- lam
    HopfGap_safe(par)$H
  }
  
  lambda0 <- uniroot(H_of_lambda, lower=lo, upper=hi)$root
  cat(sprintf("\nHopf bracket found: [%.2f, %.2f]\n", lo, hi))
  cat(sprintf("lambda_Hopf ≈ %.6f\n", lambda0))
  
  lambda_set <- c(max(1, lambda0 - 5), lambda0, lambda0 + 5)
  names(lambda_set) <- c("below","near","above")
}

cat("\n--- Lambda set to simulate ---\n")
print(lambda_set)

## ---------------------------------------------------------
## 8) Initial conditions: high employment, good omega, LOW debt
## ---------------------------------------------------------
fp_base <- fixed_point_solver(par0)

## Keep initial e within your desired range:
e0     <- clamp01(min(0.95, max(0.75, fp_base$e_star + 0.01)))
omega0 <- clamp01(fp_base$omega_star - 0.01)

## Low debt initial condition:
d0 <- 0.05  # small positive (log(d0) exists); tells "leverage builds"

state0 <- c(x = logit(e0), y = logit(omega0), z = log(d0))

times <- seq(0, 500, by = 0.05)

run_one <- function(lam, tag) {
  par <- par0
  par$lambda <- lam
  
  sol <- ode(
    y = state0,
    times = times,
    func = goodwin_minsky_logstates,
    parms = par,
    method = "lsoda",
    rtol = 1e-8,
    atol = 1e-10
  )
  
  df <- as.data.frame(sol) %>%
    mutate(
      e     = inv_logit(x),
      omega = inv_logit(y),
      d     = exp(z),
      
      pi    = 1 - omega,
      i_eff = i_eff(d, par),
      r     = (pi - i_eff * d) / par$sigma,
      kappa = kappa(r, par),
      
      gY    = kappa / par$sigma - par$delta,
      gy    = gY - par$beta,                 # per-capita growth rate
      
      dt    = c(0, diff(time)),
      y_pc  = exp(cumsum(gy * dt)),          # per-capita output index, y_pc(0)=1
      
      lambda = lam,
      tag = tag
    ) %>%
    select(time, tag, lambda,
           e, omega, d,
           pi, i_eff, r, kappa, gY, gy, y_pc)
  
  df
}

df_all <- bind_rows(
  run_one(lambda_set["below"], "below"),
  run_one(lambda_set["near"],  "near"),
  run_one(lambda_set["above"], "above")
)

## ---------------------------------------------------------
## 9) Plots
## ---------------------------------------------------------
df_core_long <- df_all %>%
  select(time, tag, e, omega, d) %>%
  pivot_longer(cols=c(e, omega, d), names_to="var", values_to="value") %>%
  mutate(var = factor(var, levels=c("e","omega","d"),
                      labels=c("Employment rate (e)", "Wage share (ω)", "Debt ratio (d)")))

p_core <- ggplot(df_core_long, aes(time, value)) +
  geom_line(linewidth=0.6) +
  facet_grid(var ~ tag, scales="free_y") +
  labs(title="Goodwin–Minsky–Keen dynamics around Hopf candidate (core states)",
       subtitle="Option 1 fragility: i(d)=i0+i1 d (financial stress can destabilize debt)",
       x="time", y="") +
  theme_minimal() +
  theme(legend.position="none", panel.spacing=unit(1.0,"lines"))

ggsave("outputs/fig02_timepaths_core.png", p_core, width=12, height=7, dpi=300)

df_fin_long <- df_all %>%
  select(time, tag, r, pi, i_eff, kappa, gY, y_pc) %>%
  pivot_longer(cols=c(r, pi, i_eff, kappa, gY, y_pc), names_to="var", values_to="value") %>%
  mutate(var = factor(var, levels=c("r","pi","i_eff","kappa","gY","y_pc"),
                      labels=c("Net profit rate r",
                               "Profit share π",
                               "Effective interest i(d)",
                               "Investment κ(r)",
                               "Output growth gY",
                               "Output per capita (index)")))

p_fin <- ggplot(df_fin_long, aes(time, value)) +
  geom_line(linewidth=0.6) +
  facet_grid(var ~ tag, scales="free_y") +
  labs(title="Goodwin–Minsky–Keen dynamics (finance block)",
       x="time", y="") +
  theme_minimal() +
  theme(legend.position="none", panel.spacing=unit(1.0,"lines"))

ggsave("outputs/fig02_timepaths_finance.png", p_fin, width=12, height=10, dpi=300)

## 3D trajectory (above)
df3 <- df_all %>% filter(tag=="above")
png("outputs/fig02_3D_e_omega_d_above.png", width=1400, height=900, res=150)
scatterplot3d(df3$e, df3$omega, df3$d,
              type="l", angle=35,
              xlab="e", ylab="ω", zlab="d",
              main="3D state-space trajectory (above Hopf candidate)")
dev.off()

## ---------------------------------------------------------
## 10) Quick regime summary
## ---------------------------------------------------------
summ <- df_all %>%
  group_by(tag) %>%
  summarise(
    e_min = min(e), e_max = max(e),
    omega_min = min(omega), omega_max = max(omega),
    d_min = min(d), d_max = max(d),
    iEff_min = min(i_eff), iEff_max = max(i_eff),
    .groups="drop"
  )

cat("\n--- Range summary by regime tag ---\n")
print(summ)

cat("\nSaved outputs in ./outputs:\n",
    "- fig02_ok_lambda.png\n",
    "- fig02_H_lambda.png\n",
    "- fig02_timepaths_core.png\n",
    "- fig02_timepaths_finance.png\n",
    "- fig02_3D_e_omega_d_above.png\n", sep = "")

############################################################
## End 02_goodwinMinsky.R
############################################################
