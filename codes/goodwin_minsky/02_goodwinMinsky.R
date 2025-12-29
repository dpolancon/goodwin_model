############################################################
## 02_simulate_near_hopf_robust.R
## Robust Hopf scan in lambda + simulation around Hopf
## (handles inadmissible fixed points by skipping them)
############################################################

suppressPackageStartupMessages({
  library(deSolve)
  library(numDeriv)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scatterplot3d)
  library(tibble) 
  library(kableExtra)
  library(knitr)
  library(fs)        # dir_create(), path_ext()
})

## -----------------------------
## 0) Utility: logit transforms
## -----------------------------
inv_logit <- function(x) 1 / (1 + exp(-x))
logit     <- function(p) log(p / (1 - p))
clamp01   <- function(u, eps = 1e-10) pmin(pmax(u, eps), 1 - eps)

## -----------------------------
## 1) Model blocks (Keen-orthodox)
## -----------------------------
Phi <- function(e, par) par$phi0 + par$phi1 * e

kappa <- function(r, par) {
  par$kmin + (par$kmax - par$kmin) / (1 + exp(-par$lambda * (r - par$r0)))
}

pi_fun <- function(omega) 1 - omega
r_fun  <- function(omega, d, par) (pi_fun(omega) - par$i * d) / par$sigma
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
  dd     <- kap - pi + par$i * d - d * gY
  
  c(de, domega, dd)
}

## -----------------------------
## 2) Fixed point solver (interior)
## -----------------------------
e_star_from_phillips <- function(par) (par$alpha - par$phi0) / par$phi1

fixed_point_solver <- function(par, tol = 1e-10) {
  g_n <- par$alpha + par$beta
  
  e_star <- e_star_from_phillips(par)
  if (!is.finite(e_star) || e_star <= 0 || e_star >= 1) {
    stop(sprintf("Inadmissible e*: %.6f", e_star))
  }
  
  k_target <- par$sigma * (g_n + par$delta)
  
  kappa_minus <- function(r) kappa(r, par) - k_target
  
  r_lo <- -2; r_hi <- 2
  for (iter in 1:60) {
    if (kappa_minus(r_lo) <= 0 && kappa_minus(r_hi) >= 0) break
    r_lo <- r_lo * 1.5
    r_hi <- r_hi * 1.5
  }
  if (!(kappa_minus(r_lo) <= 0 && kappa_minus(r_hi) >= 0)) {
    stop("Could not bracket r* (kappa never hits target).")
  }
  
  r_star <- uniroot(kappa_minus, lower = r_lo, upper = r_hi, tol = tol)$root
  
  if (g_n <= 0) stop("Need g_n > 0.")
  
  d_star <- (k_target - par$sigma * r_star) / g_n
  omega_star <- 1 - par$i * d_star - par$sigma * r_star
  
  if (!(is.finite(d_star) && d_star > 0)) stop(sprintf("Inadmissible d*: %.6f", d_star))
  if (!(omega_star > 0 && omega_star < 1)) stop(sprintf("Inadmissible omega*: %.6f", omega_star))
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
    residuals = res
  )
}

## -----------------------------
## 3) Local Hopf diagnostic: H = a1*a2 - a3
## -----------------------------
routh_hurwitz_from_J <- function(J) {
  trJ <- sum(diag(J))
  a1  <- -trJ
  
  M11 <- det(J[2:3, 2:3])
  M22 <- det(J[c(1,3), c(1,3)])
  M33 <- det(J[1:2, 1:2])
  a2  <- M11 + M22 + M33
  
  a3  <- -det(J)
  
  H <- a1 * a2 - a3
  list(a1 = a1, a2 = a2, a3 = a3, HopfGap = H)
}

local_diagnostics <- function(par) {
  fp <- fixed_point_solver(par)
  x_star <- c(fp$e_star, fp$omega_star, fp$d_star)
  
  J <- jacobian(func = function(x) F_vec(x, par), x = x_star)
  ev <- eigen(J, only.values = TRUE)$values
  rh <- routh_hurwitz_from_J(J)
  
  list(fp = fp, J = J, eigenvalues = ev, rh = rh)
}

HopfGap_safe <- function(par) {
  tmp <- try(local_diagnostics(par), silent = TRUE)
  if (inherits(tmp, "try-error")) return(list(ok = FALSE, H = NA_real_, msg = as.character(tmp)))
  list(ok = TRUE, H = tmp$rh$HopfGap, fp = tmp$fp, ev = tmp$eigenvalues)
}

## -----------------------------
## 4) ODE in transformed states (x,y,z)
## -----------------------------
goodwin_minsky_keen_logstates <- function(t, state, par) {
  x <- state["x"]
  y <- state["y"]
  z <- state["z"]
  
  e     <- inv_logit(x)
  omega <- inv_logit(y)
  d     <- exp(z)
  
  vec <- F_vec(c(e, omega, d), par)
  de <- vec[1]; domega <- vec[2]; dd <- vec[3]
  
  dx <- de     / (e * (1 - e))
  dy <- domega / (omega * (1 - omega))
  dz <- dd / d
  
  list(c(dx, dy, dz))
}

## -----------------------------
## 5) Baseline calibration (feasible high-employment regime)
## -----------------------------
e_star_target <- 0.92     # puts cycles in the 0.75–0.95 neighborhood
phi1 <- 0.25
phi0 <- 0.02 - phi1 * e_star_target

par0 <- list(
  alpha = 0.02,
  beta  = 0.015,
  sigma = 3,
  delta = 0.015,
  i     = 0.05,
  
  ## Phillips (disciplined to target e*)
  phi1  = phi1,                       # try 0.05–0.08
  phi0  = phi0,
  
  ## Investment function (instability lever)
  kmin   = 0.00,
  kmax   = 0.60,
  lambda = 50,
  r0     = 0.05
)


## -----------------------------
## 6) Robust scan in lambda
## -----------------------------
lam_grid <- seq(10, 140, by = 1)

scan <- lapply(lam_grid, function(lam) {
  par <- par0
  par$lambda <- lam
  out <- HopfGap_safe(par)
  
  tibble(
    lambda = lam,
    ok = out$ok,
    H = out$H,
    msg = ifelse(out$ok, NA_character_, out$msg),
    e_star = if (out$ok) out$fp$e_star else NA_real_,
    omega_star = if (out$ok) out$fp$omega_star else NA_real_,
    d_star = if (out$ok) out$fp$d_star else NA_real_
  )
}) %>% bind_rows()

cat("\n--- Scan summary ---\n")
print(scan %>% count(ok))

## Plot admissible region + H(lambda)
dir.create("outputs", showWarnings = FALSE)

p_ok <- ggplot(scan, aes(x = lambda, y = as.numeric(ok))) +
  geom_line(linewidth = 0.6) +
  scale_y_continuous(breaks = c(0,1), labels = c("inadmissible", "admissible")) +
  labs(title = "Admissible interior fixed point across λ",
       x = "Investment aggressiveness (λ)", y = "") +
  theme_minimal()

ggsave("outputs/fig02_admissible_lambda.png", p_ok, width = 9, height = 3.5, dpi = 300)

p_H <- ggplot(scan %>% filter(ok), aes(lambda, H)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  labs(title = "Hopf functional H(λ) on admissible region",
       subtitle = "H = a1*a2 - a3 (local Hopf boundary if H crosses 0)",
       x = "λ", y = "H") +
  theme_minimal()

ggsave("outputs/fig02_H_lambda.png", p_H, width = 9, height = 4.5, dpi = 300)

## -----------------------------
## 7) Find a Hopf bracket inside admissible region
## -----------------------------
scan_ok <- scan %>% filter(ok) %>% arrange(lambda)

## locate sign change in H
idx <- which(sign(scan_ok$H[-1]) * sign(scan_ok$H[-nrow(scan_ok)]) < 0)

if (length(idx) == 0) {
  cat("\nNo sign change in H(λ) on admissible region.\n")
  cat("Meaning: with this calibration, you do not cross a local Hopf by moving λ alone.\n")
  cat("Try shifting i, r0, kmax, sigma, or phi0 to move the boundary.\n")
  
  lambda0 <- par0$lambda
  lambda_set <- c(lambda0 - 5, lambda0, lambda0 + 5)
  names(lambda_set) <- c("below", "near", "above")
} else {
  k <- idx[1]
  lo <- scan_ok$lambda[k]
  hi <- scan_ok$lambda[k + 1]
  
  H_of_lambda <- function(lam) {
    par <- par0
    par$lambda <- lam
    HopfGap_safe(par)$H
  }
  
  lambda0 <- uniroot(H_of_lambda, lower = lo, upper = hi)$root
  cat(sprintf("\nHopf bracket found: [%.2f, %.2f]\n", lo, hi))
  cat(sprintf("lambda_Hopf ≈ %.6f\n", lambda0))
  
  lambda_set <- c(max(1, lambda0 - 5), lambda0, lambda0 + 5)
  names(lambda_set) <- c("above", "near", "below")
}

cat("\n--- Lambda set to simulate ---\n")
print(lambda_set)

## -----------------------------
## -----------------------------
## 8) Initial conditions (low debt, high employment, feasible ω)
## -----------------------------
fp_base <- fixed_point_solver(par0)

e0     <- clamp01(fp_base$e_star + 0.01)      # small displacement
omega0 <- clamp01(fp_base$omega_star - 0.01)  # small displacement
d0     <- 0                                # low initial leverage (history builds)

state0 <- c(x = logit(e0), y = logit(omega0), z = log(d0))


if (inherits(fp_base, "try-error")) {
  base_row <- scan_ok %>% slice(1)
  par_tmp <- par0; par_tmp$lambda <- base_row$lambda
  fp_base <- fixed_point_solver(par_tmp)
  cat(sprintf("\nBaseline fp inadmissible. Using λ=%.2f for initial point.\n", base_row$lambda))
} else {
  par_tmp <- par0
}

e0     <- clamp01(fp_base$e_star + 0.02)
omega0 <- clamp01(fp_base$omega_star - 0.02)
d0     <- max(1e-6, fp_base$d_star * 1.05)

state0 <- c(x = logit(e0), y = logit(omega0), z = log(d0))
times <- seq(0, 500, by = 0.05)

run_one <- function(lam, tag) {
  par <- par0
  par$lambda <- lam
  
  sol <- ode(
    y = state0,
    times = times,
    func = goodwin_minsky_keen_logstates,
    parms = par,
    method = "lsoda",
    rtol = 1e-8,
    atol = 1e-10
  )
  
  as.data.frame(sol) |>
    mutate(
      e     = inv_logit(x),
      omega = inv_logit(y),
      d     = exp(z),
      
      pi    = 1 - omega,
      r     = (pi - par$i * d) / par$sigma,
      kappa = kappa(r, par),
      
      gY    = kappa / par$sigma - par$delta,
      gy    = gY - par$beta,
      
      dt    = c(0, diff(time)),
      y_pc  = exp(cumsum(gy * dt)),   # per-capita output index (y_pc(0) = 1)
      
      lambda = lam,
      tag = tag
    ) |>
    select(time, tag, lambda,
           e, omega, d,
           r, pi, kappa, gY, gy, y_pc)
}


df_all <- bind_rows(
  run_one(lambda_set["above"], "above"),
  run_one(lambda_set["near"],  "near"),
  run_one(lambda_set["below"], "below"))

## Stacked time paths: core
df_core_long <- df_all |>
  select(time, tag, e, omega, d) |>
  pivot_longer(
    cols = c(e, omega, d),
    names_to = "var",
    values_to = "value"
  ) |>
  mutate(
    var = factor(
      var,
      levels = c("e", "omega", "d"),
      labels = c(
        "Employment rate (e)",
        "Wage share (ω)",
        "Debt ratio (d)")))

p_core <- ggplot(df_core_long, aes(time, value)) +
  geom_line(linewidth = 0.6) +
  facet_grid(var ~ tag, scales = "free_y") +
  labs(title = "Keen-orthodox dynamics around Hopf candidate (core states)",
       x = "time", y = "") +
  theme_minimal() +
  theme(legend.position = "none", panel.spacing = unit(1.0, "lines"))

ggsave("outputs/fig02_timepaths_core.png", p_core, width = 12, height = 7, dpi = 300)

## Stacked time paths: finance
df_fin_long <- df_all |>
  select(time, tag, r, pi, kappa, gY, y_pc) |>
  pivot_longer(
    cols = c(r, pi, kappa, gY, y_pc),
    names_to = "var",
    values_to = "value"
  ) |>  mutate(var = factor(var,
            levels = c("r", "pi", "kappa", "gY", "y_pc"),
            labels = c(
              "Net profit rate r",
              "Profit share π",
              "Investment κ(r)",
              "Output growth rate gY",
              "Output per capita")))


p_fin <- ggplot(df_fin_long, aes(time, value)) +
  geom_line(linewidth = 0.6) +
  facet_grid(var ~ tag, scales = "free_y") +
  labs(title = "Keen-orthodox dynamics around Hopf candidate (finance block)",
       x = "time", y = "") +
  theme_minimal() +
  theme(legend.position = "none", panel.spacing = unit(1.0, "lines"))

ggsave("outputs/fig02_timepaths_finance.png", p_fin, width = 12, height = 9, dpi = 300)

## 3D trajectory (above)
df3 <- df_all %>% filter(tag == "above")
png("outputs/fig02_3D_e_omega_d_above.png", width = 1400, height = 900, res = 150)
scatterplot3d(df3$e, df3$omega, df3$d,
              type = "l", angle = 35,
              xlab = "e", ylab = "ω", zlab = "d",
              main = "3D state-space trajectory (above Hopf candidate)")
dev.off()

cat("\nSaved outputs in ./outputs:\n",
    "- fig02_admissible_lambda.png\n",
    "- fig02_H_lambda.png\n",
    "- fig02_timepaths_core.png\n",
    "- fig02_timepaths_finance.png\n",
    "- fig02_3D_e_omega_d_above.png\n", sep = "")
############################################################


############################################################
## Helper: table_as_is()  (patched)
############################################################



calibration_tbl <- tribble(
  ~parameter, ~symbol, ~value, ~description,
  "Productivity growth",        "α", 0.02,  "Labor productivity growth rate",
  "Population growth",          "β", 0.01,  "Labor force growth rate",
  "Capital-output ratio",       "σ", 3.00,  "Inverse of output-capital ratio",
  "Depreciation",               "δ", 0.00,  "Capital depreciation rate",
  "Interest rate",              "i", 0.05,  "Interest on outstanding debt",
  
  "Phillips slope",             "φ₁", 0.06, "Sensitivity of wage inflation to employment",
  "Phillips intercept",         "φ₀", -0.0328,
  "Phillips curve intercept calibrated to e* = 0.88",
  
  "Investment floor",           "κ_min", 0.00, "Lower bound on investment rate",
  "Investment ceiling",         "κ_max", 0.60, "Upper bound on investment rate",
  "Investment aggressiveness",  "λ", 60,    "Slope of logistic investment function",
  "Target profit rate",         "r₀", 0.05, "Center of investment response"
)

table_as_is <- function(
    x,
    file,
    caption = NULL,
    label = NULL,
    format = NULL
) {
  stopifnot(is.data.frame(x))
  
  # Ensure output directory exists
  fs::dir_create(fs::path_dir(file))
  
  # Infer format from file extension if not supplied
  if (is.null(format)) {
    ext <- fs::path_ext(file)
    format <- switch(
      ext,
      "tex" = "latex",
      "html" = "html",
      "md" = "markdown",
      "qmd" = "markdown",
      "pdf" = "latex",
      stop("Cannot infer table format from file extension.")
    )
  }
  
  # Build table strictly as-is
  tab <- knitr::kable(
    x,
    format = format,
    caption = caption,
    booktabs = (format == "latex"),
    longtable = (format == "latex"),
    align = "l"
  )
  
  # Optional label (LaTeX / HTML only)
  if (!is.null(label) && format %in% c("latex", "html")) {
    tab <- kableExtra::kable_styling(tab)
    attr(tab, "label") <- label
  }
  
  # Write to disk
  writeLines(as.character(tab), con = file)
  
  invisible(tab)
}


#### Export table #####

table_as_is(
  calibration_tbl,
  file = "outputs/calibration_table.tex",
  caption = "Baseline calibration parameters",
  label = "tab:calibration"
)

