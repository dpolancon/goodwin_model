# ============================================================
# R/hooks/simulate_model.R
# Hook for Stage 5
# Must define: simulate_model(row, rF_sim, t_end, dt, x0=NULL)
# Returns a tibble with columns: time, e, omega, d
# ============================================================

simulate_model <- function(row, rF_sim, t_end, dt, x0 = NULL) {
  if (!requireNamespace("deSolve", quietly = TRUE)) {
    stop("simulate_model() requires the 'deSolve' package. Install it via install.packages('deSolve').")
  }
  if (!requireNamespace("tibble", quietly = TRUE)) stop("simulate_model() requires 'tibble'.")
  if (!requireNamespace("dplyr", quietly = TRUE))  stop("simulate_model() requires 'dplyr'.")
  
  # row arrives as a 1-row tibble/data.frame from pmap; coerce safely
  if (is.data.frame(row)) row <- as.list(row[1, , drop = FALSE])
  if (!is.list(row)) stop("simulate_model(): 'row' must be list-like (a 1-row tibble/data.frame is fine).")
  
  logistic <- function(x) 1 / (1 + exp(-x))
  
  # Pull global defaults if you have them, else fallback
  p <- if (exists("par_base", inherits = TRUE)) get("par_base", inherits = TRUE) else list(
    kappa_min = 0.02,
    kappa_max = 0.25,
    kappa0    = 0.10,
    kappa1    = 30.0,
    phi3      = 8.0,
    phi4      = 1.0,
    phi0      = -0.02,
    alpha     = 0.02,
    phi1_min  = 0.10,
    phi1_max  = 5.00
  )
  
  # Candidate overrides
  if (!is.null(row$kappa_max) && is.finite(row$kappa_max)) p$kappa_max <- as.numeric(row$kappa_max)
  if (!is.null(row$kappa_min) && is.finite(row$kappa_min)) p$kappa_min <- as.numeric(row$kappa_min)
  
  # Core candidate parameters (must exist)
  sigma <- as.numeric(row$sigma)
  g_n   <- as.numeric(row$g_n)
  i     <- as.numeric(row$i)
  delta <- as.numeric(row$delta)
  
  psi  <- as.numeric(row$psi)
  phi2 <- as.numeric(row$phi2)
  
  stopifnot(
    is.finite(sigma), is.finite(g_n), is.finite(i), is.finite(delta),
    is.finite(psi), is.finite(phi2),
    g_n > 0, sigma > 0
  )
  
  kappa_fun <- function(r) {
    p$kappa_min + (p$kappa_max - p$kappa_min) * logistic(p$kappa1 * (r - p$kappa0))
  }
  Z_fun      <- function(d, f) logistic(p$phi3 * ((d - 1) + p$phi4 * (f - 1)))
  lambda_fun <- function(r, rF, psi) logistic(psi * (r - rF))
  
  # ---- Choose phi1 (constant) ----
  e_target <- if (exists("cfg", inherits = TRUE) &&
                  !is.null(get("cfg", inherits = TRUE)$targets$e_target)) {
    get("cfg", inherits = TRUE)$targets$e_target
  } else {
    0.94
  }
  
  phi1 <- NA_real_
  if (!is.null(row$phi1_endog) && is.finite(row$phi1_endog)) {
    phi1 <- as.numeric(row$phi1_endog)
  } else if (!is.null(row$Z_star) && is.finite(row$Z_star)) {
    phi1 <- (p$alpha - p$phi0 + phi2 * as.numeric(row$Z_star)) / e_target
  } else {
    phi1 <- 1.0
  }
  if (!is.finite(phi1)) phi1 <- 1.0
  
  # ---- Initial conditions ----
  if (is.null(x0)) {
    e0 <- if (!is.null(row$e_star) && is.finite(row$e_star)) as.numeric(row$e_star) else e_target
    w0 <- if (!is.null(row$omega_star) && is.finite(row$omega_star)) as.numeric(row$omega_star) else 0.65
    d0 <- if (!is.null(row$d_star) && is.finite(row$d_star)) as.numeric(row$d_star) else 0.5
    x0 <- c(e = e0, omega = w0, d = d0)
  } else {
    x0 <- as.numeric(x0)
    if (is.null(names(x0))) stop("x0 must be a named vector with e, omega, d.")
    if (!all(c("e", "omega", "d") %in% names(x0))) stop("x0 must be a named vector with e, omega, d.")
  }
  
  # ---- ODE system (reduced 3D) ----
  rhs <- function(t, state, parms) {
    e     <- state[["e"]]
    omega <- state[["omega"]]
    d     <- state[["d"]]
    
    r <- (1 - omega - i * d) / sigma
    
    kappa <- kappa_fun(r)
    g     <- kappa / sigma - delta
    
    lam <- lambda_fun(r, rF = rF_sim, psi = psi)
    lam <- min(max(lam, 1e-10), 1 - 1e-10)
    
    iotaF <- r * (lam / (1 - lam))
    f     <- iotaF / g_n
    
    Z <- Z_fun(d, f)
    
    de     <- (g - g_n) * e
    domega <- omega * (p$phi0 + phi1 * e - p$alpha - phi2 * Z)
    dd     <- kappa - (1 - omega) + i * d - d * g
    
    list(c(de, domega, dd))
  }
  
  times <- seq(0, t_end, by = dt)
  
  sol <- deSolve::ode(
    y = x0,
    times = times,
    func = rhs,
    parms = NULL,
    method = "lsoda"
  )
  
  # deSolve already returns `time` as a column; just select by string to avoid tidyselect junk
  sol <- as.data.frame(sol)
  tibble::as_tibble(sol) %>%
    dplyr::select("time", "e", "omega", "d")
}
