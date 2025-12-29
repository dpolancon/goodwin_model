############################################################
## reduced_goodwin_minsky_3D_integrated.R
## 3D reduced Goodwin–Minsky with quasi-steady financialization
##
## Features integrated (from your workflow upgrades):
##  - Closed-form interior steady state + admissibility gates
##  - Analytic Jacobian at SS + eigenvalues
##  - Routh–Hurwitz coefficients + Hopf functional H
##  - Hopf scan over rF_bar + bracketed roots (uniroot refine)
##  - Regime comparisons around Hopf:
##      * below/above root1, between roots, above root2, baseline
##      * COMMON initial condition across regimes (optional, default ON)
##  - Behavior classifier (fixed point vs damped vs limit cycle vs runaway/clamp)
##  - Reproducible outputs:
##      * CSV tables, plots, sessionInfo.txt, README.md
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

set.seed(123)

## -----------------------------
## 0) Helpers
## -----------------------------
inv_logit <- function(x) 1 / (1 + exp(-x))

require_params <- function(par, keys) {
  miss <- setdiff(keys, names(par))
  if (length(miss) > 0) stop("Missing parameters: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

# Accept legacy names but standardize internally (no silent renaming: legacy stored)
normalize_par <- function(par) {
  par <- as.list(par)
  
  legacy <- list()
  
  # Investment logistic: old names -> new names
  if (!is.null(par$lambda) && is.null(par$kappa_lambda)) {
    par$kappa_lambda <- par$lambda
    legacy$lambda <- par$lambda
  }
  if (!is.null(par$r0) && is.null(par$kappa_r0)) {
    par$kappa_r0 <- par$r0
    legacy$r0 <- par$r0
  }
  
  # Buffett discipline mapping if present (old code)
  if (!is.null(par$eta) && is.null(par$phi2)) {
    par$phi2 <- par$eta
    legacy$eta <- par$eta
  }
  if (!is.null(par$psi) && is.null(par$phi3)) {
    par$phi3 <- par$psi
    legacy$psi <- par$psi
  }
  
  # Legacy dbar parked (not in locked kernel)
  if (!is.null(par$dbar)) legacy$dbar <- par$dbar
  
  par$legacy_params <- legacy
  par
}

make_state0_from_ss <- function(ss, perturb = c(e=0.01, omega=-0.01, d=0.10)) {
  c(
    e     = ss$e     + perturb[["e"]],
    omega = ss$omega + perturb[["omega"]],
    d     = ss$d     + perturb[["d"]]
  )
}

## -----------------------------
## 1) Model blocks (LOCKED)
## -----------------------------
pi_fun <- function(omega) 1 - omega

r_fun <- function(omega, d, par) {
  (pi_fun(omega) - par$i * d) / par$sigma
}

# Investment logistic
kappa_fun <- function(r, par) {
  par$kappa_min + (par$kappa_max - par$kappa_min) /
    (1 + exp(-par$kappa_lambda * (r - par$kappa_r0)))
}

kappa_r_fun <- function(r, par) {
  s <- inv_logit(par$kappa_lambda * (r - par$kappa_r0))
  (par$kappa_max - par$kappa_min) * par$kappa_lambda * s * (1 - s)
}

# Growth
g_fun <- function(r, par) kappa_fun(r, par) / par$sigma - par$delta
g_r_fun <- function(r, par) kappa_r_fun(r, par) / par$sigma

# Payout share theta(r): logistic in (i - r)
theta_fun <- function(r, par) {
  x <- par$theta_k * ((par$i - r) - par$theta_0)
  par$theta_min + (par$theta_max - par$theta_min) * inv_logit(x)
}
theta_r_fun <- function(r, par) {
  x <- par$theta_k * ((par$i - r) - par$theta_0)
  s <- inv_logit(x)
  (par$theta_max - par$theta_min) * s * (1 - s) * (-par$theta_k)
}

# Portfolio tilt lambda(r; rF_tilde): logistic in (rF_tilde - r)
lambda_tilt_fun <- function(r, par) {
  x <- par$lam_k * ((par$rF_tilde - r) - par$lam_0)
  par$lam_min + (par$lam_max - par$lam_min) * inv_logit(x)
}
lambda_r_fun <- function(r, par) {
  x <- par$lam_k * ((par$rF_tilde - r) - par$lam_0)
  s <- inv_logit(x)
  (par$lam_max - par$lam_min) * s * (1 - s) * (-par$lam_k)
}

# Quasi-steady financialization fbar(r; rF_bar, rF_tilde)
fbar_fun <- function(r, par) {
  lam <- lambda_tilt_fun(r, par)
  th  <- theta_fun(r, par)
  num <- lam * (1 - th) * r
  den <- par$rF_bar - par$tau_F - g_fun(r, par)
  - num / den
}

fbar_r_fun <- function(r, par) {
  lam   <- lambda_tilt_fun(r, par)
  lam_r <- lambda_r_fun(r, par)
  th    <- theta_fun(r, par)
  th_r  <- theta_r_fun(r, par)
  
  num <- lam * (1 - th) * r
  den <- par$rF_bar - par$tau_F - g_fun(r, par)
  
  num_r <- lam_r * (1 - th) * r + lam * (-th_r) * r + lam * (1 - th) * 1
  den_r <- - g_r_fun(r, par)
  
  # f = - num/den
  - (num_r * den - num * den_r) / (den^2)
}

# Discipline kernel Z(d,f)
Z_fun <- function(d, f, par) {
  x <- par$phi3 * ((d - 1) + par$phi4 * (f - 1))
  inv_logit(x)
}
Z_d_fun <- function(d, f, par) {
  Z <- Z_fun(d, f, par)
  par$phi3 * Z * (1 - Z)
}
Z_f_fun <- function(d, f, par) {
  Z <- Z_fun(d, f, par)
  par$phi3 * par$phi4 * Z * (1 - Z)
}

# Composed discipline Ztilde(omega, d)
Ztilde_partials <- function(omega, d, par) {
  r_omega <- -1 / par$sigma
  r_d     <- -par$i / par$sigma
  
  r <- r_fun(omega, d, par)
  f <- fbar_fun(r, par)
  fb_r <- fbar_r_fun(r, par)
  
  Zd <- Z_d_fun(d, f, par)
  Zf <- Z_f_fun(d, f, par)
  
  Z_omega <- Zf * fb_r * r_omega
  Z_dtot  <- Zd + Zf * fb_r * r_d
  
  list(r=r, f=f, fb_r=fb_r, Z=Z_fun(d, f, par),
       Z_omega=Z_omega, Z_d=Z_dtot,
       r_omega=r_omega, r_d=r_d)
}

## -----------------------------
## 2) Interior steady state (LOCKED)
## -----------------------------
steady_state_reduced3D <- function(par, eps=1e-12, denom_eps=1e-6) {
  par <- normalize_par(par)
  
  require_params(par, c(
    "alpha","beta","delta","sigma","i",
    "kappa_min","kappa_max","kappa_lambda","kappa_r0",
    "phi0","phi1","phi2","phi3","phi4",
    "rF_bar","rF_tilde","tau_F",
    "theta_min","theta_max","theta_k","theta_0",
    "lam_min","lam_max","lam_k","lam_0"
  ))
  
  g_n <- par$alpha + par$beta
  if (g_n <= 0) return(list(ok=FALSE, reason="g_n <= 0"))
  
  k_star <- par$sigma * (g_n + par$delta)
  
  if (k_star <= par$kappa_min + eps || k_star >= par$kappa_max - eps) {
    return(list(ok=FALSE, reason="kappa_star outside (kappa_min, kappa_max)"))
  }
  
  s <- (k_star - par$kappa_min) / (par$kappa_max - par$kappa_min)
  if (s <= eps || s >= 1 - eps) return(list(ok=FALSE, reason="s not in (0,1)"))
  
  r_star <- par$kappa_r0 + (1 / par$kappa_lambda) * log(s / (1 - s))
  
  d_star <- (k_star - par$sigma * r_star) / g_n
  if (d_star <= eps) return(list(ok=FALSE, reason="d* <= 0"))
  
  omega_star <- 1 - par$i * d_star - par$sigma * r_star
  if (omega_star <= eps || omega_star >= 1 - eps) return(list(ok=FALSE, reason="omega* not in (0,1)"))
  
  g_star <- g_fun(r_star, par)
  denom  <- par$rF_bar - par$tau_F - g_star
  if (!is.finite(denom) || abs(denom) < denom_eps) {
    return(list(ok=FALSE, reason="denominator too close to 0 in fbar"))
  }
  
  f_star <- fbar_fun(r_star, par)
  Z_star <- Z_fun(d_star, f_star, par)
  
  if (par$phi1 <= 0) return(list(ok=FALSE, reason="phi1 must be > 0"))
  e_star <- (par$alpha - par$phi0 + par$phi2 * Z_star) / par$phi1
  if (e_star <= eps || e_star >= 1 - eps) return(list(ok=FALSE, reason="e* not in (0,1)"))
  
  list(ok=TRUE,
       e=e_star, omega=omega_star, d=d_star,
       r=r_star, kappa=k_star, g=g_star, g_n=g_n,
       denom=denom,
       theta=theta_fun(r_star, par),
       lambda_tilt=lambda_tilt_fun(r_star, par),
       f=f_star, Z=Z_star,
       par=par)
}

## -----------------------------
## 3) Jacobian + RH/Hopf + eigenvalues
## -----------------------------
jacobian_reduced3D_at_ss <- function(ss, par) {
  par <- normalize_par(par)
  e <- ss$e; omega <- ss$omega; d <- ss$d
  
  # composed discipline derivatives
  rpart <- Ztilde_partials(omega, d, par)
  r <- rpart$r
  
  kap_r <- kappa_r_fun(r, par)
  g_r   <- g_r_fun(r, par)
  
  r_omega <- rpart$r_omega
  r_d     <- rpart$r_d
  Z_omega <- rpart$Z_omega
  Z_d     <- rpart$Z_d
  
  # Eq1: e_dot = (g - g_n)*e, at SS g=g_n => J11=0
  J12 <- e * g_r * r_omega
  J13 <- e * g_r * r_d
  
  # Eq2: omega_dot = omega*(phi0 + phi1 e - alpha - phi2 Ztilde)
  # at SS bracket=0 => J22 = omega * (-phi2 * Z_omega), J23=omega*(-phi2*Z_d)
  J21 <- omega * par$phi1
  J22 <- omega * (-par$phi2 * Z_omega)
  J23 <- omega * (-par$phi2 * Z_d)
  
  # Eq3: d_dot = kappa(r) - (1-omega) + i d - d g(r)
  J32 <- 1 + r_omega * (kap_r - d * g_r)
  J33 <- kap_r * r_d + par$i - ss$g_n - d * g_r * r_d
  
  J <- matrix(c(
    0,   J12, J13,
    J21, J22, J23,
    0,   J32, J33
  ), nrow=3, byrow=TRUE)
  
  list(J=J,
       parts=list(J12=J12,J13=J13,J21=J21,J22=J22,J23=J23,J32=J32,J33=J33),
       aux=list(r=r, kap_r=kap_r, g_r=g_r,
                Z=rpart$Z, f=rpart$f, denom=ss$denom,
                Z_omega=Z_omega, Z_d=Z_d))
}

rh_hopf_3x3 <- function(J) {
  a1 <- -sum(diag(J))
  a2 <- det(J[-1,-1]) + det(J[-2,-2]) + det(J[-3,-3])
  a3 <- -det(J)
  H  <- a1*a2 - a3
  stable_RH <- (a1 > 0) && (a2 > 0) && (a3 > 0) && (H > 0)
  list(a1=a1, a2=a2, a3=a3, H=H, stable_RH=stable_RH)
}

analyze_at_par <- function(par, denom_eps=1e-6) {
  ss <- steady_state_reduced3D(par, denom_eps=denom_eps)
  if (!isTRUE(ss$ok)) return(ss)
  
  Jp <- jacobian_reduced3D_at_ss(ss, ss$par)
  rh <- rh_hopf_3x3(Jp$J)
  eig <- eigen(Jp$J)$values
  
  list(ok=TRUE, ss=ss, J=Jp$J, J_parts=Jp$parts, aux=Jp$aux, rh=rh, eigenvalues=eig)
}

## -----------------------------
## 4) Hopf scan over rF_bar (bracket + uniroot refine)
## -----------------------------
hopf_scan <- function(par,
                      rF_bar_grid = seq(0.005, 0.20, by=0.001),
                      denom_eps = 1e-6) {
  
  par <- normalize_par(par)
  
  Hvals <- rep(NA_real_, length(rF_bar_grid))
  a1v   <- rep(NA_real_, length(rF_bar_grid))
  a2v   <- rep(NA_real_, length(rF_bar_grid))
  a3v   <- rep(NA_real_, length(rF_bar_grid))
  okRH  <- rep(FALSE, length(rF_bar_grid))
  
  for (j in seq_along(rF_bar_grid)) {
    par$rF_bar <- rF_bar_grid[j]
    out <- analyze_at_par(par, denom_eps=denom_eps)
    if (!isTRUE(out$ok)) next
    
    rh <- out$rh
    if (rh$a1 > 0 && rh$a2 > 0 && rh$a3 > 0) {
      Hvals[j] <- rh$H
      a1v[j] <- rh$a1; a2v[j] <- rh$a2; a3v[j] <- rh$a3
      okRH[j] <- TRUE
    }
  }
  
  idx <- which(okRH & is.finite(Hvals))
  roots <- tibble(lower=numeric(0), upper=numeric(0), root=numeric(0))
  
  if (length(idx) >= 2) {
    for (k in 2:length(idx)) {
      i1 <- idx[k-1]; i2 <- idx[k]
      H1 <- Hvals[i1]; H2 <- Hvals[i2]
      if (is.finite(H1) && is.finite(H2) && H1*H2 < 0) {
        a <- rF_bar_grid[i1]; b <- rF_bar_grid[i2]
        
        fH <- function(rFbar) {
          par2 <- par
          par2$rF_bar <- rFbar
          out2 <- analyze_at_par(par2, denom_eps=denom_eps)
          if (!isTRUE(out2$ok)) return(NA_real_)
          rh2 <- out2$rh
          if (!(rh2$a1 > 0 && rh2$a2 > 0 && rh2$a3 > 0)) return(NA_real_)
          rh2$H
        }
        
        Ha <- fH(a); Hb <- fH(b)
        if (is.finite(Ha) && is.finite(Hb) && Ha*Hb < 0) {
          rt <- uniroot(fH, lower=a, upper=b)$root
          roots <- bind_rows(roots, tibble(lower=a, upper=b, root=rt))
        }
      }
    }
  }
  
  list(grid=rF_bar_grid, H=Hvals, a1=a1v, a2=a2v, a3=a3v, okRH=okRH, roots=roots)
}

## -----------------------------
## 5) RHS + simulation (lsoda + guardrails)
## -----------------------------
rhs_reduced3D <- function(t, state, par) {
  par <- normalize_par(par)
  
  e <- state[["e"]]
  w <- state[["omega"]]
  d <- state[["d"]]
  
  d_floor <- if (!is.null(par$d_floor)) par$d_floor else 1e-9
  d_safe  <- max(d, d_floor)
  e_safe  <- min(max(e, 1e-9), 1 - 1e-9)
  w_safe  <- min(max(w, 1e-9), 1 - 1e-9)
  
  r <- r_fun(w_safe, d_safe, par)
  kap <- kappa_fun(r, par)
  g   <- kap / par$sigma - par$delta
  g_n <- par$alpha + par$beta
  
  th  <- theta_fun(r, par)
  lam <- lambda_tilt_fun(r, par)
  den <- par$rF_bar - par$tau_F - g
  f   <- - (lam * (1 - th) * r) / den
  Zt  <- Z_fun(d_safe, f, par)
  
  de <- (g - g_n) * e_safe
  dw <- w_safe * (par$phi0 + par$phi1 * e_safe - par$alpha - par$phi2 * Zt)
  dd <- kap - (1 - w_safe) + par$i * d_safe - d_safe * g
  
  if (d <= d_floor && dd < 0) dd <- 0
  
  y_pc <- exp(par$alpha * t) * e_safe
  
  list(
    c(de, dw, dd),
    c(r=r, kappa=kap, g=g, g_n=g_n, theta=th, lambda_tilt=lam,
      denom=den, f=f, Ztilde=Zt, y_pc=y_pc)
  )
}

simulate_model <- function(par, state0=NULL,
                           t_end=500, dt=0.05,
                           perturb=c(e=0.01, omega=-0.01, d=0.10)) {
  
  par <- normalize_par(par)
  
  if (is.null(state0)) {
    ss <- steady_state_reduced3D(par)
    if (!isTRUE(ss$ok)) stop("No admissible steady state: ", ss$reason)
    state0 <- make_state0_from_ss(ss, perturb=perturb)
  }
  
  times <- seq(0, t_end, by=dt)
  
  sol <- ode(
    y = state0,
    times = times,
    func = rhs_reduced3D,
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
## 6) Plotting
## -----------------------------
plot_states <- function(df, tag, out_dir="outputs") {
  p <- df |>
    select(time, e, omega, d) |>
    pivot_longer(-time, names_to="var", values_to="value") |>
    ggplot(aes(x=time, y=value)) +
    geom_line() +
    facet_wrap(~var, scales="free_y", ncol=1) +
    labs(title=paste0("States over time (", tag, ")"), x="t", y=NULL)
  
  ggsave(file.path(out_dir, paste0("states_", tag, ".png")),
         p, width=7, height=7, dpi=160)
  p
}

plot_finance <- function(df, tag, out_dir="outputs") {
  p <- df |>
    select(time, r, kappa, g, theta, lambda_tilt, denom, f, Ztilde) |>
    pivot_longer(-time, names_to="var", values_to="value") |>
    ggplot(aes(x=time, y=value)) +
    geom_line() +
    facet_wrap(~var, scales="free_y", ncol=1) +
    labs(title=paste0("Finance/discipline block (", tag, ")"), x="t", y=NULL)
  
  ggsave(file.path(out_dir, paste0("finance_", tag, ".png")),
         p, width=7, height=9, dpi=160)
  p
}

plot_phases <- function(df, tag, out_dir="outputs") {
  p1 <- ggplot(df, aes(x=e, y=omega)) + geom_path() +
    labs(title=paste0("Phase: e–omega (", tag, ")"), x="e", y="omega")
  p2 <- ggplot(df, aes(x=omega, y=d)) + geom_path() +
    labs(title=paste0("Phase: omega–d (", tag, ")"), x="omega", y="d")
  p3 <- ggplot(df, aes(x=e, y=d)) + geom_path() +
    labs(title=paste0("Phase: e–d (", tag, ")"), x="e", y="d")
  
  ggsave(file.path(out_dir, paste0("phase_e_omega_", tag, ".png")), p1, width=6, height=4, dpi=160)
  ggsave(file.path(out_dir, paste0("phase_omega_d_", tag, ".png")), p2, width=6, height=4, dpi=160)
  ggsave(file.path(out_dir, paste0("phase_e_d_", tag, ".png")), p3, width=6, height=4, dpi=160)
  
  list(e_omega=p1, omega_d=p2, e_d=p3)
}

plot_3d_traj <- function(df, tag, out_dir="outputs", thin=5) {
  df2 <- df |> slice(seq(1, n(), by=thin))
  png(file.path(out_dir, paste0("traj3d_", tag, ".png")), width=900, height=700)
  scatterplot3d(df2$e, df2$omega, df2$d,
                type="l",
                xlab="e", ylab="omega", zlab="d",
                main=paste0("3D trajectory (", tag, ")"))
  dev.off()
  invisible(TRUE)
}

## -----------------------------
## 7) Simulation behavior classifier
## -----------------------------
classify_dynamics <- function(df, tail_frac=0.25) {
  n <- nrow(df)
  i0 <- max(1, floor((1 - tail_frac) * n))
  tail <- df[i0:n, , drop=FALSE]
  
  eps <- 1e-6
  hit_e     <- any(tail$e <= 0 + eps | tail$e >= 1 - eps)
  hit_omega <- any(tail$omega <= 0 + eps | tail$omega >= 1 - eps)
  hit_d0    <- any(tail$d <= 0 + eps)
  
  sd_e <- sd(tail$e, na.rm=TRUE)
  sd_w <- sd(tail$omega, na.rm=TRUE)
  sd_d <- sd(tail$d, na.rm=TRUE)
  sd_max <- max(sd_e, sd_w, sd_d)
  
  half <- floor(nrow(tail)/2)
  tail1 <- tail[1:half, , drop=FALSE]
  tail2 <- tail[(half+1):nrow(tail), , drop=FALSE]
  
  amp1 <- max(abs(tail1$omega - mean(tail1$omega, na.rm=TRUE)), na.rm=TRUE)
  amp2 <- max(abs(tail2$omega - mean(tail2$omega, na.rm=TRUE)), na.rm=TRUE)
  amp_ratio <- if (is.finite(amp1) && amp1 > 0) amp2 / amp1 else NA_real_
  
  if (hit_e || hit_omega || hit_d0) {
    cls <- "boundary_crash_or_clamp"
  } else if (is.finite(sd_max) && sd_max < 1e-5) {
    cls <- "fixed_point_convergence"
  } else if (is.finite(amp_ratio) && amp_ratio > 1.1) {
    cls <- "diverging_oscillation"
  } else if (is.finite(amp_ratio) && amp_ratio < 0.9) {
    cls <- "damped_oscillation"
  } else {
    cls <- "limit_cycle_or_neutral"
  }
  
  tibble(
    tail_sd_e = sd_e, tail_sd_omega = sd_w, tail_sd_d = sd_d,
    amp_ratio = amp_ratio,
    hit_bounds = hit_e || hit_omega || hit_d0,
    class = cls
  )
}

## -----------------------------
## 8) Bundle runner: Hopf scan + regime comparisons
## -----------------------------
run_bundle <- function(par,
                       out_dir="outputs/reduced3D",
                       t_end=500, dt=0.05,
                       hopf_grid=seq(0.005, 0.20, by=0.001),
                       hopf_eps=0.002,
                       denom_eps=1e-6,
                       common_state0 = TRUE,
                       ref_tag_preference = c("between_roots", "above_root1", "baseline")) {
  
  par <- normalize_par(par)
  dir_create(out_dir)
  
  # session info
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  
  # baseline diagnostics
  base <- analyze_at_par(par, denom_eps=denom_eps)
  if (!isTRUE(base$ok)) stop("Baseline inadmissible: ", base$reason)
  
  ss_tbl <- tibble(
    e=base$ss$e, omega=base$ss$omega, d=base$ss$d,
    r=base$ss$r, kappa=base$ss$kappa, g=base$ss$g, g_n=base$ss$g_n,
    denom=base$ss$denom,
    theta=base$ss$theta, lambda_tilt=base$ss$lambda_tilt,
    f=base$ss$f, Z=base$ss$Z,
    rF_bar=par$rF_bar, rF_tilde=par$rF_tilde, tau_F=par$tau_F
  )
  write.csv(ss_tbl, file.path(out_dir, "steady_state.csv"), row.names=FALSE)
  
  write.csv(as.data.frame(base$J) |> mutate(row=paste0("row",1:3)) |> relocate(row),
            file.path(out_dir, "jacobian.csv"), row.names=FALSE)
  
  write.csv(tibble(a1=base$rh$a1, a2=base$rh$a2, a3=base$rh$a3, H=base$rh$H, stable_RH=base$rh$stable_RH),
            file.path(out_dir, "rh_hopf.csv"), row.names=FALSE)
  
  write.csv(tibble(eigen_real=Re(base$eigenvalues), eigen_imag=Im(base$eigenvalues)),
            file.path(out_dir, "eigenvalues.csv"), row.names=FALSE)
  
  # Hopf scan
  scan <- hopf_scan(par, rF_bar_grid=hopf_grid, denom_eps=denom_eps)
  scan_tbl <- tibble(rF_bar=scan$grid, H=scan$H, a1=scan$a1, a2=scan$a2, a3=scan$a3, okRH=scan$okRH)
  write.csv(scan_tbl, file.path(out_dir, "hopf_scan.csv"), row.names=FALSE)
  write.csv(scan$roots, file.path(out_dir, "hopf_roots.csv"), row.names=FALSE)
  
  pH <- ggplot(scan_tbl |> filter(okRH, is.finite(H)), aes(x=rF_bar, y=H)) +
    geom_line() + geom_hline(yintercept=0) +
    labs(title="Hopf functional H over rF_bar (RH-positive region)", x="rF_bar", y="H")
  ggsave(file.path(out_dir, "hopf_H_curve.png"), pH, width=7, height=4, dpi=160)
  
  # Build regimes around Hopf roots
  pars_list <- list(baseline = par)
  
  if (nrow(scan$roots) >= 1) {
    r1 <- scan$roots$root[1]
    pars_list$below_root1 <- modifyList(par, list(rF_bar = r1 - hopf_eps))
    pars_list$above_root1 <- modifyList(par, list(rF_bar = r1 + hopf_eps))
  }
  if (nrow(scan$roots) >= 2) {
    r1 <- scan$roots$root[1]
    r2 <- scan$roots$root[2]
    pars_list$between_roots <- modifyList(par, list(rF_bar = 0.5*(r1 + r2)))
    pars_list$above_root2   <- modifyList(par, list(rF_bar = r2 + hopf_eps))
  }
  
  # Choose reference initial condition (for clean comparisons)
  state0_ref <- NULL
  ref_tag <- NULL
  if (common_state0) {
    for (cand in ref_tag_preference) {
      if (!is.null(pars_list[[cand]])) {
        out_ref <- analyze_at_par(pars_list[[cand]], denom_eps=denom_eps)
        if (isTRUE(out_ref$ok)) {
          state0_ref <- make_state0_from_ss(out_ref$ss)
          ref_tag <- cand
          break
        }
      }
    }
    if (is.null(state0_ref)) {
      # fallback: baseline SS
      state0_ref <- make_state0_from_ss(base$ss)
      ref_tag <- "baseline"
    }
  }
  
  # Regime summary table
  regime_summary <- tibble(
    tag = character(), rF_bar = numeric(), ok = logical(),
    H = numeric(), a1 = numeric(), a2 = numeric(), a3 = numeric(),
    stable_RH = logical(),
    tail_sd_e = numeric(), tail_sd_omega = numeric(), tail_sd_d = numeric(),
    amp_ratio = numeric(), hit_bounds = logical(), class = character()
  )
  
  # Run each regime
  for (tg in names(pars_list)) {
    par_t <- normalize_par(pars_list[[tg]])
    out_t <- analyze_at_par(par_t, denom_eps=denom_eps)
    
    if (!isTRUE(out_t$ok)) {
      regime_summary <- bind_rows(regime_summary, tibble(
        tag=tg, rF_bar=par_t$rF_bar, ok=FALSE,
        H=NA_real_, a1=NA_real_, a2=NA_real_, a3=NA_real_, stable_RH=FALSE,
        tail_sd_e=NA_real_, tail_sd_omega=NA_real_, tail_sd_d=NA_real_,
        amp_ratio=NA_real_, hit_bounds=NA, class=NA_character_
      ))
      next
    }
    
    # Simulation (common initial condition if requested)
    df <- simulate_model(
      par_t,
      state0 = if (common_state0) state0_ref else NULL,
      t_end = t_end, dt = dt
    )
    
    # Save simulation data
    write.csv(df, file.path(out_dir, paste0("sim_", tg, ".csv")), row.names=FALSE)
    
    # Plots
    plot_states(df, tg, out_dir)
    plot_finance(df, tg, out_dir)
    plot_phases(df, tg, out_dir)
    plot_3d_traj(df, tg, out_dir, thin=5)
    
    # Classifier
    cls <- classify_dynamics(df, tail_frac=0.25)
    
    rh <- out_t$rh
    regime_summary <- bind_rows(regime_summary, tibble(
      tag=tg, rF_bar=par_t$rF_bar, ok=TRUE,
      H=rh$H, a1=rh$a1, a2=rh$a2, a3=rh$a3, stable_RH=rh$stable_RH,
      tail_sd_e=cls$tail_sd_e, tail_sd_omega=cls$tail_sd_omega, tail_sd_d=cls$tail_sd_d,
      amp_ratio=cls$amp_ratio, hit_bounds=cls$hit_bounds, class=cls$class
    ))
  }
  
  write.csv(regime_summary, file.path(out_dir, "regime_summary.csv"), row.names=FALSE)
  print(regime_summary)
  
  # README
  readme <- c(
    "# Reduced 3D Goodwin–Minsky (quasi-steady financialization) outputs",
    "",
    "## What this run does",
    "- Computes admissible interior steady state (closed-form) + saves steady_state.csv",
    "- Builds Jacobian at SS, RH/Hopf functional H, eigenvalues",
    "- Scans rF_bar for Hopf crossings (RH-positivity gated) and refines roots",
    "- Runs simulations for regimes around Hopf roots + baseline",
    "- Uses a common initial condition across regimes if common_state0=TRUE",
    "- Classifies tail behavior (fixed point / damped / limit cycle / runaway/clamp)",
    "",
    "## Key files",
    "- steady_state.csv, jacobian.csv, rh_hopf.csv, eigenvalues.csv",
    "- hopf_scan.csv, hopf_roots.csv, hopf_H_curve.png",
    "- regime_summary.csv",
    "- sim_<tag>.csv and plots: states_*, finance_*, phase_*, traj3d_*",
    "- sessionInfo.txt",
    "",
    "## Notes",
    paste0("- common_state0 = ", common_state0, " (reference tag: ", ifelse(is.null(ref_tag), "NA", ref_tag), ")")
  )
  writeLines(readme, file.path(out_dir, "README.md"))
  
  invisible(list(base=base, scan=scan, regimes=regime_summary, state0_ref=state0_ref, ref_tag=ref_tag))
}

############################################################
## 9) Baseline parameters (old example + new defaults)
############################################################
par <- list(
  # macro
  alpha = 0.02, beta = 0.01,
  delta = 0.02, sigma = 3,
  i = 0.04,
  
  # investment logistic
  kappa_min = 0.00, kappa_max = 0.30,
  kappa_lambda = 20, kappa_r0 = 0.04,
  
  # wage dynamics
  phi0 = -0.06, phi1 = 0.10,
  
  # discipline
  phi2 = 0.02,   # intensity
  phi3 = 10,     # steepness
  phi4 = 1.0,    # loading on (f-1)
  
  # financialization sweep lever
  rF_bar = 0.08,
  rF_tilde = 0.06,
  tau_F = 0.01,
  
  # payout theta(r)
  theta_min = 0.10, theta_max = 0.90,
  theta_k = 20, theta_0 = 0.00,
  
  # portfolio tilt lambda(r)
  lam_min = 0.10, lam_max = 0.90,
  lam_k = 20, lam_0 = 0.00,
  
  # numerics
  d_floor = 1e-9
)

############################################################
## 10) Execute
############################################################
bundle <- run_bundle(
  par,
  out_dir = "outputs/reduced3D",
  t_end = 500,
  dt = 0.05,
  hopf_grid = seq(0.005, 0.20, by = 0.001),
  hopf_eps = 0.002,
  common_state0 = TRUE
)

cat("\nDONE. See outputs/reduced3D/\n")


## ============================================================
## 11) Post-run: load sims, build SFC-style auxiliaries, compare regimes
## ============================================================

load_simulations <- function(out_dir, tags) {
  dfs <- lapply(tags, function(tg) {
    f <- file.path(out_dir, paste0("sim_", tg, ".csv"))
    if (!file.exists(f)) return(NULL)
    read.csv(f) |>
      as_tibble() |>
      mutate(tag = tg)
  })
  bind_rows(dfs)
}

add_aux_sfc <- function(df, par) {
  par <- normalize_par(par)
  
  df |>
    mutate(
      # core derived
      pi = 1 - omega,
      r_check = (pi - par$i * d) / par$sigma,  # should match r column closely
      g_check = kappa / par$sigma - par$delta, # should match g
      
      # debt ratio identity check:
      # ddot = kappa - (1-omega) + i d - d g
      d_dot = kappa - (1 - omega) + par$i * d - d * g,
      # residual of the identity (should be ~0 if columns consistent):
      sfc_resid = (d_dot - (kappa - (1 - omega) + par$i * d - d * g)),
      
      # financialization mechanics
      num_f = lambda_tilt * (1 - theta) * r,
      den_f = denom,
      f_check = - num_f / den_f,              # should match f
      f_err = f - f_check,
      
      # discipline argument
      z_arg = (d - 1) + par$phi4 * (f - 1)
    )
}

## ---- A) Overlay plots: states compared across regimes
plot_compare_states <- function(df_all, out_dir) {
  p <- df_all |>
    select(time, tag, e, omega, d) |>
    pivot_longer(cols = c(e, omega, d), names_to = "var", values_to = "value") |>
    ggplot(aes(x = time, y = value, group = tag)) +
    geom_line(aes(linetype = tag)) +
    facet_wrap(~var, scales = "free_y", ncol = 1) +
    labs(title = "Regime comparison: states over time", x = "t", y = NULL)
  
  ggsave(file.path(out_dir, "compare_states.png"), p, width = 8, height = 8, dpi = 160)
  p
}

## ---- B) Overlay plots: key auxiliaries compared across regimes
plot_compare_key_aux <- function(df_all, out_dir) {
  p <- df_all |>
    select(time, tag, r, g, kappa, denom, f, Ztilde, theta, lambda_tilt) |>
    pivot_longer(-c(time, tag), names_to = "var", values_to = "value") |>
    ggplot(aes(x = time, y = value, group = tag)) +
    geom_line(aes(linetype = tag)) +
    facet_wrap(~var, scales = "free_y", ncol = 2) +
    labs(title = "Regime comparison: key auxiliaries", x = "t", y = NULL)
  
  ggsave(file.path(out_dir, "compare_aux_key.png"), p, width = 10, height = 8, dpi = 160)
  p
}

## ---- C) “Everything panel”: include SFC checks + f mechanics
plot_compare_full_aux <- function(df_all, out_dir) {
  vars <- c(
    "pi","r","r_check","g","g_check","kappa",
    "theta","lambda_tilt","denom",
    "num_f","f","f_check","f_err",
    "Ztilde","z_arg",
    "d_dot","sfc_resid"
  )
  keep <- intersect(vars, names(df_all))
  
  p <- df_all |>
    select(any_of(c("time","tag", keep))) |>
    pivot_longer(-c(time, tag), names_to = "var", values_to = "value") |>
    ggplot(aes(x = time, y = value, group = tag)) +
    geom_line(aes(linetype = tag)) +
    facet_wrap(~var, scales = "free_y", ncol = 3) +
    labs(title = "Regime comparison: full auxiliary + SFC checks", x = "t", y = NULL)
  
  ggsave(file.path(out_dir, "compare_aux_full.png"), p, width = 14, height = 10, dpi = 160)
  p
}

## ---- D) Regime report table: bounds, denom risk, SFC error, tail stats
regime_report <- function(df_all, denom_floor = 1e-3) {
  df_all |>
    group_by(tag) |>
    summarize(
      n = n(),
      t_end = max(time, na.rm = TRUE),
      
      # bound hits (using same eps logic)
      hit_e = any(e <= 1e-6 | e >= 1 - 1e-6, na.rm = TRUE),
      hit_omega = any(omega <= 1e-6 | omega >= 1 - 1e-6, na.rm = TRUE),
      hit_d = any(d <= 1e-6, na.rm = TRUE),
      
      # denom danger (this is usually what kills your “financialization” runs)
      denom_min_abs = min(abs(denom), na.rm = TRUE),
      denom_near0 = any(abs(denom) < denom_floor, na.rm = TRUE),
      
      # sanity checks
      r_max_abs_err = max(abs(r - r_check), na.rm = TRUE),
      g_max_abs_err = max(abs(g - g_check), na.rm = TRUE),
      f_max_abs_err = max(abs(f - f_check), na.rm = TRUE),
      
      # SFC identity residual (should be ~0, but keep it anyway)
      sfc_resid_max = max(abs(sfc_resid), na.rm = TRUE),
      
      # “volatility” snapshot
      sd_e = sd(e, na.rm = TRUE),
      sd_omega = sd(omega, na.rm = TRUE),
      sd_d = sd(d, na.rm = TRUE),
      
      .groups = "drop"
    )
}



## ============================================================
## 12) Call these after run_bundle()
## ============================================================
# Example usage (after bundle <- run_bundle(...)):
tags <- bundle$regimes$tag
df_all <- load_simulations("outputs/reduced3D", tags)
df_all <- add_aux_sfc(df_all, par)   # uses your baseline par for constants
write.csv(df_all, file.path("outputs/reduced3D", "all_regimes_long.csv"), row.names = FALSE)
#
plot_compare_states(df_all, "outputs/reduced3D")
plot_compare_key_aux(df_all, "outputs/reduced3D")
plot_compare_full_aux(df_all, "outputs/reduced3D")

rep <- regime_report(df_all, denom_floor = 1e-3)
write.csv(rep, file.path("outputs/reduced3D", "regime_report.csv"), row.names = FALSE)
print(rep)
