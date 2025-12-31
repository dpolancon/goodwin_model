# ============================================================
# fast_calibration_sigma_fixed_targetOmega_v2.R
# FIX: make omega_target feasible (avoid d*<0) with sigma=2.25, g_n=0.05
# by adjusting i upward ONLY if needed, then recalibrate kappa0 (closed-form).
# Then scan (rF, psi, phi2) and set phi1 endogenously to hit e_target.
#
# Outputs to:
# outputs/wealth_goodwin/grid_search/best_results_targetOmega/
#   - results.csv
#   - results_augmented.csv
#   - hopf_roots.csv
#   - best_candidates_top50.csv
#   - plots/*.png
#   - sessionInfo.txt
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

# ----------------------------
# Paths
# ----------------------------
out_dir  <- "outputs/wealth_goodwin/grid_search/best_results_targetOmega"
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Targets (advanced-economy-ish)
# ----------------------------
TARGETS <- list(
  e_target = 0.94,
  e_min = 0.88, e_max = 0.98,
  
  omega_target = 0.65,
  omega_min = 0.62, omega_max = 0.70,
  
  d_max = 6.0,
  f_max = 10.0,
  
  phi1_min = 0.30,
  phi1_max = 1.50
)

# ----------------------------
# Utility
# ----------------------------
clamp01 <- function(x) pmin(1, pmax(0, x))
rf_tilde_from_rf <- function(rF, psi) -clamp01(psi) * rF

score_band <- function(x, lo, hi) {
  ifelse(is.na(x), Inf,
         ifelse(x < lo, (lo - x)^2,
                ifelse(x > hi, (x - hi)^2, 0)))
}

# ----------------------------
# Logistic bounded + inverse
# ----------------------------
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
# Model base parameters
# ----------------------------
par_base <- list(
  # macro
  sigma = 2.25,  # FIXED as requested
  delta = 0.05,
  
  # Make g_n = alpha + beta = 0.05 as requested
  alpha = 0.03,
  beta  = 0.02,
  
  # Interest (will be increased ONLY if omega_target is infeasible)
  i     = 0.03,
  
  # wage-share block
  phi0 = 0.00,
  phi1 = 0.60,  # overwritten endogenously
  phi2 = 0.70,
  phi3 = 5.0,
  phi4 = 1.0,
  
  # investment share kappa(r): bounded logistic
  kappa_min = 0.01,
  kappa_max = 0.35,
  kappa0    = 0.05,   # WILL BE RECALIBRATED to hit omega_target
  kappa1    = 40.0,
  
  # portfolio tilt lambda(rF_eff - r - lam0)
  lam0 = 0.00,
  lam1 = 40.0,
  
  # finance
  rF = 0.08,
  rF_tilde = 0.00,
  psi = 0.20
)

# ----------------------------
# Functions
# ----------------------------
kappa_fun <- function(r, p) logistic_bounded(r, p$kappa_min, p$kappa_max, p$kappa0, p$kappa1)
kappa_r   <- function(r, p) logistic_bounded_d(r, p$kappa_min, p$kappa_max, p$kappa0, p$kappa1)
g_fun     <- function(r, p) kappa_fun(r, p) / p$sigma - p$delta

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
# Backbone steady-state (given kappa0)
# ----------------------------
steady_state_backbone <- function(p) {
  gn <- p$alpha + p$beta
  if (!is.finite(gn) || gn <= 0) return(list(ok = FALSE, reason = "gn_nonpositive"))
  
  kappa_star <- p$sigma * (gn + p$delta)
  if (!(kappa_star > p$kappa_min && kappa_star < p$kappa_max)) {
    return(list(ok = FALSE, reason = "kappa_star_out_of_bounds"))
  }
  
  r_star <- logistic_bounded_inv(kappa_star, p$kappa_min, p$kappa_max, p$kappa0, p$kappa1)
  if (!is.finite(r_star)) return(list(ok = FALSE, reason = "r_star_nonfinite"))
  
  d_star <- (kappa_star - p$sigma * r_star) / gn
  omega_star <- 1 - p$i * d_star - p$sigma * r_star
  
  admiss <- TRUE
  reasons <- character(0)
  if (!is.finite(d_star) || !is.finite(omega_star)) { admiss <- FALSE; reasons <- c(reasons, "nonfinite_backbone") }
  if (d_star < 0) { admiss <- FALSE; reasons <- c(reasons, "d_star_negative") }
  if (omega_star <= 0 || omega_star >= 1) { admiss <- FALSE; reasons <- c(reasons, "omega_star_outside_0_1") }
  
  list(
    ok = admiss,
    reason = if (admiss) "OK" else paste(reasons, collapse = ";"),
    gn = gn,
    kappa_star = kappa_star,
    r_star = r_star,
    d_star = d_star,
    omega_star = omega_star
  )
}

# ----------------------------
# Closed-form solve for r* from omega_target (backbone identities)
# omega = 1 - i d - sigma r
# d = (kappa* - sigma r)/gn = sigma(gn+delta-r)/gn
# => omega = 1 - (i*sigma/gn)(gn+delta) + sigma r (i/gn - 1)
# => r = [omega - 1 + (i*sigma/gn)(gn+delta)] / [sigma (i/gn - 1)]
# ----------------------------
r_from_omega_closed <- function(sigma, gn, delta, i, omega) {
  denom <- sigma * (i/gn - 1)
  if (abs(denom) < 1e-10) return(NA_real_)
  (omega - 1 + (i*sigma/gn)*(gn + delta)) / denom
}

# ----------------------------
# Ensure omega_target is feasible with d>=0 and r>r_min by adjusting i upward if needed.
# This is NOT a new feature, it's a feasibility correction for your target tuple.
# ----------------------------
enforce_feasible_i_for_targets <- function(p, omega_target, d_max, r_min = 0.01) {
  gn <- p$alpha + p$beta
  sig <- p$sigma
  del <- p$delta
  
  omega_min_at_d0 <- 1 - sig * (gn + del)
  
  # If omega_target < omega_min_at_d0 AND gn > i, impossible with d>=0.
  # Need i > gn. Also must keep r* >= r_min.
  i0 <- p$i
  
  if (omega_target >= omega_min_at_d0) {
    # already feasible with d>=0 under gn>i case; keep i
    return(p)
  }
  
  # We need i>gn. We'll find minimal i that gives r*(omega_target) >= r_min
  # AND d*(omega_target) <= d_max AND d>=0
  # Using closed-form relationships for r and d.
  f_ok <- function(i_try) {
    if (!is.finite(i_try) || i_try <= gn + 1e-6) return(FALSE)
    r_try <- r_from_omega_closed(sig, gn, del, i_try, omega_target)
    if (!is.finite(r_try)) return(FALSE)
    
    kappa_star <- sig * (gn + del)
    d_try <- (kappa_star - sig * r_try) / gn
    
    # constraints
    if (!(d_try >= 0 && d_try <= d_max)) return(FALSE)
    if (!(r_try >= r_min)) return(FALSE)
    
    TRUE
  }
  
  # If current i works (rare here), keep it.
  if (f_ok(i0)) return(p)
  
  # Search i in a reasonable band (you can widen if you want)
  # Start from max(i0, gn+0.001) upward.
  i_lo <- max(i0, gn + 0.001)
  i_hi <- 0.50  # hard cap; if we hit this, your targets are inconsistent with the backbone.
  
  # Coarse grid to find any feasible i
  grid_i <- seq(i_lo, i_hi, by = 0.001)
  ok_idx <- which(vapply(grid_i, f_ok, logical(1)))
  
  if (length(ok_idx) == 0) {
    stop(paste0(
      "Targets infeasible with current backbone even after i search. ",
      "Given sigma=", sig, ", gn=", gn, ", delta=", del,
      ", omega_target=", omega_target,
      ", need either higher sigma, different omega_target, or different (gn,delta) assumptions."
    ))
  }
  
  p$i <- grid_i[min(ok_idx)]
  p
}

# ----------------------------
# Recalibrate kappa0 (closed-form) to hit omega_target given feasible i
# ----------------------------
recalibrate_kappa0_for_omega_closed <- function(p, omega_target) {
  gn <- p$alpha + p$beta
  
  kappa_star <- p$sigma * (gn + p$delta)
  if (!(kappa_star > p$kappa_min && kappa_star < p$kappa_max)) {
    stop("kappa_star not in (kappa_min, kappa_max). Adjust sigma/gn/delta/bounds.")
  }
  
  r_star <- r_from_omega_closed(p$sigma, gn, p$delta, p$i, omega_target)
  if (!is.finite(r_star)) stop("Closed-form r_star failed (likely i ~ gn).")
  
  # Convert r_star into kappa0 via logistic inversion algebra:
  # r = kappa0 - (1/kappa1) * ln((b-a)/(kappa_star-a) - 1)
  # => kappa0 = r + (1/kappa1)*ln((b-a)/(kappa_star-a) - 1)
  term <- (p$kappa_max - p$kappa_min) / (kappa_star - p$kappa_min) - 1
  if (!(is.finite(term) && term > 0)) stop("Bad inversion term; check bounds.")
  
  p$kappa0 <- r_star + (1 / p$kappa1) * log(term)
  
  bb <- steady_state_backbone(p)
  if (!isTRUE(bb$ok)) {
    stop(paste0("kappa0 solved but backbone infeasible: ", bb$reason))
  }
  
  list(p = p, bb = bb)
}

# ----------------------------
# Full steady-state with finance block + endogenous phi1 (targets e*)
# Option 1: f* = iotaF*/g_n, iotaF* = lambda* sF*, sF* = r*/(1-lambda*)
# ----------------------------
steady_state_full <- function(p, bb, phi2, e_target) {
  gn <- bb$gn
  r_star <- bb$r_star
  d_star <- bb$d_star
  omega_star <- bb$omega_star
  
  lam_star <- lambda_fun(r_star, p)
  if (!is.finite(lam_star) || lam_star >= 0.999) {
    return(list(ok = FALSE, reason = "lambda_too_close_to_1"))
  }
  
  sF_star <- r_star / (1 - lam_star)
  iotaF_star <- lam_star * sF_star
  f_star <- iotaF_star / gn
  
  p$phi2 <- phi2
  Z_star <- Z_fun(d_star, f_star, p)
  
  numer <- (p$alpha - p$phi0 + phi2 * Z_star)
  if (!is.finite(numer) || numer <= 0) {
    return(list(ok = FALSE, reason = "phi1_numer_nonpositive"))
  }
  
  phi1 <- numer / e_target
  e_star <- numer / phi1
  
  admiss <- TRUE
  reasons <- character(0)
  if (!all(is.finite(c(f_star, Z_star, phi1, e_star)))) { admiss <- FALSE; reasons <- c(reasons, "nonfinite_state") }
  if (f_star < 0) { admiss <- FALSE; reasons <- c(reasons, "f_star_negative") }
  if (phi1 < TARGETS$phi1_min || phi1 > TARGETS$phi1_max) { admiss <- FALSE; reasons <- c(reasons, "phi1_out_of_bounds") }
  
  list(
    ok = admiss,
    reason = if (admiss) "OK" else paste(reasons, collapse = ";"),
    r_star = r_star,
    d_star = d_star,
    omega_star = omega_star,
    lam_star = lam_star,
    sF_star = sF_star,
    iotaF_star = iotaF_star,
    f_star = f_star,
    Z_star = Z_star,
    phi1 = phi1,
    e_star = e_star
  )
}

# ----------------------------
# Jacobian at steady state
# ----------------------------
jacobian_at_ss <- function(ss, bb, p, phi2) {
  gn <- bb$gn
  r  <- ss$r_star
  d  <- ss$d_star
  e  <- ss$e_star
  om <- ss$omega_star
  f  <- ss$f_star
  phi1 <- ss$phi1
  
  r_om <- -1 / p$sigma
  r_d  <- -p$i / p$sigma
  
  kap_rv <- kappa_r(r, p)
  g_r   <- kap_rv / p$sigma
  
  sF <- ss$sF_star
  lam <- ss$lam_star
  lam_rv <- lambda_r(r, p)
  
  denom <- g_fun(r, p)
  if (abs(denom) < 1e-10) denom <- 1e-10
  
  f_r <- sF * (lam_rv * denom - lam * g_r) / (denom^2)
  
  Zd <- Z_d_fun(d, f, p)
  Zf <- Z_f_fun(d, f, p)
  
  Z_om <- Zf * f_r * r_om
  Z_dv <- Zd + Zf * f_r * r_d
  
  J11 <- 0
  J12 <- e * g_r * r_om
  J13 <- e * g_r * r_d
  
  J21 <- om * phi1
  J22 <- -om * phi2 * Z_om
  J23 <- -om * phi2 * Z_dv
  
  J31 <- 0
  J32 <- 1 + kap_rv * r_om - d * g_r * r_om
  J33 <- kap_rv * r_d + p$i - gn - d * g_r * r_d
  
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
  
  tibble(a1 = a1, a2 = a2, a3 = a3, H = H, tr = tr, det = detJ)
}

# ============================================================
# 0) Make targets feasible + recalibrate kappa0
# ============================================================
par_base <- enforce_feasible_i_for_targets(
  par_base,
  omega_target = TARGETS$omega_target,
  d_max = TARGETS$d_max,
  r_min = 0.01
)

cal <- recalibrate_kappa0_for_omega_closed(par_base, omega_target = TARGETS$omega_target)
par_base <- cal$p
bb0 <- cal$bb

message("Backbone calibrated:")
message("  sigma = ", par_base$sigma)
message("  g_n   = ", bb0$gn)
message("  i     = ", par_base$i)
message("  kappa0= ", par_base$kappa0)
message("  r*    = ", bb0$r_star, "  d* = ", bb0$d_star, "  omega* = ", bb0$omega_star)

# ============================================================
# 1) Grids (keep small & informative)
# ============================================================
rF_grid   <- seq(0.03, 0.12, by = 0.003)
psi_grid  <- seq(0.00, 0.40, by = 0.02)
phi2_grid <- seq(0.50, 0.90, by = 0.05)

# ============================================================
# 2) Run search
# ============================================================
rows <- vector("list", length(rF_grid) * length(psi_grid) * length(phi2_grid))
k <- 0L

for (psi in psi_grid) {
  for (rF in rF_grid) {
    p <- par_base
    p$psi <- clamp01(psi)
    p$rF  <- rF
    p$rF_tilde <- rf_tilde_from_rf(rF, p$psi)
    rF_eff <- p$rF + p$rF_tilde
    
    for (phi2 in phi2_grid) {
      k <- k + 1L
      
      ss <- steady_state_full(p, bb0, phi2 = phi2, e_target = TARGETS$e_target)
      
      if (!isTRUE(ss$ok)) {
        rows[[k]] <- tibble(
          sigma = p$sigma, gn = bb0$gn, i = p$i, kappa0 = p$kappa0,
          psi = p$psi, rF = rF, rF_tilde = p$rF_tilde, rF_eff = rF_eff,
          phi2 = phi2, phi1 = NA_real_,
          feasible = FALSE, reason = ss$reason,
          r_star = bb0$r_star, d_star = bb0$d_star, omega_star = bb0$omega_star,
          e_star = NA_real_, f_star = NA_real_, Z_star = NA_real_,
          maxReEig = NA_real_, maxImEig = NA_real_,
          a1 = NA_real_, a2 = NA_real_, a3 = NA_real_, H = NA_real_
        )
        next
      }
      
      J <- jacobian_at_ss(ss, bb0, p, phi2 = phi2)
      ev <- eigen(J, only.values = TRUE)$values
      rh <- rh_coeffs(J)
      
      rows[[k]] <- tibble(
        sigma = p$sigma, gn = bb0$gn, i = p$i, kappa0 = p$kappa0,
        psi = p$psi, rF = rF, rF_tilde = p$rF_tilde, rF_eff = rF_eff,
        phi2 = phi2, phi1 = ss$phi1,
        feasible = TRUE, reason = "OK",
        r_star = ss$r_star, d_star = ss$d_star, omega_star = ss$omega_star,
        e_star = ss$e_star, f_star = ss$f_star, Z_star = ss$Z_star,
        maxReEig = max(Re(ev)),
        maxImEig = max(abs(Im(ev))),
        a1 = rh$a1, a2 = rh$a2, a3 = rh$a3, H = rh$H
      )
    }
  }
}

res <- bind_rows(rows)

# ============================================================
# 3) Augment classifications + score
# ============================================================
res_aug <- res %>%
  mutate(
    stable = feasible & !is.na(maxReEig) & (maxReEig < 0),
    oscillatory = feasible & !is.na(maxImEig) & (maxImEig > 1e-6),
    rh_ok = feasible & !is.na(a1) & !is.na(a2) & !is.na(a3) & (a1 > 0) & (a2 > 0) & (a3 > 0),
    period = ifelse(oscillatory, 2*pi/maxImEig, NA_real_),
    
    econ_ok = feasible &
      !is.na(e_star) & !is.na(omega_star) & !is.na(d_star) & !is.na(f_star) & !is.na(phi1) &
      (e_star >= TARGETS$e_min) & (e_star <= TARGETS$e_max) &
      (omega_star >= TARGETS$omega_min) & (omega_star <= TARGETS$omega_max) &
      (d_star >= 0) & (d_star <= TARGETS$d_max) &
      (f_star >= 0) & (f_star <= TARGETS$f_max) &
      (phi1 >= TARGETS$phi1_min) & (phi1 <= TARGETS$phi1_max),
    
    regime = case_when(
      !feasible ~ "infeasible",
      stable & oscillatory ~ "stable_oscillatory",
      stable ~ "stable_real",
      !stable & oscillatory ~ "unstable_oscillatory",
      TRUE ~ "unstable_real"
    ),
    
    score = ifelse(!feasible, Inf,
                   score_band(e_star, TARGETS$e_min, TARGETS$e_max) +
                     score_band(omega_star, TARGETS$omega_min, TARGETS$omega_max) +
                     score_band(d_star, 0, TARGETS$d_max) +
                     score_band(f_star, 0, TARGETS$f_max) +
                     score_band(phi1, TARGETS$phi1_min, TARGETS$phi1_max) +
                     ifelse(!is.na(maxReEig) & maxReEig < 0, 0, 0.5)
    )
  )

# ============================================================
# 4) Hopf boundary in rF for each (psi, phi2) within rh_ok & econ_ok
# ============================================================
hopf_roots <- res_aug %>%
  filter(rh_ok, econ_ok) %>%
  arrange(rF) %>%
  group_by(psi, phi2) %>%
  summarise(
    hopf_rF = {
      Hvec <- H
      rvec <- rF
      ok <- which(!is.na(Hvec) & is.finite(Hvec))
      if (length(ok) < 2) NA_real_ else {
        H2 <- Hvec[ok]; r2 <- rvec[ok]
        idx <- which(sign(H2[-1]) != sign(H2[-length(H2)]))
        if (length(idx) == 0) NA_real_ else {
          ii <- idx[1]
          r1 <- r2[ii]; r2b <- r2[ii+1]
          H1 <- H2[ii]; H2b <- H2[ii+1]
          r1 + (0 - H1) * (r2b - r1) / (H2b - H1)
        }
      }
    },
    .groups = "drop"
  )

# ============================================================
# 5) Best candidates
# ============================================================
best <- res_aug %>%
  filter(rh_ok, econ_ok) %>%
  arrange(score, maxReEig) %>%
  slice_head(n = 50)

# ============================================================
# 6) Export
# ============================================================
write_csv(res,     file.path(out_dir, "results.csv"))
write_csv(res_aug, file.path(out_dir, "results_augmented.csv"))
write_csv(hopf_roots, file.path(out_dir, "hopf_roots.csv"))
write_csv(best,    file.path(out_dir, "best_candidates_top50.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# ============================================================
# 7) Quick plots
# ============================================================
theme_set(theme_minimal(base_size = 12))

p1 <- ggplot(res_aug, aes(x = rF_eff, y = phi2, color = econ_ok)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c("FALSE"="grey70","TRUE"="black")) +
  labs(title = "econ_ok points in (rF_eff, phi2) with sigma fixed + omega-targeted backbone",
       x = "rF_eff", y = "phi2", color = "econ_ok")
ggsave(file.path(plot_dir, "econok_rFeff_vs_phi2.png"), p1, width = 7, height = 5, dpi = 220)

p2 <- res_aug %>%
  count(regime) %>%
  ggplot(aes(x = reorder(regime, -n), y = n)) +
  geom_col() +
  labs(title = "Regime counts", x = "regime", y = "count") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(plot_dir, "regime_counts.png"), p2, width = 7, height = 5, dpi = 220)

if (nrow(hopf_roots %>% filter(!is.na(hopf_rF))) > 0) {
  p3 <- hopf_roots %>%
    filter(!is.na(hopf_rF)) %>%
    ggplot(aes(x = psi, y = phi2, fill = hopf_rF)) +
    geom_tile() +
    labs(title = "Hopf boundary: hopf_rF over (psi, phi2) [sigma fixed + omega-targeted backbone]",
         x = "psi", y = "phi2", fill = "hopf_rF")
  ggsave(file.path(plot_dir, "hopf_tiles_psi_phi2.png"), p3, width = 7, height = 5, dpi = 220)
}

p4 <- res_aug %>%
  filter(feasible) %>%
  ggplot(aes(x = r_star, y = omega_star)) +
  geom_point(alpha = 0.5, size = 0.8) +
  labs(title = "Pinning check: omega* vs r* (sigma fixed, omega targeted)",
       x = "r_star", y = "omega_star")
ggsave(file.path(plot_dir, "pinning_omega_vs_rstar.png"), p4, width = 7, height = 5, dpi = 220)

message("Done. Wrote outputs to: ", out_dir)
message("Key files: results_augmented.csv, hopf_roots.csv, best_candidates_top50.csv, plots/*.png")
