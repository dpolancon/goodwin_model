
# ============================================================
# staged_grid_search_consolidated_v2.R
# Wealth-Goodwin staged grid search (Stage 1-3)
# - Always exports attempted grids + failure reasons
# - Never fails on empty rows (schema preallocation)
# - Conditional progression Stage1 -> Stage2 -> Stage3
# - NEW: kappa_max scan integrated into Stage 1 backbone
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(tibble)
})

# ----------------------------
# Paths
# ----------------------------
base_dir <- "outputs/wealth_goodwin/grid_search"

dir_stage1 <- file.path(base_dir, "stage1_backbone")
dir_stage2 <- file.path(base_dir, "stage2_finance")
dir_stage3 <- file.path(base_dir, "stage3_stability")

dir.create(dir_stage1, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_stage2, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_stage3, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Targets
# ----------------------------
omega_target <- 0.65
omega_band   <- c(0.62, 0.70)

e_target <- 0.94
e_band   <- c(0.88, 0.98)

# ----------------------------
# Minimal grids (as requested)
# ----------------------------
sigma_grid <- c(2.25, 2.35, 2.50)
gn_grid    <- seq(0.05, 0.08, length.out = 4)
i_grid     <- seq(0.025, 0.040, length.out = 4)
delta_grid <- seq(0.03, 0.06,  length.out = 4)

# NEW: scan kappa_max to avoid hard-choking feasibility
kappa_max_grid <- c(0.25, 0.35, 0.45)

# Stage 2 scans (tight)
rF_grid   <- seq(0.02, 0.12, length.out = 9)
psi_grid  <- c(5, 10, 15, 20)
phi2_grid <- seq(0.0, 3.0, length.out = 7)

# Hopf scan along rF around each candidate
hopf_rF_span <- 0.05
hopf_n_grid  <- 31

# ----------------------------
# Base parameters (REPLACE with your sourced params if available)
# IMPORTANT: only kappa_max is scanned; everything else stays the same
# ----------------------------
par_base <- list(
  # kappa(r): logistic bounds + center/steepness
  kappa_min = 0.02,
  kappa_max = 0.25,   # baseline (overridden in scan)
  kappa0    = 0.10,
  kappa1    = 30.0,
  
  # Z(d,f)
  phi3 = 8.0,
  phi4 = 1.0,
  
  # omega dynamics
  phi0  = -0.02,
  alpha = 0.02,
  
  # phi1 bounds (endogenized Stage 2)
  phi1_min = 0.10,
  phi1_max = 5.00
)

# ----------------------------
# Utilities
# ----------------------------
logistic <- function(x) 1 / (1 + exp(-x))

kappa_fun <- function(r, p) {
  p$kappa_min + (p$kappa_max - p$kappa_min) * logistic(p$kappa1 * (r - p$kappa0))
}

kappa_prime <- function(r, p) {
  s <- logistic(p$kappa1 * (r - p$kappa0))
  (p$kappa_max - p$kappa_min) * p$kappa1 * s * (1 - s)
}

kappa_inv <- function(kappa_star, p) {
  y <- (kappa_star - p$kappa_min) / (p$kappa_max - p$kappa_min)
  y <- pmin(pmax(y, 1e-12), 1 - 1e-12)
  p$kappa0 + (1 / p$kappa1) * log(y / (1 - y))
}

Z_fun <- function(d, f, p) {
  logistic(p$phi3 * ((d - 1) + p$phi4 * (f - 1)))
}

lambda_fun <- function(r, rF, psi) logistic(psi * (r - rF))

collapse_reasons <- function(reasons) {
  reasons <- unique(reasons[!is.na(reasons) & reasons != ""])
  if (length(reasons) == 0) "ok" else paste(reasons, collapse = "|")
}

# ============================================================
# Stage 1: Backbone feasibility (with kappa_max scan)
# ============================================================
cat("\n====================\nStage 1: Backbone (with kappa_max scan)\n====================\n")

stage1_grid <- tidyr::crossing(
  sigma     = sigma_grid,
  g_n       = gn_grid,
  i         = i_grid,
  delta     = delta_grid,
  kappa_max = kappa_max_grid
) %>%
  mutate(try_id = row_number())

stage1_schema <- tibble(
  try_id = integer(),
  sigma = double(), g_n = double(), i = double(), delta = double(),
  kappa_min = double(), kappa_max = double(),
  kappa_star = double(), r_star = double(),
  d_star = double(), omega_star = double(),
  backbone_ok = logical(),
  reason = character()
)

stage1_res <- pmap_dfr(stage1_grid, function(sigma, g_n, i, delta, kappa_max, try_id) {
  
  # local params (same model, just scanning kappa_max)
  p <- par_base
  p$kappa_max <- kappa_max
  
  reasons <- character()
  
  kappa_star <- sigma * (g_n + delta)
  if (!is.finite(kappa_star)) reasons <- c(reasons, "nonfinite_kappa_star")
  
  if (is.finite(kappa_star) && (kappa_star <= p$kappa_min || kappa_star >= p$kappa_max)) {
    reasons <- c(reasons, "kappa_star_out_of_bounds")
  }
  
  r_star <- NA_real_
  d_star <- NA_real_
  omega_star <- NA_real_
  
  if (!("kappa_star_out_of_bounds" %in% reasons) && is.finite(kappa_star)) {
    
    r_star <- kappa_inv(kappa_star, p)
    if (!is.finite(r_star)) reasons <- c(reasons, "nonfinite_r_star")
    
    # d* = (kappa* - sigma*r*) / g_n
    if (is.finite(r_star) && is.finite(g_n) && g_n > 0) {
      d_star <- (kappa_star - sigma * r_star) / g_n
    } else {
      reasons <- c(reasons, "bad_gn_for_d_star")
    }
    
    if (!is.finite(d_star)) reasons <- c(reasons, "nonfinite_d_star")
    if (is.finite(d_star) && d_star < 0) reasons <- c(reasons, "d_star_negative")
    
    # omega* = 1 - i d* - sigma r*
    if (is.finite(d_star) && is.finite(r_star)) {
      omega_star <- 1 - i * d_star - sigma * r_star
    }
    
    if (!is.finite(omega_star)) reasons <- c(reasons, "nonfinite_omega_star")
    if (is.finite(omega_star) && (omega_star <= 0 || omega_star >= 1)) reasons <- c(reasons, "omega_outside_0_1")
    if (is.finite(omega_star) && !(omega_star >= omega_band[1] && omega_star <= omega_band[2])) {
      reasons <- c(reasons, "omega_outside_target_band")
    }
  }
  
  backbone_ok <- (length(reasons) == 0)
  
  tibble(
    try_id = try_id,
    sigma = sigma, g_n = g_n, i = i, delta = delta,
    kappa_min = p$kappa_min, kappa_max = p$kappa_max,
    kappa_star = kappa_star,
    r_star = r_star,
    d_star = d_star,
    omega_star = omega_star,
    backbone_ok = backbone_ok,
    reason = collapse_reasons(reasons)
  )
}) %>%
  bind_rows(stage1_schema)

# Counters
n_try <- nrow(stage1_res)
n_backbone_ok <- sum(stage1_res$backbone_ok, na.rm = TRUE)

cat("Stage 1 counts:\n")
cat("  n_try         =", n_try, "\n")
cat("  n_backbone_ok =", n_backbone_ok, "\n")

# Exports (same core deliverables)
write_csv(stage1_res, file.path(dir_stage1, "stage1_backbone.csv"))

stage1_reason_counts <- stage1_res %>%
  separate_rows(reason, sep = "\\|") %>%
  count(reason, sort = TRUE)

write_csv(stage1_reason_counts, file.path(dir_stage1, "failure_reasons.csv"))

# NEW: reason counts by kappa_max
stage1_reason_by_kmax <- stage1_res %>%
  separate_rows(reason, sep = "\\|") %>%
  count(kappa_max, reason, sort = TRUE)

write_csv(stage1_reason_by_kmax, file.path(dir_stage1, "failure_reasons_by_kappa_max.csv"))

# NEW: backbone rate by kappa_max
stage1_backbone_rate <- stage1_res %>%
  group_by(kappa_max) %>%
  summarise(
    n = n(),
    n_ok = sum(backbone_ok, na.rm = TRUE),
    ok_rate = n_ok / n,
    .groups = "drop"
  )

write_csv(stage1_backbone_rate, file.path(dir_stage1, "backbone_rate_by_kappa_max.csv"))

# Plots (same + new)
p1 <- ggplot(stage1_res, aes(x = r_star, y = omega_star, color = backbone_ok)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~kappa_max) +
  labs(title = "Stage 1: omega* vs r* (by kappa_max)", x = "r*", y = "omega*") +
  theme_minimal()

p2 <- ggplot(stage1_res, aes(x = r_star, y = d_star, color = backbone_ok)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~kappa_max) +
  labs(title = "Stage 1: d* vs r* (by kappa_max)", x = "r*", y = "d*") +
  theme_minimal()

p3 <- ggplot(stage1_res, aes(x = omega_star, fill = backbone_ok)) +
  geom_histogram(bins = 40, alpha = 0.7) +
  geom_vline(xintercept = omega_band, linetype = "dashed") +
  facet_wrap(~kappa_max) +
  labs(title = "Stage 1: omega* histogram (by kappa_max)", x = "omega*", y = "count") +
  theme_minimal()

p4 <- ggplot(stage1_res, aes(x = d_star, fill = backbone_ok)) +
  geom_histogram(bins = 40, alpha = 0.7) +
  facet_wrap(~kappa_max) +
  labs(title = "Stage 1: d* histogram (by kappa_max)", x = "d*", y = "count") +
  theme_minimal()

p4b <- ggplot(stage1_backbone_rate, aes(x = factor(kappa_max), y = ok_rate)) +
  geom_col() +
  labs(title = "Stage 1: backbone OK rate by kappa_max", x = "kappa_max", y = "OK rate") +
  theme_minimal()

# Tile plot for binding constraints by kappa_max (top reasons only)
top_reasons <- stage1_reason_counts %>% filter(reason != "ok") %>% slice_head(n = 8) %>% pull(reason)
p4c <- stage1_reason_by_kmax %>%
  filter(reason %in% top_reasons) %>%
  ggplot(aes(x = factor(kappa_max), y = reason, fill = n)) +
  geom_tile() +
  labs(title = "Stage 1: reason counts by kappa_max (top reasons)", x = "kappa_max", y = "reason") +
  theme_minimal()

ggsave(file.path(dir_stage1, "omega_vs_r.png"), p1, width = 10, height = 6, dpi = 160)
ggsave(file.path(dir_stage1, "d_vs_r.png"),     p2, width = 10, height = 6, dpi = 160)
ggsave(file.path(dir_stage1, "hist_omega.png"), p3, width = 10, height = 6, dpi = 160)
ggsave(file.path(dir_stage1, "hist_d.png"),     p4, width = 10, height = 6, dpi = 160)
ggsave(file.path(dir_stage1, "backbone_rate_by_kappa_max.png"), p4b, width = 8, height = 5, dpi = 160)
ggsave(file.path(dir_stage1, "reasons_by_kappa_max_tile.png"),  p4c, width = 10, height = 5, dpi = 160)

# Stage 1 diagnosis (top 3 constraints + minimal tweaks)
cat("\nStage 1 top binding constraints:\n")
print(stage1_reason_counts %>% slice_head(n = 8))

top3 <- stage1_reason_counts %>% slice_head(n = 3) %>% pull(reason)
cat("\nStage 1 minimal tweak suggestions (within same model):\n")
if ("kappa_star_out_of_bounds" %in% top3) {
  cat("  - kappa* out of bounds: kappa_max too tight relative to sigma*(g_n+delta). Raise kappa_max scan range,\n")
  cat("    or narrow (sigma,g_n,delta) to keep kappa* inside bounds.\n")
}
if ("d_star_negative" %in% top3) {
  cat("  - d*<0: you need r* <= kappa*/sigma. If r* is too high for your kappa* region, shift kappa0 or kappa1\n")
  cat("    (same logistic form) so the implied r* is lower for the same kappa*.\n")
}
if ("omega_outside_target_band" %in% top3) {
  cat("  - omega miss: omega* = 1 - i d* - sigma r*. Once d* is feasible, tune omega via sigma/i and by selecting\n")
  cat("    candidates with suitable r* (which is controlled by kappa mapping + kappa*).\n")
}

# ============================================================
# Stage 2: Finance / discipline calibration (conditional on Stage 1 survivors)
# ============================================================
cat("\n====================\nStage 2: Finance/Discipline\n====================\n")

stage1_ok <- stage1_res %>% filter(backbone_ok)

stage2_schema <- tibble(
  try_id = integer(),
  cand_id = integer(),
  sigma = double(), g_n = double(), i = double(), delta = double(),
  kappa_min = double(), kappa_max = double(),
  kappa_star = double(), r_star = double(), d_star = double(), omega_star = double(),
  rF = double(), psi = double(), phi2 = double(),
  lambda_star = double(), f_star = double(), Z_star = double(),
  phi1_endog = double(),
  econ_ok = logical(),
  reason = character()
)

if (nrow(stage1_ok) == 0) {
  cat("Stage 2: no backbone survivors. Exporting empty stage2 schema + diagnosis.\n")
  write_csv(stage2_schema, file.path(dir_stage2, "stage2_finance_scan.csv"))
  write_csv(tibble(reason = "no_stage1_survivors", n = 0), file.path(dir_stage2, "failure_reasons.csv"))
} else {
  
  stage2_grid <- tidyr::crossing(
    cand_id = seq_len(nrow(stage1_ok)),
    rF = rF_grid,
    psi = psi_grid,
    phi2 = phi2_grid
  )
  
  stage2_res <- pmap_dfr(stage2_grid, function(cand_id, rF, psi, phi2) {
    
    row <- stage1_ok[cand_id, ]
    
    # local params (carry kappa_max forward for later stability derivatives)
    p <- par_base
    p$kappa_max <- row$kappa_max
    
    reasons <- character()
    
    r_star <- row$r_star
    d_star <- row$d_star
    g_n    <- row$g_n
    
    lambda_star <- lambda_fun(r_star, rF = rF, psi = psi)
    if (!is.finite(lambda_star)) reasons <- c(reasons, "nonfinite_lambda_star")
    if (is.finite(lambda_star) && (lambda_star < 1e-4 || lambda_star > 1 - 1e-4)) reasons <- c(reasons, "lambda_saturation")
    
    # f* identity:
    # iotaF* = r* * lambda*/(1-lambda*)
    # f*     = iotaF* / g_n
    f_star <- NA_real_
    if (is.finite(lambda_star) && is.finite(r_star) && is.finite(g_n) && g_n > 0) {
      if (lambda_star >= 1) {
        reasons <- c(reasons, "lambda_ge_1")
      } else {
        iotaF_star <- r_star * (lambda_star / (1 - lambda_star))
        f_star <- iotaF_star / g_n
      }
    } else {
      reasons <- c(reasons, "bad_inputs_for_f_star")
    }
    
    if (!is.finite(f_star)) reasons <- c(reasons, "nonfinite_f_star")
    if (is.finite(f_star) && f_star < 0) reasons <- c(reasons, "f_star_negative")
    
    Z_star <- NA_real_
    if (is.finite(d_star) && is.finite(f_star)) {
      Z_star <- Z_fun(d_star, f_star, p)
    } else {
      reasons <- c(reasons, "bad_inputs_for_Z")
    }
    if (!is.finite(Z_star)) reasons <- c(reasons, "nonfinite_Z_star")
    
    # Endogenize phi1 to hit e_target:
    # 0 = phi0 + phi1*e - alpha - phi2*Z  => phi1 = (alpha - phi0 + phi2*Z)/e_target
    phi1_endog <- NA_real_
    if (is.finite(Z_star) && e_target > 0) {
      phi1_endog <- (p$alpha - p$phi0 + phi2 * Z_star) / e_target
    } else {
      reasons <- c(reasons, "bad_inputs_for_phi1")
    }
    
    if (!is.finite(phi1_endog)) reasons <- c(reasons, "nonfinite_phi1")
    if (is.finite(phi1_endog) && (phi1_endog < p$phi1_min || phi1_endog > p$phi1_max)) {
      reasons <- c(reasons, "phi1_out_of_bounds")
    }
    
    econ_ok <- (length(reasons) == 0)
    
    tibble(
      try_id = row$try_id,
      cand_id = cand_id,
      sigma = row$sigma, g_n = row$g_n, i = row$i, delta = row$delta,
      kappa_min = row$kappa_min, kappa_max = row$kappa_max,
      kappa_star = row$kappa_star, r_star = row$r_star, d_star = row$d_star, omega_star = row$omega_star,
      rF = rF, psi = psi, phi2 = phi2,
      lambda_star = lambda_star, f_star = f_star, Z_star = Z_star,
      phi1_endog = phi1_endog,
      econ_ok = econ_ok,
      reason = collapse_reasons(reasons)
    )
  }) %>%
    bind_rows(stage2_schema)
  
  n_stage2_try <- nrow(stage2_res)
  n_econ_ok <- sum(stage2_res$econ_ok, na.rm = TRUE)
  
  cat("Stage 2 counts:\n")
  cat("  n_try     =", n_stage2_try, "\n")
  cat("  n_econ_ok =", n_econ_ok, "\n")
  
  write_csv(stage2_res, file.path(dir_stage2, "stage2_finance_scan.csv"))
  
  stage2_reason_counts <- stage2_res %>%
    separate_rows(reason, sep = "\\|") %>%
    count(reason, sort = TRUE)
  
  write_csv(stage2_reason_counts, file.path(dir_stage2, "failure_reasons.csv"))
  
  candidates_top <- stage2_res %>%
    mutate(omega_gap = abs(omega_star - omega_target),
           lam_mid  = abs(lambda_star - 0.5)) %>%
    arrange(desc(econ_ok), omega_gap, lam_mid) %>%
    slice_head(n = 50)
  
  write_csv(candidates_top, file.path(dir_stage2, "candidates_top.csv"))
  
  # Plots
  p5 <- ggplot(stage2_res, aes(x = rF, y = phi2, color = econ_ok)) +
    geom_point(alpha = 0.6) +
    facet_grid(psi ~ kappa_max) +
    labs(title = "Stage 2: econ_ok cloud in (rF, phi2) by psi and kappa_max", x = "rF", y = "phi2") +
    theme_minimal()
  
  p6 <- ggplot(stage2_res, aes(x = lambda_star, fill = econ_ok)) +
    geom_histogram(bins = 40, alpha = 0.7) +
    facet_wrap(~kappa_max) +
    labs(title = "Stage 2: lambda* distribution (by kappa_max)", x = "lambda*", y = "count") +
    theme_minimal()
  
  p7 <- ggplot(stage2_res, aes(x = lambda_star, y = f_star, color = econ_ok)) +
    geom_point(alpha = 0.6) +
    facet_wrap(~kappa_max) +
    labs(title = "Stage 2: f* vs lambda* (by kappa_max)", x = "lambda*", y = "f*") +
    theme_minimal()
  
  ggsave(file.path(dir_stage2, "econ_ok_cloud.png"), p5, width = 12, height = 8, dpi = 160)
  ggsave(file.path(dir_stage2, "lambda_dist.png"),   p6, width = 10, height = 6, dpi = 160)
  ggsave(file.path(dir_stage2, "f_vs_lambda.png"),   p7, width = 10, height = 6, dpi = 160)
  
  cat("\nStage 2 top binding constraints:\n")
  print(stage2_reason_counts %>% slice_head(n = 10))
}

# ============================================================
# Stage 3: Stability / Hopf (equilibrium-consistent + stricter Hopf)
# - Baseline stability recomputes (lambda*, f*, Z*, phi1_endog)
# - Hopf scan recomputes equilibrium objects for each rFj
# - Strict bracketing: H_left * H_right < 0, away from ~0
# - Requires RH positivity on both sides
# - Requires maxRe eigenvalue sign flip (stable <-> unstable) across bracket
# - Refines root with uniroot(H)
# - Verifies at root: |maxRe| small AND complex pair exists
# - Builds Hopf surface summaries + "3D-ish" ggplots (projection: color/contour)
# - Exports extra plots + LaTeX tables (plain tabular)
# - Prints manifest of all outputs at end
# ============================================================
cat("\n====================\nStage 3: Stability/Hopf (strict boundary)\n====================\n")

# ----------------------------
# Tolerances / sanity gates
# ----------------------------
H_bracket_eps   <- 1e-10   # avoid treating near-zero as a "flip"
eig_re_tol      <- 1e-5    # |maxRe| at root
eig_im_tol      <- 1e-5    # require complex part > tol
maxRe_flip_tol  <- 1e-6    # fuzz in sign flip test
uniroot_tol     <- 1e-8
uniroot_maxiter <- 100

# ----------------------------
# Plain LaTeX tabular writer (no extra packages)
# ----------------------------
latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&#_{}$])", "\\\\\\1", x, perl = TRUE)
  x
}

df_to_latex_tabular <- function(df, caption = NULL, label = NULL, align = NULL) {
  stopifnot(is.data.frame(df))
  if (is.null(align)) align <- paste0(rep("l", ncol(df)), collapse = "")
  header <- paste(latex_escape(names(df)), collapse = " & ")
  body <- apply(df, 1, function(row) paste(latex_escape(row), collapse = " & "))
  body <- paste0(body, " \\\\")
  out <- c(
    "\\begin{table}[H]",
    "\\centering",
    if (!is.null(caption)) paste0("\\caption{", caption, "}"),
    if (!is.null(label))   paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{", align, "}"),
    "\\hline",
    paste0(header, " \\\\"),
    "\\hline",
    body,
    "\\hline",
    "\\end{tabular}",
    "\\end{table}"
  )
  paste(out, collapse = "\n")
}

# ----------------------------
# Schemas
# ----------------------------
stage3_schema <- tibble(
  try_id = integer(), cand_id = integer(),
  sigma = double(), g_n = double(), i = double(), delta = double(),
  kappa_min = double(), kappa_max = double(),
  r_star = double(), d_star = double(), omega_star = double(),
  rF = double(), psi = double(), phi2 = double(),
  r_implied = double(),
  lambda_star = double(), f_star = double(), Z_star = double(),
  phi1_endog = double(),
  c1 = double(), c2 = double(), c3 = double(), H = double(),
  maxReEig = double(), maxImEig = double(),
  stable = logical(),
  RH_ok = logical(),
  has_complex = logical(),
  reason = character()
)

hopf_schema <- tibble(
  try_id = integer(), cand_id = integer(),
  sigma = double(), g_n = double(), i = double(), delta = double(),
  kappa_max = double(),
  psi = double(), phi2 = double(),
  omega_star = double(), d_star = double(),
  rF0 = double(),
  rF_left = double(), rF_right = double(),
  H_left = double(), H_right = double(),
  maxRe_left = double(), maxRe_right = double(),
  stable_left = logical(), stable_right = logical(),
  RH_left = logical(), RH_right = logical(),
  rF_root = double(),
  H_root = double(),
  maxRe_root = double(), maxIm_root = double(),
  lambda_root = double(), f_root = double(), Z_root = double(), phi1_root = double(),
  dmaxRe_drF = double(),
  hopf_ok = logical(),
  reason = character()
)

# ----------------------------
# Guard rails: need stage2_res
# ----------------------------
if (!exists("stage2_res")) {
  
  cat("Stage 3: stage2_res not in env (Stage 2 did not run). Exporting empty outputs.\n")
  write_csv(stage3_schema, file.path(dir_stage3, "stage3_stability.csv"))
  write_csv(hopf_schema,  file.path(dir_stage3, "hopf_roots.csv"))
  write_csv(tibble(reason = "no_stage2_run", n = 0), file.path(dir_stage3, "failure_reasons.csv"))
  
} else {
  
  stage2_ok <- stage2_res %>% filter(econ_ok)
  
  if (nrow(stage2_ok) == 0) {
    
    cat("Stage 3: no econ_ok candidates. Exporting empty outputs.\n")
    write_csv(stage3_schema, file.path(dir_stage3, "stage3_stability.csv"))
    write_csv(hopf_schema,  file.path(dir_stage3, "hopf_roots.csv"))
    write_csv(tibble(reason = "no_stage2_econ_ok", n = 0), file.path(dir_stage3, "failure_reasons.csv"))
    
  } else {
    
    # ------------------------------------------------------------
    # Core helper: recompute equilibrium objects at rFj + compute stability objects
    # ------------------------------------------------------------
    compute_at_rF <- function(row, rFj) {
      
      if (is.data.frame(row)) row <- as.list(row)
      if (!is.list(row))      row <- as.list(row)
      
      p <- par_base
      p$kappa_max <- row$kappa_max
      
      sigma <- row$sigma
      g_n   <- row$g_n
      i     <- row$i
      delta <- row$delta
      psi   <- row$psi
      phi2  <- row$phi2
      
      e_star <- e_target
      omega  <- row$omega_star
      d      <- row$d_star
      
      # implied r from backbone identity
      r <- (1 - omega - i*d) / sigma
      if (!is.finite(r)) return(list(ok=FALSE, reason="nonfinite_r"))
      
      kappa <- kappa_fun(r, p)
      kp    <- kappa_prime(r, p)
      if (!is.finite(kappa) || !is.finite(kp)) return(list(ok=FALSE, reason="nonfinite_kappa_or_kp"))
      
      g  <- kappa / sigma - delta
      gr <- kp / sigma
      
      # lambda at rFj
      lam   <- lambda_fun(r, rFj, psi)
      denom <- 1 - lam
      if (!is.finite(lam) || !is.finite(denom) || abs(denom) < 1e-12) {
        return(list(ok=FALSE, reason="lambda_denominator_near_zero_or_nonfinite"))
      }
      
      # f* identity: iotaF = r * lam/(1-lam), f = iotaF/g_n
      q     <- lam / denom
      iotaF <- r * q
      f     <- iotaF / g_n
      if (!is.finite(f)) return(list(ok=FALSE, reason="nonfinite_f_star"))
      
      # Z(d,f)
      Z <- Z_fun(d, f, p)
      if (!is.finite(Z)) return(list(ok=FALSE, reason="nonfinite_Z_star"))
      
      # endogenize phi1 at THIS rFj to hit e_target
      phi1 <- (p$alpha - p$phi0 + phi2 * Z) / e_target
      if (!is.finite(phi1)) return(list(ok=FALSE, reason="nonfinite_phi1_endog"))
      
      # derivatives for Jacobian
      lam_r   <- psi * lam * (1 - lam)     # d lambda / d r
      q_r     <- lam_r / (denom^2)
      iotaF_r <- q + r*q_r
      f_r     <- iotaF_r / g_n
      if (!is.finite(f_r)) return(list(ok=FALSE, reason="nonfinite_f_r"))
      
      Z_s <- Z * (1 - Z)
      
      dr_domega <- -1 / sigma
      dr_dd     <- -i / sigma
      
      # s = phi3[(d-1)+phi4(f-1)]
      ds_domega <- p$phi3 * (p$phi4 * f_r * dr_domega)
      ds_dd     <- p$phi3 * (1 + p$phi4 * f_r * dr_dd)
      
      Z_omega <- Z_s * ds_domega
      Z_d     <- Z_s * ds_dd
      
      # Jacobian for (e, omega, d)
      J11 <- (g - g_n)
      J12 <- e_star * gr * dr_domega
      J13 <- e_star * gr * dr_dd
      
      J21 <- omega * phi1
      J22 <- omega * (-phi2 * Z_omega)
      J23 <- omega * (-phi2 * Z_d)
      
      J31 <- 0
      J32 <- 1 + kp * dr_domega - d * gr * dr_domega
      J33 <- i - g + kp * dr_dd - d * gr * dr_dd
      
      J <- matrix(c(J11,J12,J13,
                    J21,J22,J23,
                    J31,J32,J33), nrow=3, byrow=TRUE)
      
      if (any(!is.finite(J))) return(list(ok=FALSE, reason="nonfinite_jacobian"))
      
      ev <- tryCatch(eigen(J)$values, error=function(e) NA_complex_)
      if (all(is.na(ev))) return(list(ok=FALSE, reason="eigen_failed"))
      
      maxRe <- max(Re(ev))
      maxIm <- max(abs(Im(ev)))
      
      # RH coefficients for det(λI - J) = λ^3 + c1 λ^2 + c2 λ + c3
      c1 <- -sum(diag(J))
      c2 <- (J11*J22 - J12*J21) + (J11*J33 - J13*J31) + (J22*J33 - J23*J32)
      c3 <- -det(J)
      H  <- c1*c2 - c3
      
      RH_ok <- is.finite(c1) && is.finite(c2) && is.finite(c3) && is.finite(H) &&
        (c1 > 0 && c2 > 0 && c3 > 0)
      
      list(
        ok=TRUE,
        r=r, lam=lam, f=f, Z=Z, phi1=phi1,
        c1=c1, c2=c2, c3=c3, H=H, RH_ok=RH_ok,
        maxRe=maxRe, maxIm=maxIm,
        stable = (is.finite(maxRe) && maxRe < 0),
        has_complex = (is.finite(maxIm) && maxIm > eig_im_tol)
      )
    }
    
    # ------------------------------------------------------------
    # Baseline stage3: compute stability for each econ_ok candidate
    # ------------------------------------------------------------
    row_lists <- split(stage2_ok, seq_len(nrow(stage2_ok)))
    
    stage3_res <- purrr::map_dfr(row_lists, function(df_row) {
      
      row <- as.list(df_row)
      out <- compute_at_rF(row, rFj = row$rF)
      
      if (!isTRUE(out$ok)) {
        return(tibble(
          try_id = row$try_id, cand_id = row$cand_id,
          sigma = row$sigma, g_n = row$g_n, i = row$i, delta = row$delta,
          kappa_min = row$kappa_min, kappa_max = row$kappa_max,
          r_star = row$r_star, d_star = row$d_star, omega_star = row$omega_star,
          rF = row$rF, psi = row$psi, phi2 = row$phi2,
          r_implied = NA_real_,
          lambda_star = NA_real_, f_star = NA_real_, Z_star = NA_real_,
          phi1_endog = NA_real_,
          c1 = NA_real_, c2 = NA_real_, c3 = NA_real_, H = NA_real_,
          maxReEig = NA_real_, maxImEig = NA_real_,
          stable = NA, RH_ok = NA, has_complex = NA,
          reason = out$reason
        ))
      }
      
      tibble(
        try_id = row$try_id, cand_id = row$cand_id,
        sigma = row$sigma, g_n = row$g_n, i = row$i, delta = row$delta,
        kappa_min = row$kappa_min, kappa_max = row$kappa_max,
        r_star = row$r_star, d_star = row$d_star, omega_star = row$omega_star,
        rF = row$rF, psi = row$psi, phi2 = row$phi2,
        r_implied = out$r,
        lambda_star = out$lam, f_star = out$f, Z_star = out$Z,
        phi1_endog = out$phi1,
        c1 = out$c1, c2 = out$c2, c3 = out$c3, H = out$H,
        maxReEig = out$maxRe, maxImEig = out$maxIm,
        stable = out$stable, RH_ok = out$RH_ok, has_complex = out$has_complex,
        reason = "ok"
      )
    }) %>%
      bind_rows(stage3_schema)
    
    write_csv(stage3_res, file.path(dir_stage3, "stage3_stability.csv"))
    
    stage3_reason_counts <- stage3_res %>%
      separate_rows(reason, sep="\\|") %>%
      count(reason, sort=TRUE)
    
    write_csv(stage3_reason_counts, file.path(dir_stage3, "failure_reasons.csv"))
    
    cat("Stage 3 counts:\n")
    cat("  n_try        =", nrow(stage3_res), "\n")
    cat("  n_stable     =", sum(stage3_res$stable, na.rm=TRUE), "\n")
    cat("  n_RH_ok      =", sum(stage3_res$RH_ok,  na.rm=TRUE), "\n")
    cat("  n_hasComplex =", sum(stage3_res$has_complex, na.rm=TRUE), "\n")
    
    # ------------------------------------------------------------
    # Hopf root detection: strict + uniroot + verification
    # ------------------------------------------------------------
    hopf_rows <- list()
    
    stage3_scan <- stage3_res %>%
      filter(reason == "ok") %>%
      filter(is.finite(H), is.finite(maxReEig))
    
    for (k in seq_len(nrow(stage3_scan))) {
      
      row_df <- stage3_scan[k, ]
      row <- as.list(row_df)
      
      rF_seq <- seq(row$rF - hopf_rF_span, row$rF + hopf_rF_span, length.out = hopf_n_grid)
      
      H_seq     <- rep(NA_real_, length(rF_seq))
      maxRe_seq <- rep(NA_real_, length(rF_seq))
      maxIm_seq <- rep(NA_real_, length(rF_seq))
      RH_seq    <- rep(NA, length(rF_seq))
      st_seq    <- rep(NA, length(rF_seq))
      cx_seq    <- rep(NA, length(rF_seq))
      
      for (j in seq_along(rF_seq)) {
        rFj <- rF_seq[j]
        tmp <- compute_at_rF(row, rFj = rFj)
        if (!isTRUE(tmp$ok)) next
        H_seq[j]     <- tmp$H
        maxRe_seq[j] <- tmp$maxRe
        maxIm_seq[j] <- tmp$maxIm
        RH_seq[j]    <- tmp$RH_ok
        st_seq[j]    <- tmp$stable
        cx_seq[j]    <- tmp$has_complex
      }
      
      for (m in seq_len(length(rF_seq) - 1)) {
        
        H1 <- H_seq[m];  H2 <- H_seq[m+1]
        if (!is.finite(H1) || !is.finite(H2)) next
        if (abs(H1) < H_bracket_eps || abs(H2) < H_bracket_eps) next
        if (H1 * H2 >= 0) next
        
        if (!(isTRUE(RH_seq[m]) && isTRUE(RH_seq[m+1]))) next
        
        R1 <- maxRe_seq[m]; R2 <- maxRe_seq[m+1]
        if (!is.finite(R1) || !is.finite(R2)) next
        
        flip_ok <- ((R1 < -maxRe_flip_tol && R2 >  maxRe_flip_tol) ||
                      (R2 < -maxRe_flip_tol && R1 >  maxRe_flip_tol) ||
                      (abs(R1) <= maxRe_flip_tol || abs(R2) <= maxRe_flip_tol))
        if (!flip_ok) next
        
        x1 <- rF_seq[m]; x2 <- rF_seq[m+1]
        
        H_of_rF <- function(x) {
          tmp <- compute_at_rF(row, rFj = x)
          if (!isTRUE(tmp$ok)) return(NA_real_)
          tmp$H
        }
        
        root <- tryCatch(
          uniroot(H_of_rF, lower = x1, upper = x2, tol = uniroot_tol, maxiter = uniroot_maxiter),
          error = function(e) NULL
        )
        if (is.null(root) || !is.finite(root$root)) next
        rF_root <- root$root
        
        tmp_root <- compute_at_rF(row, rFj = rF_root)
        if (!isTRUE(tmp_root$ok)) next
        
        # verification at root
        hopf_ok <- TRUE
        if (!(is.finite(tmp_root$maxRe) && abs(tmp_root$maxRe) <= eig_re_tol)) hopf_ok <- FALSE
        if (!isTRUE(tmp_root$has_complex)) hopf_ok <- FALSE
        
        # transversality proxy: finite difference slope d(maxRe)/drF around root
        h <- 1e-4
        tmp_a <- compute_at_rF(row, rFj = rF_root - h)
        tmp_b <- compute_at_rF(row, rFj = rF_root + h)
        dmaxRe_drF <- NA_real_
        if (isTRUE(tmp_a$ok) && isTRUE(tmp_b$ok) && is.finite(tmp_a$maxRe) && is.finite(tmp_b$maxRe)) {
          dmaxRe_drF <- (tmp_b$maxRe - tmp_a$maxRe) / (2*h)
        }
        
        hopf_rows[[length(hopf_rows) + 1]] <- tibble(
          try_id = row$try_id, cand_id = row$cand_id,
          sigma = row$sigma, g_n = row$g_n, i = row$i, delta = row$delta,
          kappa_max = row$kappa_max,
          psi = row$psi, phi2 = row$phi2,
          omega_star = row$omega_star, d_star = row$d_star,
          rF0 = row$rF,
          rF_left = x1, rF_right = x2,
          H_left = H1, H_right = H2,
          maxRe_left = R1, maxRe_right = R2,
          stable_left = st_seq[m], stable_right = st_seq[m+1],
          RH_left = RH_seq[m], RH_right = RH_seq[m+1],
          rF_root = rF_root,
          H_root = tmp_root$H,
          maxRe_root = tmp_root$maxRe, maxIm_root = tmp_root$maxIm,
          lambda_root = tmp_root$lam, f_root = tmp_root$f, Z_root = tmp_root$Z, phi1_root = tmp_root$phi1,
          dmaxRe_drF = dmaxRe_drF,
          hopf_ok = hopf_ok,
          reason = if (hopf_ok) "ok" else "failed_root_verification"
        )
      }
    }
    
    hopf_out <- if (length(hopf_rows) == 0) hopf_schema else bind_rows(hopf_rows, hopf_schema)
    write_csv(hopf_out, file.path(dir_stage3, "hopf_roots.csv"))
    
    cat("\nStage 3 Hopf roots found (strict):\n")
    cat("  total rows    =", sum(is.finite(hopf_out$rF_root), na.rm=TRUE), "\n")
    cat("  verified Hopf =", sum(is.finite(hopf_out$rF_root) & hopf_out$hopf_ok, na.rm=TRUE), "\n")
    
    # convenience
    hopf_verified <- hopf_out %>% filter(is.finite(rF_root), hopf_ok)
    
    # ------------------------------------------------------------
    # Plots: baseline tiles + boundary
    # ------------------------------------------------------------
    p8 <- ggplot(stage3_res, aes(x = rF, y = phi2, fill = stable)) +
      geom_tile() +
      facet_grid(psi ~ kappa_max) +
      labs(title = "Stage 3: stability tiles (by psi and kappa_max)", x = "rF", y = "phi2") +
      theme_minimal()
    ggsave(file.path(dir_stage3, "stability_tiles.png"), p8, width = 12, height = 8, dpi = 160)
    
    p8b <- ggplot(stage3_res, aes(x = rF, y = phi2, fill = RH_ok)) +
      geom_tile() +
      facet_grid(psi ~ kappa_max) +
      labs(title = "Stage 3: RH positivity tiles (c1,c2,c3>0)", x = "rF", y = "phi2") +
      theme_minimal()
    ggsave(file.path(dir_stage3, "RH_ok_tiles.png"), p8b, width = 12, height = 8, dpi = 160)
    
    p8c <- ggplot(stage3_res, aes(x = rF, y = phi2, fill = maxReEig)) +
      geom_tile() +
      facet_grid(psi ~ kappa_max) +
      labs(title = "Stage 3: maxRe(eigs) heatmap", x = "rF", y = "phi2") +
      theme_minimal()
    ggsave(file.path(dir_stage3, "maxRe_tiles.png"), p8c, width = 12, height = 8, dpi = 160)
    
    p9 <- stage3_res %>%
      mutate(H_sign = case_when(
        !is.finite(H) ~ NA_character_,
        H > 0 ~ "H>0",
        H < 0 ~ "H<0",
        TRUE ~ "H=0"
      )) %>%
      ggplot(aes(x = rF, y = phi2, fill = H_sign)) +
      geom_tile() +
      facet_grid(psi ~ kappa_max) +
      labs(title = "Stage 3: sign(H) tiles", x = "rF", y = "phi2") +
      theme_minimal()
    ggsave(file.path(dir_stage3, "hopf_tiles.png"), p9, width = 12, height = 8, dpi = 160)
    
    p10 <- hopf_verified %>%
      ggplot(aes(x = rF_root, y = phi2)) +
      geom_point(alpha = 0.7) +
      facet_grid(psi ~ kappa_max) +
      labs(title = "Stage 3: verified Hopf boundary points in (rF_root, phi2)", x = "rF_root", y = "phi2") +
      theme_minimal()
    ggsave(file.path(dir_stage3, "hopf_boundary_rFroot_phi2.png"), p10, width = 12, height = 8, dpi = 160)
    
    p10b <- hopf_verified %>%
      ggplot(aes(x = rF_root)) +
      geom_histogram(bins = 40, alpha = 0.7) +
      labs(title = "Stage 3: verified rF_root distribution", x = "rF_root", y = "count") +
      theme_minimal()
    ggsave(file.path(dir_stage3, "hopf_root_hist.png"), p10b, width = 10, height = 5, dpi = 160)
    
    p10c <- hopf_verified %>%
      ggplot(aes(x = rF_root, y = kappa_max)) +
      geom_point(alpha = 0.7) +
      labs(title = "Stage 3: verified Hopf roots (rF_root vs kappa_max)", x = "rF_root", y = "kappa_max") +
      theme_minimal()
    ggsave(file.path(dir_stage3, "hopf_roots_scatter.png"), p10c, width = 10, height = 5, dpi = 160)
    
    p10d <- hopf_verified %>%
      ggplot(aes(x = phi2, y = psi, color = dmaxRe_drF)) +
      geom_point(alpha = 0.8) +
      facet_grid(sigma ~ kappa_max) +
      labs(title = "Transversality proxy: d(maxRe)/drF at verified Hopf roots", x = "phi2", y = "psi", color = "slope") +
      theme_minimal()
    ggsave(file.path(dir_stage3, "hopf_transversality_map.png"), p10d, width = 12, height = 8, dpi = 160)
    
    # ------------------------------------------------------------
    # Hopf surface construction (median rF_root by grid cell)
    # ------------------------------------------------------------
    # The surface is naturally on (sigma, kappa_max, psi, phi2) -> rF_root
    hopf_surface <- hopf_verified %>%
      group_by(sigma, kappa_max, psi, phi2) %>%
      summarise(
        n_roots = n(),
        rF_root_med = median(rF_root, na.rm = TRUE),
        rF_root_p25 = quantile(rF_root, 0.25, na.rm = TRUE),
        rF_root_p75 = quantile(rF_root, 0.75, na.rm = TRUE),
        .groups = "drop"
      )
    
    write_csv(hopf_surface, file.path(dir_stage3, "hopf_surface.csv"))
    
    # Surface heatmap tiles
    pS1 <- hopf_surface %>%
      ggplot(aes(x = phi2, y = psi, fill = rF_root_med)) +
      geom_tile() +
      facet_grid(sigma ~ kappa_max) +
      labs(title = "Hopf surface (median rF_root): tiles over (phi2, psi)", x = "phi2", y = "psi", fill = "rF_root") +
      theme_minimal()
    ggsave(file.path(dir_stage3, "hopf_surface_tiles.png"), pS1, width = 12, height = 8, dpi = 160)
    
    # Count tiles: where surface is well-supported
    pS2 <- hopf_surface %>%
      ggplot(aes(x = phi2, y = psi, fill = n_roots)) +
      geom_tile() +
      facet_grid(sigma ~ kappa_max) +
      labs(title = "Hopf surface support: count of verified roots per cell", x = "phi2", y = "psi", fill = "n") +
      theme_minimal()
    ggsave(file.path(dir_stage3, "hopf_surface_counts.png"), pS2, width = 12, height = 8, dpi = 160)
    
    # ------------------------------------------------------------
    # "3D" ggplots (projection): z encoded as color + contours + size
    # ------------------------------------------------------------
    p_surface_cloud <- hopf_surface %>%
      ggplot(aes(x = phi2, y = psi)) +
      geom_point(aes(color = rF_root_med, size = n_roots), alpha = 0.8) +
      facet_grid(sigma ~ kappa_max) +
      labs(
        title = "Hopf surface cloud (projection): color = median rF_root, size = n_roots",
        x = "phi2",
        y = "psi",
        color = "rF_root",
        size = "n"
      ) +
      theme_minimal()
    ggsave(file.path(dir_stage3, "hopf_surface_cloud.png"), p_surface_cloud, width = 12, height = 8, dpi = 160)
    
    p_surface_contour <- hopf_surface %>%
      ggplot(aes(x = phi2, y = psi)) +
      geom_tile(aes(fill = rF_root_med)) +
      geom_contour(aes(z = rF_root_med), bins = 12, alpha = 0.8) +
      facet_grid(sigma ~ kappa_max) +
      labs(
        title = "Hopf surface contours: median rF_root over (phi2, psi)",
        x = "phi2",
        y = "psi",
        fill = "rF_root"
      ) +
      theme_minimal()
    ggsave(file.path(dir_stage3, "hopf_surface_contours.png"), p_surface_contour, width = 12, height = 8, dpi = 160)
    
    # One PNG per (kappa_max, sigma) slice for close reading
    slice_keys <- hopf_surface %>%
      distinct(kappa_max, sigma) %>%
      arrange(kappa_max, sigma)
    
    for (kk in seq_len(nrow(slice_keys))) {
      km <- slice_keys$kappa_max[kk]
      sg <- slice_keys$sigma[kk]
      
      df_slice <- hopf_surface %>% filter(kappa_max == km, sigma == sg)
      if (nrow(df_slice) == 0) next
      
      p_slice <- ggplot(df_slice, aes(x = phi2, y = psi)) +
        geom_tile(aes(fill = rF_root_med)) +
        geom_contour(aes(z = rF_root_med), bins = 12, alpha = 0.8) +
        geom_point(aes(size = n_roots), alpha = 0.6) +
        labs(
          title = paste0("Hopf surface slice: kappa_max=", km, ", sigma=", sg),
          x = "phi2",
          y = "psi",
          fill = "rF_root",
          size = "n"
        ) +
        theme_minimal()
      
      ggsave(
        filename = file.path(dir_stage3, paste0("hopf_surface_slice_kmax_", km, "_sigma_", sg, ".png")),
        plot = p_slice,
        width = 8, height = 6, dpi = 160
      )
    }
    
    # ------------------------------------------------------------
    # Tables: Stage 3 summary + Hopf summary
    # ------------------------------------------------------------
    tab_stage3 <- stage3_res %>%
      filter(reason == "ok") %>%
      group_by(psi, kappa_max) %>%
      summarise(
        n = n(),
        n_stable = sum(stable, na.rm=TRUE),
        share_stable = n_stable / n,
        n_RH_ok = sum(RH_ok, na.rm=TRUE),
        share_RH_ok = n_RH_ok / n,
        .groups = "drop"
      ) %>%
      arrange(psi, kappa_max)
    
    tex_stage3 <- df_to_latex_tabular(
      tab_stage3,
      caption = "Stage 3 summary by $(\\psi,\\kappa_{\\max})$.",
      label   = "tab:stage3_summary",
      align   = "rrrrrr"
    )
    writeLines(tex_stage3, file.path(dir_stage3, "table_stage3_summary.tex"))
    write_csv(tab_stage3, file.path(dir_stage3, "table_stage3_summary.csv"))
    
    tab_hopf <- hopf_verified %>%
      group_by(psi, kappa_max) %>%
      summarise(
        n_verified = n(),
        rF_root_med = median(rF_root, na.rm=TRUE),
        rF_root_p25 = quantile(rF_root, 0.25, na.rm=TRUE),
        rF_root_p75 = quantile(rF_root, 0.75, na.rm=TRUE),
        slope_med  = median(dmaxRe_drF, na.rm=TRUE),
        .groups = "drop"
      ) %>%
      arrange(psi, kappa_max)
    
    tex_hopf <- df_to_latex_tabular(
      tab_hopf,
      caption = "Verified Hopf root summary by $(\\psi,\\kappa_{\\max})$ (median/quantiles).",
      label   = "tab:hopf_summary",
      align   = "rrrrrr"
    )
    writeLines(tex_hopf, file.path(dir_stage3, "table_hopf_summary.tex"))
    write_csv(tab_hopf, file.path(dir_stage3, "table_hopf_summary.csv"))
    
    # ------------------------------------------------------------
    # Manifest: list every output under base_dir
    # ------------------------------------------------------------
    manifest <- tibble(
      path = list.files(base_dir, recursive = TRUE, full.names = TRUE)
    ) %>%
      mutate(
        bytes = file.size(path),
        ext   = tools::file_ext(path)
      ) %>%
      arrange(path)
    
    write_csv(manifest, file.path(base_dir, "manifest_outputs.csv"))
    
    cat("\n====================\nOUTPUT MANIFEST (base_dir)\n====================\n")
    print(manifest)
    cat("\nSaved manifest to:\n  ", file.path(base_dir, "manifest_outputs.csv"), "\n", sep="")
    cat("Saved Stage 3 LaTeX tables to:\n",
        "  ", file.path(dir_stage3, "table_stage3_summary.tex"), "\n",
        "  ", file.path(dir_stage3, "table_hopf_summary.tex"), "\n", sep="")
  }
}

# ============================================================
# Stage 4 + Stage 5 (Consolidated)
# - Reads Stage 3 outputs: stage3_stability.csv, hopf_roots.csv
# - Stage 4: soft scoring + (optional) gating + shortlists + plots + LaTeX
# - Stage 5: builds sim plan around Hopf roots, runs ODE sims (USER HOOK),
#           computes simple cycle diagnostics, exports shortlist + plots + LaTeX
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(purrr)
})

cat("\n====================\nStage 4+5: start\n====================\n")

# ----------------------------
# Paths
# ----------------------------
base_dir   <- "outputs/wealth_goodwin/grid_search"
dir_stage3 <- file.path(base_dir, "stage3_stability")
dir_stage4 <- file.path(base_dir, "stage4_scoring")
dir_stage5 <- file.path(base_dir, "stage5_simulation")

dir.create(dir_stage4, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_stage5, showWarnings = FALSE, recursive = TRUE)

stage3_file <- file.path(dir_stage3, "stage3_stability.csv")
hopf_file   <- file.path(dir_stage3, "hopf_roots.csv")

stopifnot(file.exists(stage3_file), file.exists(hopf_file))

# ----------------------------
# Helpers
# ----------------------------
to_bool <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), NA, x))
  if (is.numeric(x)) return(ifelse(is.na(x), NA, x != 0))
  if (is.character(x)) {
    y <- tolower(trimws(x))
    return(ifelse(y %in% c("true","t","1","yes","y"), TRUE,
                  ifelse(y %in% c("false","f","0","no","n"), FALSE, NA)))
  }
  rep(NA, length(x))
}

# Count peaks: simple local maxima count after smoothing-free rule
count_peaks <- function(y) {
  y <- as.numeric(y)
  if (length(y) < 5 || all(!is.finite(y))) return(NA_integer_)
  y <- y[is.finite(y)]
  if (length(y) < 5) return(NA_integer_)
  dy <- diff(y)
  s  <- sign(dy)
  # peak when slope goes + to -
  sum(s[-length(s)] > 0 & s[-1] < 0, na.rm = TRUE)
}

# Basic amplitude
amp <- function(y) {
  y <- as.numeric(y)
  if (all(!is.finite(y))) return(NA_real_)
  y <- y[is.finite(y)]
  if (length(y) == 0) return(NA_real_)
  max(y) - min(y)
}

# Safe extraction of rF0 from hopf file (NO case_when; avoids "object rF not found")
add_rF0 <- function(hopf_df) {
  hopf2 <- hopf_df
  if ("rF0" %in% names(hopf2)) {
    hopf2 <- hopf2 %>% mutate(rF0 = suppressWarnings(as.numeric(rF0)))
  } else if ("rF" %in% names(hopf2)) {
    hopf2 <- hopf2 %>% mutate(rF0 = suppressWarnings(as.numeric(rF)))
  } else {
    hopf2 <- hopf2 %>% mutate(rF0 = NA_real_)
  }
  hopf2
}

# ----------------------------
# Read Stage 3 inputs
# ----------------------------
stage3 <- readr::read_csv(stage3_file, show_col_types = FALSE)
hopf   <- readr::read_csv(hopf_file,   show_col_types = FALSE)

stage3_ok <- stage3 %>%
  filter(reason == "ok") %>%
  mutate(
    RH_ok       = to_bool(RH_ok),
    stable      = to_bool(stable),
    has_complex = to_bool(has_complex)
  ) %>%
  filter(is.finite(try_id), is.finite(cand_id), is.finite(rF)) %>%
  filter(is.finite(omega_star), is.finite(sigma), is.finite(lambda_star))

hopf2 <- hopf %>%
  mutate(
    hopf_ok = to_bool(hopf_ok),
    rF_root = suppressWarnings(as.numeric(rF_root))
  ) %>%
  add_rF0()

hopf_verified <- hopf2 %>%
  filter(is.finite(rF_root), is.finite(rF0)) %>%
  filter((hopf_ok == TRUE) | (reason == "ok"))

# Aggregate hopf info per candidate (try_id,cand_id,rF0)
hopf_by_cand <- hopf_verified %>%
  group_by(try_id, cand_id, rF0) %>%
  summarise(
    n_roots     = n(),
    rF_root_med = median(rF_root, na.rm = TRUE),
    rF_root_p25 = quantile(rF_root, 0.25, na.rm = TRUE),
    rF_root_p75 = quantile(rF_root, 0.75, na.rm = TRUE),
    maxIm_root_med = if ("maxIm_root" %in% names(.)) median(maxIm_root, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    hopf_complex = is.finite(maxIm_root_med) & (maxIm_root_med > 1e-5)
  )

df <- stage3_ok %>%
  left_join(hopf_by_cand, by = c("try_id","cand_id", "rF" = "rF0")) %>%
  mutate(
    has_hopf = is.finite(rF_root_med) & (n_roots > 0)
  )

cat("Input counts:\n")
cat("  stage3_ok rows     =", nrow(stage3_ok), "\n")
cat("  hopf_verified rows =", nrow(hopf_verified), "\n")
cat("  merged rows        =", nrow(df), "\n")

# ============================================================
# STAGE 4: Scoring + Shortlists
# ============================================================
cat("\n====================\nStage 4: scoring\n====================\n")

# ---- Targets / tolerances (edit here) ----
omega_pref    <- c(0.62, 0.65)
omega_hard_hi <- 0.65
omega_hard_lo <- 0.60
omega_center  <- 0.64
gate_omega_max <- 0.67

sigma_pref_max <- 2.35
sigma_soft_max <- 2.50

lambda_center <- 0.50
lambda_soft   <- c(0.20, 0.80)

rF_root_pref  <- c(0.01, 0.12)
rF_root_soft  <- c(0.00, 0.15)

# ---- Gates (strict by default, but we export relaxed too) ----
gate_RH_required      <- TRUE
gate_complex_required <- TRUE
gate_hopf_required    <- FALSE  # keep FALSE unless you want ONLY hopf-verified

# ---- Weights (edit here) ----
w_omega  <- 8.0
w_sigma  <- 1.0
w_lambda <- 1.2
w_hopf   <- 2.5
w_dist   <- 2.0
w_RH     <- 50.0
w_cx     <- 25.0

df_scored <- df %>%
  mutate(
    # omega penalty: asymmetric, punish >0.65 hard
    pen_omega = case_when(
      omega_star >= omega_pref[1] & omega_star <= omega_pref[2] ~ 0,
      omega_star > omega_hard_hi ~ 5 + 50*(omega_star - omega_hard_hi),
      omega_star < omega_hard_lo ~ 2 + 20*(omega_hard_lo - omega_star),
      TRUE ~ abs(omega_star - omega_center)
    ),
    
    pen_sigma = case_when(
      sigma <= sigma_pref_max ~ 0,
      sigma <= sigma_soft_max ~ (sigma - sigma_pref_max),
      TRUE ~ 2 + 5*(sigma - sigma_soft_max)
    ),
    
    pen_lambda = case_when(
      lambda_star >= lambda_soft[1] & lambda_star <= lambda_soft[2] ~ abs(lambda_star - lambda_center),
      TRUE ~ 2 + 5*abs(lambda_star - lambda_center)
    ),
    
    pen_hopf = case_when(
      has_hopf & rF_root_med >= rF_root_pref[1] & rF_root_med <= rF_root_pref[2] ~ 0,
      has_hopf & rF_root_med >= rF_root_soft[1] & rF_root_med <= rF_root_soft[2] ~ abs(rF_root_med - median(rF_root_pref)),
      has_hopf ~ 2 + 5*abs(rF_root_med - median(rF_root_pref)),
      TRUE ~ 10
    ),
    
    pen_dist = case_when(
      has_hopf ~ abs(rF - rF_root_med),
      TRUE ~ 1
    ),
    
    # IMPORTANT: vector-safe booleans (NO isTRUE() here)
    pen_RH = ifelse(RH_ok == TRUE, 0, 1),
    pen_cx = ifelse(has_complex == TRUE, 0, 1),
    
    score = w_omega*pen_omega + w_sigma*pen_sigma + w_lambda*pen_lambda +
      w_hopf*pen_hopf + w_dist*pen_dist + w_RH*pen_RH + w_cx*pen_cx
  ) %>%
  arrange(score)

write_csv(df_scored, file.path(dir_stage4, "stage4_scored_candidates.csv"))

# Diagnostics: how many pass each gate individually
gate_diag <- tibble(
  n_total = nrow(df_scored),
  n_RH_ok = sum(df_scored$RH_ok == TRUE, na.rm = TRUE),
  n_complex = sum(df_scored$has_complex == TRUE, na.rm = TRUE),
  n_hopf = sum(df_scored$has_hopf == TRUE, na.rm = TRUE),
  n_omega_le_cap = sum(df_scored$omega_star <= gate_omega_max, na.rm = TRUE)
)
write_csv(gate_diag, file.path(dir_stage4, "gate_diagnostics.csv"))
print(gate_diag)

# Strict gated shortlist
df_gated_strict <- df_scored %>%
  filter(if (gate_RH_required) RH_ok == TRUE else TRUE) %>%
  filter(if (gate_complex_required) has_complex == TRUE else TRUE) %>%
  filter(if (gate_hopf_required) has_hopf == TRUE else TRUE) %>%
  filter(is.finite(omega_star) & omega_star <= gate_omega_max)

# Relaxed gated shortlist (fallback): only omega cap + hopf preferred in score
df_gated_relaxed <- df_scored %>%
  filter(is.finite(omega_star) & omega_star <= gate_omega_max)

write_csv(df_gated_strict,  file.path(dir_stage4, "stage4_scored_candidates_GATED_STRICT.csv"))
write_csv(df_gated_relaxed, file.path(dir_stage4, "stage4_scored_candidates_GATED_RELAXED.csv"))

cat("Stage 4 gated sizes:\n")
cat("  strict  =", nrow(df_gated_strict), "\n")
cat("  relaxed =", nrow(df_gated_relaxed), "\n")

pick_cols <- c(
  "try_id","cand_id",
  "sigma","g_n","i","delta","kappa_max",
  "psi","phi2","rF",
  "omega_star","d_star",
  "lambda_star","f_star","Z_star",
  "stable","RH_ok","has_complex",
  "has_hopf","n_roots","rF_root_med","rF_root_p25","rF_root_p75",
  "pen_omega","pen_sigma","pen_lambda","pen_hopf","pen_dist",
  "score"
)

short_top50_all <- df_scored %>% select(any_of(pick_cols)) %>% slice_head(n = 50)
short_top50_rel <- df_gated_relaxed %>% select(any_of(pick_cols)) %>% slice_head(n = 50)
short_top50_str <- df_gated_strict %>% select(any_of(pick_cols)) %>% slice_head(n = 50)

write_csv(short_top50_all, file.path(dir_stage4, "stage4_shortlist_top50_ALL.csv"))
write_csv(short_top50_rel, file.path(dir_stage4, "stage4_shortlist_top50_GATED_RELAXED.csv"))
write_csv(short_top50_str, file.path(dir_stage4, "stage4_shortlist_top50_GATED_STRICT.csv"))

# Plots (simple ggplot + ggsave)
df_plot <- df_scored %>% mutate(score_decile = ntile(score, 10) %>% as.factor())

p1 <- ggplot(df_plot, aes(x = sigma, y = omega_star, color = score_decile)) +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = omega_pref, linetype = "dashed") +
  labs(title = "Stage 4: omega* vs sigma (score deciles)", x = "sigma", y = "omega*") +
  theme_minimal()
ggsave(file.path(dir_stage4, "score_scatter_omega_sigma.png"), p1, width = 10, height = 6, dpi = 160)

p2 <- ggplot(df_plot %>% filter(has_hopf), aes(x = rF_root_med, y = phi2, color = score_decile)) +
  geom_point(alpha = 0.75) +
  facet_grid(psi ~ kappa_max) +
  geom_vline(xintercept = rF_root_pref, linetype = "dashed") +
  labs(title = "Stage 4: Hopf boundary locations (rF_root_med vs phi2)",
       x = "rF_root (median)", y = "phi2") +
  theme_minimal()
ggsave(file.path(dir_stage4, "score_hopf_boundary_phi2.png"), p2, width = 12, height = 8, dpi = 160)

p3 <- ggplot(df_plot, aes(x = score)) +
  geom_histogram(bins = 50, alpha = 0.75) +
  labs(title = "Stage 4: score distribution", x = "score (lower=better)", y = "count") +
  theme_minimal()
ggsave(file.path(dir_stage4, "score_hist.png"), p3, width = 10, height = 5, dpi = 160)

p4 <- ggplot(df_plot, aes(x = omega_star)) +
  geom_histogram(bins = 40, alpha = 0.75) +
  geom_vline(xintercept = omega_pref, linetype = "dashed") +
  labs(title = "Stage 4: omega* distribution", x = "omega*", y = "count") +
  theme_minimal()
ggsave(file.path(dir_stage4, "omega_hist.png"), p4, width = 10, height = 5, dpi = 160)

# LaTeX table (top 20, relaxed gated) rounded to 2 decimals
# Requires: \usepackage{booktabs} and if using [H], \usepackage{float}
short_tex <- df_gated_relaxed %>%
  select(try_id, cand_id, psi, kappa_max, sigma, omega_star, rF, has_hopf, rF_root_med, score) %>%
  slice_head(n = 20) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  mutate(has_hopf = ifelse(has_hopf == TRUE, "T", "F"))

cols  <- names(short_tex)
align <- paste(rep("l", length(cols)), collapse = "")
header <- paste(cols, collapse = " & ")
rows <- apply(short_tex, 1, function(r) paste(r, collapse = " & "))
rows <- paste0(rows, " \\\\")

tex_lines <- c(
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Stage 4 shortlist (top 20 by score; rounded to 2 decimals).}",
  "\\label{tab:stage4_shortlist_top20}",
  paste0("\\begin{tabular}{", align, "}"),
  "\\toprule",
  paste0(header, " \\\\"),
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(tex_lines, file.path(dir_stage4, "shortlist_top20.tex"))

# ============================================================
# STAGE 5: Hopf-centered simulation sweep + cycle diagnostics
# ============================================================
cat("\n====================\nStage 5: simulation sweep\n====================\n")

# ---- SIM SETTINGS (edit here) ----
sim_t_end <- 800
sim_dt    <- 0.05
burn_frac <- 0.50   # fraction of time to drop as burn-in
rF_eps    <- 0.005  # simulate around Hopf boundary: root +/- eps
n_plan_candidates <- 30  # how many candidates (from Stage 4) to plan sims for

# ---- USER HOOK: define your simulation function here ----
# Must return a data.frame/tibble with columns: time, e, omega, d
# And it must accept (row, rF_sim, t_end, dt).
#
# Example skeleton:
# simulate_model <- function(row, rF_sim, t_end, dt) {
#   # ... use deSolve::ode or your own integrator ...
#   # return tibble(time=..., e=..., omega=..., d=...)
# }
#
# For now we stop if not defined.
if (!exists("simulate_model")) {
  cat("\nERROR: You must define simulate_model(row, rF_sim, t_end, dt) before Stage 5 can run.\n")
  cat("Create it above this block (or source it), then re-run.\n\n")
} else {
  
  # Candidate pool: use relaxed gated (omega cap) and prioritize hopf existence in score
  cand_pool <- df_gated_relaxed %>%
    arrange(score) %>%
    slice_head(n = n_plan_candidates) %>%
    mutate(pool_rank = row_number())
  
  # Build simulation plan:
  # - If hopf exists: rF = root-eps, root, root+eps
  # - Else: rF = baseline rF (just to see behavior)
  plan <- cand_pool %>%
    transmute(
      pool_rank,
      try_id, cand_id,
      psi, kappa_max, sigma, g_n, i, delta,
      phi2,
      omega_star, d_star, lambda_star,
      rF_base = rF,
      has_hopf, rF_root_med
    ) %>%
    mutate(
      rF_root_med = ifelse(has_hopf == TRUE, rF_root_med, NA_real_)
    ) %>%
    group_by(pool_rank, try_id, cand_id) %>%
    group_modify(~{
      x <- .x[1,]
      if (isTRUE(x$has_hopf) && is.finite(x$rF_root_med)) {
        bind_rows(
          mutate(x, plan_case="root_minus", rF_sim = x$rF_root_med - rF_eps),
          mutate(x, plan_case="root",       rF_sim = x$rF_root_med),
          mutate(x, plan_case="root_plus",  rF_sim = x$rF_root_med + rF_eps)
        )
      } else {
        mutate(x, plan_case="baseline", rF_sim = x$rF_base)
      }
    }) %>%
    ungroup() %>%
    mutate(sim_id = row_number())
  
  write_csv(plan, file.path(dir_stage5, "stage5_sim_plan.csv"))
  cat("Stage 5 plan rows =", nrow(plan), "\n")
  
  # Run sims + diagnostics
  sim_metrics <- purrr::map_dfr(seq_len(nrow(plan)), function(j) {
    
    rowj <- plan[j, ]
    out <- tryCatch(
      simulate_model(row = rowj, rF_sim = rowj$rF_sim, t_end = sim_t_end, dt = sim_dt),
      error = function(e) NULL
    )
    
    if (is.null(out)) {
      return(tibble(
        sim_id = rowj$sim_id,
        try_id = rowj$try_id, cand_id = rowj$cand_id,
        plan_case = rowj$plan_case,
        rF_sim = rowj$rF_sim,
        ok = FALSE,
        reason = "simulate_model_error",
        amp_e = NA_real_, amp_omega = NA_real_, amp_d = NA_real_,
        peaks_omega = NA_integer_, peaks_e = NA_integer_,
        mean_omega = NA_real_
      ))
    }
    
    out <- as_tibble(out)
    if (!all(c("time","e","omega","d") %in% names(out))) {
      return(tibble(
        sim_id = rowj$sim_id,
        try_id = rowj$try_id, cand_id = rowj$cand_id,
        plan_case = rowj$plan_case,
        rF_sim = rowj$rF_sim,
        ok = FALSE,
        reason = "missing_required_cols_time_e_omega_d",
        amp_e = NA_real_, amp_omega = NA_real_, amp_d = NA_real_,
        peaks_omega = NA_integer_, peaks_e = NA_integer_,
        mean_omega = NA_real_
      ))
    }
    
    # burn-in drop
    t_cut <- min(out$time, na.rm = TRUE) + burn_frac * (max(out$time, na.rm = TRUE) - min(out$time, na.rm = TRUE))
    out2 <- out %>% filter(time >= t_cut)
    
    tibble(
      sim_id = rowj$sim_id,
      try_id = rowj$try_id, cand_id = rowj$cand_id,
      plan_case = rowj$plan_case,
      rF_sim = rowj$rF_sim,
      ok = TRUE,
      reason = "ok",
      amp_e     = amp(out2$e),
      amp_omega = amp(out2$omega),
      amp_d     = amp(out2$d),
      peaks_omega = count_peaks(out2$omega),
      peaks_e     = count_peaks(out2$e),
      mean_omega  = mean(out2$omega, na.rm = TRUE)
    )
  })
  
  write_csv(sim_metrics, file.path(dir_stage5, "stage5_sim_metrics.csv"))
  
  # Cycle shortlist criteria (edit thresholds here)
  # - want some oscillation: omega amplitude not tiny + several peaks
  omega_amp_min <- 0.01
  peaks_min     <- 3
  omega_mean_hi <- 0.70  # keep mean omega from going insane
  
  sim_short <- sim_metrics %>%
    filter(ok == TRUE) %>%
    filter(is.finite(amp_omega), is.finite(peaks_omega), is.finite(mean_omega)) %>%
    mutate(
      cycle_like = (amp_omega >= omega_amp_min) & (peaks_omega >= peaks_min) & (mean_omega <= omega_mean_hi)
    ) %>%
    arrange(desc(cycle_like), desc(amp_omega), desc(peaks_omega))
  
  write_csv(sim_short, file.path(dir_stage5, "stage5_cycle_shortlist.csv"))
  
  # Plots
  pM1 <- ggplot(sim_short, aes(x = amp_omega)) +
    geom_histogram(bins = 40, alpha = 0.75) +
    labs(title = "Stage 5: omega amplitude distribution (post burn-in)", x = "amp(omega)", y = "count") +
    theme_minimal()
  ggsave(file.path(dir_stage5, "omega_amplitude_hist.png"), pM1, width = 10, height = 5, dpi = 160)
  
  pM2 <- ggplot(sim_short, aes(x = peaks_omega, y = amp_omega)) +
    geom_point(alpha = 0.75) +
    labs(title = "Stage 5: cycle strength (peaks vs amplitude)", x = "peaks(omega)", y = "amp(omega)") +
    theme_minimal()
  ggsave(file.path(dir_stage5, "peaks_vs_amp_omega.png"), pM2, width = 10, height = 6, dpi = 160)
  
  pM3 <- ggplot(sim_short, aes(x = rF_sim, y = amp_omega)) +
    geom_point(alpha = 0.75) +
    facet_wrap(~ plan_case) +
    labs(title = "Stage 5: amp(omega) by rF_sim (cases)", x = "rF_sim", y = "amp(omega)") +
    theme_minimal()
  ggsave(file.path(dir_stage5, "amp_by_rF_case.png"), pM3, width = 12, height = 6, dpi = 160)
  
  # LaTeX table: top 15 cycle-like runs
  # Requires: \usepackage{booktabs} and if using [H], \usepackage{float}
  top_tex <- sim_short %>%
    filter(cycle_like == TRUE) %>%
    slice_head(n = 15) %>%
    select(sim_id, try_id, cand_id, plan_case, rF_sim, amp_omega, peaks_omega, mean_omega) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  cols2  <- names(top_tex)
  align2 <- paste(rep("l", length(cols2)), collapse = "")
  header2 <- paste(cols2, collapse = " & ")
  rows2 <- apply(top_tex, 1, function(r) paste(r, collapse = " & "))
  rows2 <- paste0(rows2, " \\\\")
  
  tex2 <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{Stage 5 cycle-like simulations (top 15).}",
    "\\label{tab:stage5_cycle_shortlist}",
    paste0("\\begin{tabular}{", align2, "}"),
    "\\toprule",
    paste0(header2, " \\\\"),
    "\\midrule",
    rows2,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(tex2, file.path(dir_stage5, "table_stage5_cycle_shortlist.tex"))
  
  # Manifest (stage4+stage5)
  manifest45 <- tibble(
    path = list.files(base_dir, recursive = TRUE, full.names = TRUE)
  ) %>%
    mutate(bytes = file.size(path), ext = tools::file_ext(path)) %>%
    arrange(path)
  
  write_csv(manifest45, file.path(base_dir, "manifest_outputs_stage4_5.csv"))
  
  cat("\nStage 5 done.\n")
  cat("  plan:    ", file.path(dir_stage5, "stage5_sim_plan.csv"), "\n", sep="")
  cat("  metrics: ", file.path(dir_stage5, "stage5_sim_metrics.csv"), "\n", sep="")
  cat("  shortlist:", file.path(dir_stage5, "stage5_cycle_shortlist.csv"), "\n", sep="")
  cat("  latex:   ", file.path(dir_stage5, "table_stage5_cycle_shortlist.tex"), "\n", sep="")
}

# Stage 4 manifest
manifest4 <- tibble(path = list.files(dir_stage4, full.names = TRUE)) %>%
  mutate(bytes = file.size(path), ext = tools::file_ext(path)) %>%
  arrange(path)
write_csv(manifest4, file.path(dir_stage4, "manifest_stage4.csv"))

cat("\n====================\nStage 4+5: complete\n====================\n")



# ============================================================
# 05_stage5_run.R
# Stage 5A: certify Hopf (refine root + transversality)
# Stage 5B: time-domain sims just below/above root + cycle metrics
# Stage 5C: robustness ring (small jitter) + survival
# Stage 5D: final ranking + exports
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(deSolve)
})

# ----------------------------
# 0) Paths + I/O
# ----------------------------
base_dir <- "outputs/wealth_goodwin/grid_search"  # <-- change if needed
stage5_dir <- file.path(base_dir, "stage5")
dir.create(stage5_dir, showWarnings = FALSE, recursive = TRUE)

dir_5A <- file.path(stage5_dir, "stage5A_hopf_cert")
dir_5B <- file.path(stage5_dir, "stage5B_sim")
dir_5C <- file.path(stage5_dir, "stage5C_robust")
dir_5D <- file.path(stage5_dir, "stage5D_final")
dir.create(dir_5A, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_5B, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_5C, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_5D, showWarnings = FALSE, recursive = TRUE)

# Candidate source (try these in order)
cand_candidates <- c(
  file.path(base_dir, "stage4_5_partition", "track_H_hopf_candidates.csv"),
  file.path(base_dir, "stage4_scoring", "shortlist_top20_GATED.csv"),
  file.path(base_dir, "stage4_scoring", "shortlist_top20.csv"),
  file.path(base_dir, "grid_results_augmented.csv"),
  file.path(base_dir, "grid_results.csv")
)
cand_file <- cand_candidates[file.exists(cand_candidates)][1]
stopifnot(!is.na(cand_file))
message("Stage 5 using candidates from: ", cand_file)

cand <- readr::read_csv(cand_file, show_col_types = FALSE)

# We require at least an id + a hopf seed (root guess). Name-flexible.
# If your column names differ, map them below.
map_id_col   <- c("id","ID","row_id","candidate_id")
map_rF_col   <- c("rF_root","rF_star","rF_hopf","rF_guess","rF")

id_col <- map_id_col[map_id_col %in% names(cand)][1]
rF_col <- map_rF_col[map_rF_col %in% names(cand)][1]

stopifnot(!is.na(id_col), !is.na(rF_col))

# Optional steady-state columns if you already computed them
map_e_col     <- c("e_star","e*","e_fp","e")
map_om_col    <- c("omega_star","omega*","omega_fp","omega")
map_d_col     <- c("d_star","d*","d_fp","d")

e_col  <- map_e_col[map_e_col %in% names(cand)][1]
om_col <- map_om_col[map_om_col %in% names(cand)][1]
d_col  <- map_d_col[map_d_col %in% names(cand)][1]

# ----------------------------
# 1) YOU plug model functions here (minimal hooks)
# ----------------------------
# You need these three functions defined somewhere in your project:
#   (i)  get_par_from_row(row) -> named numeric vector/list of parameters
#   (ii) fixed_point(par)      -> named numeric state c(e=?, omega=?, d=?)
#   (iii) jacobian(par, xstar) -> 3x3 numeric matrix at fixed point
#
# Either source them:
# source("R/model_core.R")
#
# Or define wrappers below that call your existing code.

get_par_from_row <- function(row) {
  # Minimal default: take all numeric columns except obvious non-params
  # If your pipeline already has a clean "par_*" naming convention, adjust here.
  drop <- c(id_col)
  out <- row %>%
    select(where(is.numeric)) %>%
    select(-any_of(drop)) %>%
    as.list()
  out <- out[!vapply(out, is.null, logical(1))]
  out
}

fixed_point <- function(par) {
  stop("Define fixed_point(par): return c(e=..., omega=..., d=...)")
}

jacobian <- function(par, xstar) {
  stop("Define jacobian(par, xstar): return 3x3 matrix")
}

ode_rhs <- function(t, state, par) {
  stop("Define ode_rhs(t, state, par): return list(c(de, domega, dd))")
}

# ----------------------------
# 2) Stage 5 knobs
# ----------------------------
eps_root   <- 0.002      # +/- around root for sim
bracket_w  <- 0.02       # bracket width for uniroot around seed
t_end      <- 800
dt         <- 0.05
burn_in    <- 300
perturb    <- c(e = 0.002, omega = -0.002, d = 0.02)

# Robustness ring
robust_reps <- 30
jitter_sd   <- 0.02      # 2% multiplicative jitter
jitter_cols <- c("phi2","kappa_max","lambda")  # change to your actual param names
robust_t_end <- 400

# Cycle detection thresholds
cycle_var_min <- 1e-5
amp_min       <- 1e-3

# ----------------------------
# 3) Helpers (small + boring, like reality)
# ----------------------------
eig_maxRe <- function(J) {
  ev <- eigen(J, only.values = TRUE)$values
  max(Re(ev))
}

hopf_indicator <- function(par, rF_value) {
  # assumes your parameter list uses rF (or r_F) naming; adjust if needed
  par2 <- par
  if ("rF" %in% names(par2)) par2[["rF"]] <- rF_value
  else if ("r_F" %in% names(par2)) par2[["r_F"]] <- rF_value
  else if ("bar_r_F" %in% names(par2)) par2[["bar_r_F"]] <- rF_value
  else par2[["rF"]] <- rF_value
  
  xstar <- fixed_point(par2)
  J <- jacobian(par2, xstar)
  eig_maxRe(J)
}

refine_root <- function(par, rF_seed, width = bracket_w) {
  a <- max(1e-6, rF_seed - width)
  b <- rF_seed + width
  fa <- hopf_indicator(par, a)
  fb <- hopf_indicator(par, b)
  
  # If no sign change, expand bracket a bit (not forever).
  k <- 0
  while (fa * fb > 0 && k < 6) {
    width <- width * 1.5
    a <- max(1e-6, rF_seed - width)
    b <- rF_seed + width
    fa <- hopf_indicator(par, a)
    fb <- hopf_indicator(par, b)
    k <- k + 1
  }
  
  if (fa * fb > 0) {
    return(list(ok = FALSE, reason = "no_bracket_sign_change", a = a, b = b, fa = fa, fb = fb, rF_root = NA_real_))
  }
  
  rr <- uniroot(function(x) hopf_indicator(par, x), interval = c(a, b))
  list(ok = TRUE, reason = "ok", a = a, b = b, fa = fa, fb = fb, rF_root = rr$root)
}

transversality_fd <- function(par, rF_root, h = 1e-4) {
  # d/d rF of maxRe eigenvalue at root (finite difference)
  f1 <- hopf_indicator(par, rF_root - h)
  f2 <- hopf_indicator(par, rF_root + h)
  (f2 - f1) / (2 * h)
}

simulate_one <- function(par, state0, t_end, dt) {
  times <- seq(0, t_end, by = dt)
  out <- deSolve::ode(
    y = state0,
    times = times,
    func = ode_rhs,
    parms = par,
    method = "lsoda"
  )
  as_tibble(out)
}

cycle_metrics <- function(sim, burn_in = burn_in, var_min = cycle_var_min, amp_min = amp_min) {
  sim2 <- sim %>% filter(time >= burn_in)
  
  # choose omega as the oscillator proxy (change if you prefer e or d)
  if (!("omega" %in% names(sim2))) return(tibble(has_cycle = FALSE, reason = "missing_omega"))
  
  x <- sim2$omega
  t <- sim2$time
  
  v <- var(x, na.rm = TRUE)
  amp <- (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) / 2
  
  if (!is.finite(v) || !is.finite(amp)) return(tibble(has_cycle = FALSE, reason = "nonfinite"))
  if (v < var_min || amp < amp_min) return(tibble(has_cycle = FALSE, reason = "low_variance_or_amp", var = v, amp = amp))
  
  # crude period: average spacing of local maxima
  dx <- diff(x)
  # local maxima indices: slope goes + to -
  peaks <- which(head(dx, -1) > 0 & tail(dx, -1) <= 0) + 1
  if (length(peaks) < 3) return(tibble(has_cycle = TRUE, reason = "few_peaks", var = v, amp = amp, period = NA_real_))
  
  periods <- diff(t[peaks])
  tibble(
    has_cycle = TRUE,
    reason = "ok",
    var = v,
    amp = amp,
    period = mean(periods, na.rm = TRUE),
    omega_mean = mean(x, na.rm = TRUE),
    omega_min = min(x, na.rm = TRUE),
    omega_max = max(x, na.rm = TRUE)
  )
}

bounded_ok <- function(sim) {
  # sanity bounds (adjust if your model allows beyond these)
  ok <- TRUE
  if ("e" %in% names(sim)) ok <- ok && all(is.finite(sim$e)) && all(sim$e >= 0) && all(sim$e <= 1.5)
  if ("omega" %in% names(sim)) ok <- ok && all(is.finite(sim$omega)) && all(sim$omega >= 0) && all(sim$omega <= 1.5)
  if ("d" %in% names(sim)) ok <- ok && all(is.finite(sim$d)) && all(sim$d >= -0.1) && all(sim$d <= 10)
  ok
}

# ----------------------------
# 4) Stage 5A: certify Hopf candidates
# ----------------------------
cert <- cand %>%
  mutate(
    candidate_id = .data[[id_col]],
    rF_seed = .data[[rF_col]]
  ) %>%
  select(candidate_id, rF_seed, everything()) %>%
  mutate(
    cert = pmap(., function(...) {
      row <- tibble(...)
      par <- get_par_from_row(row)
      
      # If steady state was stored, you can skip fixed_point in the indicator,
      # but we keep it unified and robust.
      res <- refine_root(par, rF_seed)
      if (!res$ok) return(tibble(ok = FALSE, reason = res$reason, rF_root = NA_real_, bracket_a = res$a, bracket_b = res$b))
      
      tr <- transversality_fd(par, res$rF_root)
      tibble(ok = TRUE, reason = "ok", rF_root = res$rF_root, bracket_a = res$a, bracket_b = res$b, transversality = tr)
    })
  ) %>%
  unnest(cert)

cert_ok  <- cert %>% filter(ok)
cert_bad <- cert %>% filter(!ok)

write_csv(cert_ok,  file.path(dir_5A, "hopf_certified.csv"))
write_csv(cert_bad, file.path(dir_5A, "hopf_rejected.csv"))

message("Stage 5A done: certified=", nrow(cert_ok), " rejected=", nrow(cert_bad))

# ----------------------------
# 5) Stage 5B: simulate below/above root + cycle metrics
# ----------------------------
sim_results <- cert_ok %>%
  mutate(
    sim = pmap(., function(...) {
      row <- tibble(...)
      par <- get_par_from_row(row)
      
      rF_root <- row$rF_root
      rF_lo <- max(1e-6, rF_root - eps_root)
      rF_hi <- rF_root + eps_root
      
      set_rF <- function(par, val) {
        if ("rF" %in% names(par)) par[["rF"]] <- val
        else if ("r_F" %in% names(par)) par[["r_F"]] <- val
        else if ("bar_r_F" %in% names(par)) par[["bar_r_F"]] <- val
        else par[["rF"]] <- val
        par
      }
      
      par_lo <- set_rF(par, rF_lo)
      par_hi <- set_rF(par, rF_hi)
      
      # initial state: fixed point + perturb
      xstar <- fixed_point(par)  # at seed rF if embedded in par
      # better: use root-specific, so override:
      xstar_lo <- fixed_point(par_lo)
      x0 <- xstar_lo + perturb
      
      sim_lo <- simulate_one(par_lo, x0, t_end = t_end, dt = dt)
      sim_hi <- simulate_one(par_hi, x0, t_end = t_end, dt = dt)
      
      met_lo <- cycle_metrics(sim_lo)
      met_hi <- cycle_metrics(sim_hi)
      
      ok_lo <- bounded_ok(sim_lo)
      ok_hi <- bounded_ok(sim_hi)
      
      # Save trajectory panels if you want (basic CSV export instead of plots here)
      id <- as.character(row$candidate_id)
      write_csv(sim_lo, file.path(dir_5B, paste0("traj_", id, "_below.csv")))
      write_csv(sim_hi, file.path(dir_5B, paste0("traj_", id, "_above.csv")))
      
      tibble(
        bounded_below = ok_lo,
        bounded_above = ok_hi,
        has_cycle_below = met_lo$has_cycle %||% FALSE,
        has_cycle_above = met_hi$has_cycle %||% FALSE,
        amp_below = met_lo$amp %||% NA_real_,
        amp_above = met_hi$amp %||% NA_real_,
        period_below = met_lo$period %||% NA_real_,
        period_above = met_hi$period %||% NA_real_,
        omega_mean_above = met_hi$omega_mean %||% NA_real_
      )
    })
  ) %>%
  unnest(sim)

write_csv(sim_results, file.path(dir_5B, "stage5B_cycle_metrics.csv"))
message("Stage 5B done: wrote trajectories + metrics.")

# ----------------------------
# 6) Stage 5C: robustness ring (jitter) on certified + cycle-above cases
# ----------------------------
targets <- sim_results %>%
  filter(bounded_above, has_cycle_above)

robust <- targets %>%
  mutate(
    robust = pmap(., function(...) {
      row <- tibble(...)
      par0 <- get_par_from_row(row)
      
      # set rF slightly above root for robustness testing
      rF_test <- row$rF_root + eps_root
      if ("rF" %in% names(par0)) par0[["rF"]] <- rF_test
      else if ("r_F" %in% names(par0)) par0[["r_F"]] <- rF_test
      else if ("bar_r_F" %in% names(par0)) par0[["bar_r_F"]] <- rF_test
      else par0[["rF"]] <- rF_test
      
      # base IC: fixed point at test rF + perturb
      xstar <- fixed_point(par0)
      x0 <- xstar + perturb
      
      reps <- seq_len(robust_reps)
      out <- map_dfr(reps, function(r) {
        par <- par0
        
        # multiplicative jitter on selected parameters (only if present)
        for (nm in jitter_cols) {
          if (nm %in% names(par)) {
            par[[nm]] <- par[[nm]] * exp(rnorm(1, mean = 0, sd = jitter_sd))
          }
        }
        
        sim <- simulate_one(par, x0, t_end = robust_t_end, dt = dt)
        okb <- bounded_ok(sim)
        met <- cycle_metrics(sim, burn_in = min(burn_in, robust_t_end/2))
        
        tibble(
          rep = r,
          bounded = okb,
          has_cycle = (met$has_cycle %||% FALSE),
          amp = met$amp %||% NA_real_,
          period = met$period %||% NA_real_
        )
      })
      
      tibble(
        survival_rate = mean(out$bounded & out$has_cycle, na.rm = TRUE),
        bounded_rate  = mean(out$bounded, na.rm = TRUE),
        cycle_rate    = mean(out$has_cycle, na.rm = TRUE),
        amp_med       = median(out$amp, na.rm = TRUE),
        period_med    = median(out$period, na.rm = TRUE)
      )
    })
  ) %>%
  unnest(robust)

write_csv(robust, file.path(dir_5C, "robustness_summary.csv"))
message("Stage 5C done: robustness summary exported.")

# ----------------------------
# 7) Stage 5D: rank + export top candidates
# ----------------------------
final <- sim_results %>%
  left_join(robust %>% select(candidate_id, survival_rate, bounded_rate, cycle_rate, amp_med, period_med),
            by = "candidate_id") %>%
  mutate(
    survival_rate = coalesce(survival_rate, 0),
    bounded_rate  = coalesce(bounded_rate, 0),
    cycle_rate    = coalesce(cycle_rate, 0),
    # simple score: prioritize robust cycles + boundedness
    score_stage5 = 3*survival_rate + 1*bounded_rate + 1*cycle_rate +
      0.2*if_else(is.finite(amp_above), pmin(amp_above, 1), 0)
  ) %>%
  arrange(desc(score_stage5))

write_csv(final, file.path(dir_5D, "stage5_final_ranked.csv"))
write_csv(final %>% slice_head(n = 20), file.path(dir_5D, "stage5_top20.csv"))

message("Stage 5D done: ranked outputs written to ", dir_5D)

# ------------------------------------------------------------
# tiny helper: null-coalescing (because R loves suffering)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

