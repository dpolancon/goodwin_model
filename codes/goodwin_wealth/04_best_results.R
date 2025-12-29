# ============================================================
# best_results_summary.R  (POST: summarize CSV outputs)
# Reads outputs/wealth_goodwin/grid_search/best_results/*.csv
# Writes summary tables + plots to .../best_results/summary/
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

base_dir <- "outputs/wealth_goodwin/grid_search/best_results"
sum_dir  <- file.path(base_dir, "summary")
dir.create(sum_dir, showWarnings = FALSE, recursive = TRUE)

# ---- load
res_file  <- file.path(base_dir, "results_augmented.csv")
hopf_file <- file.path(base_dir, "hopf_roots.csv")
best_file <- file.path(base_dir, "best_candidates_top50.csv")

stopifnot(file.exists(res_file))
res  <- readr::read_csv(res_file, show_col_types = FALSE)
hopf <- if (file.exists(hopf_file)) readr::read_csv(hopf_file, show_col_types = FALSE) else tibble()
best <- if (file.exists(best_file)) readr::read_csv(best_file, show_col_types = FALSE) else tibble()

# ---- 1) headline stats
headline <- tibble(
  n_total     = nrow(res),
  n_feasible  = sum(res$feasible %in% TRUE, na.rm = TRUE),
  n_rh_ok     = sum(res$rh_ok %in% TRUE, na.rm = TRUE),
  n_econ_ok   = sum(res$econ_ok %in% TRUE, na.rm = TRUE),
  share_feas  = n_feasible / n_total,
  share_rh_ok = n_rh_ok / n_total,
  share_econ  = n_econ_ok / n_total
)
write_csv(headline, file.path(sum_dir, "headline.csv"))
print(headline)

# ---- 2) regime counts
reg_counts <- res %>%
  count(regime, sort = TRUE) %>%
  mutate(share = n / sum(n))
write_csv(reg_counts, file.path(sum_dir, "regime_counts.csv"))
print(reg_counts)

# ---- 3) key ranges among econ_ok
ranges_econ <- res %>%
  filter(econ_ok) %>%
  summarise(
    n = n(),
    across(c(rF, rF_eff, psi, phi2, phi1, e_star, omega_star, d_star, f_star, maxReEig, maxImEig, H, period),
           list(min = ~suppressWarnings(min(.x, na.rm = TRUE)),
                med = ~suppressWarnings(median(.x, na.rm = TRUE)),
                max = ~suppressWarnings(max(.x, na.rm = TRUE))),
           .names = "{.col}_{.fn}")
  )
write_csv(ranges_econ, file.path(sum_dir, "ranges_econ_ok.csv"))
print(ranges_econ)

# ---- 4) shortlist: top 20 by score + closeness to Hopf (|H|)
top20 <- res %>%
  filter(econ_ok, rh_ok) %>%
  mutate(absH = abs(H)) %>%
  arrange(score, absH, maxReEig) %>%
  slice_head(n = 20) %>%
  select(sigma, psi, rF, rF_eff, phi2, phi1, e_star, omega_star, d_star, f_star,
         stable, oscillatory, maxReEig, maxImEig, H, absH, period, score)
write_csv(top20, file.path(sum_dir, "top20_candidates.csv"))
print(top20)

# ---- 5) “near Hopf” subset (small |H|) for inspection
near_hopf <- res %>%
  filter(econ_ok, rh_ok, is.finite(H)) %>%
  mutate(absH = abs(H)) %>%
  arrange(absH) %>%
  slice_head(n = 100)
write_csv(near_hopf, file.path(sum_dir, "near_hopf_top100.csv"))

# ---- 6) plots
theme_set(theme_minimal(base_size = 12))

# (A) econ_ok region map
pA <- ggplot(res, aes(x = rF_eff, y = phi2, color = econ_ok)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c("FALSE"="grey70","TRUE"="black")) +
  labs(title = "econ_ok region in (rF_eff, phi2)", x = "rF_eff", y = "phi2", color = "econ_ok")
ggsave(file.path(sum_dir, "econok_cloud_rFeff_phi2.png"), pA, width = 7, height = 5, dpi = 220)

# (B) stability types among econ_ok
pB <- res %>%
  filter(econ_ok) %>%
  ggplot(aes(x = rF_eff, y = phi2, color = regime)) +
  geom_point(alpha = 0.7, size = 0.9) +
  labs(title = "Regimes within econ_ok", x = "rF_eff", y = "phi2", color = "regime")
ggsave(file.path(sum_dir, "regimes_within_econok.png"), pB, width = 7, height = 5, dpi = 220)

# (C) abs(H) heat (proxy Hopf closeness) over psi x phi2 (econ_ok only)
pC <- res %>%
  filter(econ_ok, rh_ok) %>%
  mutate(absH = abs(H)) %>%
  ggplot(aes(x = psi, y = phi2, fill = absH)) +
  geom_tile() +
  labs(title = "Hopf proximity: |H| over (psi, phi2) [econ_ok & rh_ok]", x = "psi", y = "phi2", fill = "|H|")
ggsave(file.path(sum_dir, "absH_tiles_psi_phi2.png"), pC, width = 7, height = 5, dpi = 220)

# (D) period vs phi2 (stable oscillatory only)
pD <- res %>%
  filter(econ_ok, stable, oscillatory) %>%
  ggplot(aes(x = phi2, y = period)) +
  geom_point(alpha = 0.6, size = 0.9) +
  labs(title = "Implied period (2π/Im λ) in stable oscillatory region", x = "phi2", y = "period")
ggsave(file.path(sum_dir, "period_vs_phi2_stable_osc.png"), pD, width = 7, height = 5, dpi = 220)

# (E) optional: hopf_roots tile (if present)
if (nrow(hopf) > 0 && all(c("psi","phi2","hopf_rF") %in% names(hopf))) {
  pE <- hopf %>%
    filter(!is.na(hopf_rF)) %>%
    ggplot(aes(x = psi, y = phi2, fill = hopf_rF)) +
    geom_tile() +
    labs(title = "Hopf boundary: hopf_rF over (psi, phi2)", x = "psi", y = "phi2", fill = "hopf_rF")
  ggsave(file.path(sum_dir, "hopf_tiles_from_hopf_roots.png"), pE, width = 7, height = 5, dpi = 220)
}

message("Summary written to: ", sum_dir)
message("Key tables: headline.csv, regime_counts.csv, ranges_econ_ok.csv, top20_candidates.csv")
message("Key plots: econok_cloud_rFeff_phi2.png, regimes_within_econok.png, absH_tiles_psi_phi2.png")
