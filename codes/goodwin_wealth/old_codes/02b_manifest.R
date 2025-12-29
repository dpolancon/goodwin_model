suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(fs)
  library(ggplot2)
})

OUT_DIR <- fs::path("outputs", "wealth_goodwin", "grid_search")

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

# ---- audit function ----
audit_plane_df <- function(df, plane_name) {
  
  dfF <- df %>% filter(ok, feasible_RH, is.finite(H))
  
  if (nrow(dfF) == 0) {
    return(tibble(
      plane = plane_name,
      n_ok = sum(df$ok, na.rm = TRUE),
      n_feasible = 0,
      H_min = NA_real_, H_q05 = NA_real_, H_med = NA_real_, H_q95 = NA_real_, H_max = NA_real_,
      share_absH_lt_1e06 = NA_real_,
      share_absH_lt_1e04 = NA_real_,
      share_absH_lt_1e03 = NA_real_,
      share_absH_lt_1e02 = NA_real_
    ))
  }
  
  tibble(
    plane = plane_name,
    n_ok = sum(df$ok, na.rm = TRUE),
    n_feasible = nrow(dfF),
    H_min = min(dfF$H),
    H_q05 = as.numeric(quantile(dfF$H, 0.05)),
    H_med = median(dfF$H),
    H_q95 = as.numeric(quantile(dfF$H, 0.95)),
    H_max = max(dfF$H),
    share_absH_lt_1e06 = mean(abs(dfF$H) < 1e-6),
    share_absH_lt_1e04 = mean(abs(dfF$H) < 1e-4),
    share_absH_lt_1e03 = mean(abs(dfF$H) < 1e-3),
    share_absH_lt_1e02 = mean(abs(dfF$H) < 1e-2)
  )
}

# ---- detect plane csvs (excluding scores) ----
csv_files <- fs::dir_ls(OUT_DIR, regexp = "\\.csv$")
bn <- basename(csv_files)

plane_files <- csv_files[grepl("^hopf_plane_", bn) & !grepl("scores", bn)]
print(plane_files)

# ---- run audits ----
plane_audits <- lapply(plane_files, function(f) {
  df <- read_csv_safe(f)
  nm <- fs::path_ext_remove(fs::path_file(f))
  audit_plane_df(df, plane_name = nm)
}) %>% bind_rows()

print(plane_audits)
write_csv(plane_audits, fs::path(OUT_DIR, "TABLE_plane_H_audit.csv"))

# ---- plot: H median + [5%,95%] within feasible region ----
if (nrow(plane_audits) > 0) {
  pH <- plane_audits %>%
    mutate(plane = factor(plane, levels = plane)) %>%
    ggplot(aes(x = plane, y = H_med)) +
    geom_point() +
    geom_errorbar(aes(ymin = H_q05, ymax = H_q95), width = 0.15) +
    coord_flip() +
    labs(
      title = "Plane audit: H median with 5–95% interval (feasible RH only)",
      x = NULL, y = "H"
    )
  
  ggsave(fs::path(OUT_DIR, "PLOT_plane_H_quantiles.png"),
         pH, width = 8.5, height = 5.2, dpi = 160)
  print(pH)
}


head(plane_audits)
tail(plane_audits)


