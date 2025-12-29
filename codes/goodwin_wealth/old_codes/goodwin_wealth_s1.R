suppressPackageStartupMessages({
  library(fs)
  library(jsonlite)
})

# ============================================================
# PROCESS 0: Ensure register_run() exists (load or create module)
# ============================================================

# 1) Try to locate an existing module in common places
candidates <- c(
  "registry_module.R",
  "registry/register_run.R",
  "registry/registry_module.R",
  "goodwin_wealth/registry_module.R",
  "wealth_goodwin/registry_module.R"
)

found <- candidates[file.exists(candidates)]
if (length(found) > 0) {
  cat("Found registry module at:", found[1], "\n")
  source(found[1])
} else {
  
  cat("No registry module found. Creating minimal registry_module.R ...\n")
  
  # 2) Write a minimal module that implements register_run()
  module_code <- '
register_run <- function(root = "registry",
                         RunID, Case, mu_name, mu_value, model_stamp,
                         params = list(), reduced = list(), audit = list(), tail = list(),
                         stage1_packet_text = "", analysis_text_md = "",
                         plot_paths = c()) {

  dir.create(root, showWarnings = FALSE, recursive = TRUE)
  runs_dir <- file.path(root, "runs")
  dir.create(runs_dir, showWarnings = FALSE, recursive = TRUE)

  run_dir <- file.path(runs_dir, RunID)
  plots_dir <- file.path(run_dir, "plots")
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

  # save packet + analysis
  packet_path <- file.path(run_dir, "run_packet.txt")
  analysis_path <- file.path(run_dir, "analysis_I_VI.md")
  writeLines(stage1_packet_text, packet_path)
  writeLines(analysis_text_md, analysis_path)

  # copy plots
  saved_plots <- list()
  if (length(plot_paths) > 0) {
    for (nm in names(plot_paths)) {
      src <- plot_paths[[nm]]
      if (is.na(src) || !file.exists(src)) next
      dst <- file.path(plots_dir, paste0(nm, ".png"))
      file.copy(src, dst, overwrite = TRUE)
      saved_plots[[nm]] <- normalizePath(dst, winslash = "/", mustWork = FALSE)
    }
  }

  # JSON artifacts
  json_path <- file.path(run_dir, "artifacts.json")
  artifacts <- list(
    RunID = RunID, Date = as.character(Sys.Date()),
    Case = Case, mu_name = mu_name, mu_value = mu_value,
    model_stamp = model_stamp,
    params = params, reduced = reduced, audit = audit, tail = tail,
    packet_path = normalizePath(packet_path, winslash = "/", mustWork = FALSE),
    analysis_path = normalizePath(analysis_path, winslash = "/", mustWork = FALSE),
    plots_dir = normalizePath(plots_dir, winslash = "/", mustWork = FALSE)
  )
  jsonlite::write_json(artifacts, json_path, pretty = TRUE, auto_unbox = TRUE)

  # manifest upsert
  manifest_path <- file.path(root, "run_manifest.csv")

  row <- data.frame(
    RunID = RunID,
    Date = as.character(Sys.Date()),
    Case = Case,
    mu_name = mu_name,
    mu_value = mu_value,
    model_stamp = model_stamp,

    alpha = params$alpha %||% NA, beta = params$beta %||% NA, delta = params$delta %||% NA,
    phi0 = params$phi0 %||% NA, phi1 = params$phi1 %||% NA,
    i = params$i %||% NA, sigma = params$sigma %||% NA, u_star = params$u_star %||% NA, lambda_u = params$lambda_u %||% NA,

    e_star = reduced$e_star %||% NA, omega_star = reduced$omega_star %||% NA, d_star = reduced$d_star %||% NA,
    H = reduced$H %||% NA, dH_dmu = reduced$dH_dmu %||% NA,

    min_e = audit$min_e %||% NA, min_omega = audit$min_omega %||% NA, min_d = audit$min_d %||% NA,
    min_u = audit$min_u %||% NA, min_xi = audit$min_xi %||% NA,
    min_Delta = audit$min_Delta %||% NA, pct_Delta_le0 = audit$pct_Delta_le0 %||% NA,
    pct_kappa_clamped = audit$pct_kappa_clamped %||% NA,

    tail_mean_e = tail$tail_mean_e %||% NA, tail_sd_e = tail$tail_sd_e %||% NA,
    tail_mean_omega = tail$tail_mean_omega %||% NA, tail_sd_omega = tail$tail_sd_omega %||% NA,
    tail_mean_d = tail$tail_mean_d %||% NA, tail_sd_d = tail$tail_sd_d %||% NA,
    tail_mean_u = tail$tail_mean_u %||% NA, tail_sd_u = tail$tail_sd_u %||% NA,
    tail_mean_xi = tail$tail_mean_xi %||% NA, tail_sd_xi = tail$tail_sd_xi %||% NA,
    period_est = tail$period_est %||% NA,

    plots_dir = normalizePath(plots_dir, winslash = "/", mustWork = FALSE),
    packet_path = normalizePath(packet_path, winslash = "/", mustWork = FALSE),
    analysis_path = normalizePath(analysis_path, winslash = "/", mustWork = FALSE),
    json_path = normalizePath(json_path, winslash = "/", mustWork = FALSE),

    stringsAsFactors = FALSE
  )

  if (file.exists(manifest_path)) {
    man <- read.csv(manifest_path, stringsAsFactors = FALSE)
    man <- man[man$RunID != RunID, , drop = FALSE]
    man <- rbind(man, row)
  } else {
    man <- row
  }
  write.csv(man, manifest_path, row.names = FALSE)

  list(run_dir = run_dir, plots_dir = plots_dir, saved_plots = saved_plots,
       packet_path = packet_path, analysis_path = analysis_path, json_path = json_path)
}

`%||%` <- function(a, b) if (is.null(a) || !is.finite(a)) b else a
'
writeLines(module_code, "registry_module.R")
source("registry_module.R")
}

# 3) Verify
stopifnot(exists("register_run"), is.function(register_run))
cat("OK: register_run() is loaded.\n")



suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
  library(fs)
  library(tibble)
  library(jsonlite)
})

# ============================================================
# PROCESS 1: archive ONE run (df) with plots + stats + registry
# ============================================================

# ---------- 0) USER INPUTS (edit these) ----------
RunID <- "R12_lambda_aboveHopf"
Case  <- "aboveHopf"        # belowHopf / aboveHopf / baseline / scan
mu_name  <- "lambda"
mu_value <- "60.0"          # or "0.85" if that's your label; keep as string ok
model_stamp <- "wealth_goodwin_hopf_cycles.R"

# Choose which simulation to archive:
# dat <- res$below
dat <- res$above

out_dir <- "wealth_goodwin"
registry_root <- "registry"

tail_frac <- 0.25
t_end_label <- max(dat$time, na.rm = TRUE)

# ---------- 1) small helpers ----------
read_if_exists <- function(path) {
  if (file.exists(path)) paste(readLines(path, warn = FALSE), collapse = "\n") else ""
}

tail_window <- function(df, frac = 0.25) {
  n <- nrow(df)
  i0 <- max(1, floor((1 - frac) * n))
  df[i0:n, , drop = FALSE]
}

summ_mean_sd <- function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))

# crude period proxy: uses peaks in e(t) in tail
period_proxy <- function(df_tail) {
  x <- df_tail$e
  t <- df_tail$time
  if (length(x) < 50) return(NA_real_)
  # local maxima
  idx <- which(diff(sign(diff(x))) == -2) + 1
  idx <- idx[idx > 2 & idx < (length(x) - 1)]
  if (length(idx) < 3) return(NA_real_)
  mean(diff(t[idx]), na.rm = TRUE)
}

# ---------- 2) stats (tail + admissibility) ----------
stopifnot(is.data.frame(dat) || inherits(dat, "tbl_df"))
stopifnot(all(c("time","e","omega","d","u","xi") %in% names(dat)))

dat_tail <- tail_window(dat, tail_frac)

tail_stats <- list(
  tail_mean_e     = summ_mean_sd(dat_tail$e)[["mean"]],
  tail_sd_e       = summ_mean_sd(dat_tail$e)[["sd"]],
  tail_mean_omega = summ_mean_sd(dat_tail$omega)[["mean"]],
  tail_sd_omega   = summ_mean_sd(dat_tail$omega)[["sd"]],
  tail_mean_d     = summ_mean_sd(dat_tail$d)[["mean"]],
  tail_sd_d       = summ_mean_sd(dat_tail$d)[["sd"]],
  tail_mean_u     = summ_mean_sd(dat_tail$u)[["mean"]],
  tail_sd_u       = summ_mean_sd(dat_tail$u)[["sd"]],
  tail_mean_xi    = summ_mean_sd(dat_tail$xi)[["mean"]],
  tail_sd_xi      = summ_mean_sd(dat_tail$xi)[["sd"]],
  period_est      = period_proxy(dat_tail)
)

audit_stats <- list(
  min_e = min(dat$e, na.rm = TRUE),
  min_omega = min(dat$omega, na.rm = TRUE),
  min_d = min(dat$d, na.rm = TRUE),
  min_u = min(dat$u, na.rm = TRUE),
  min_xi = min(dat$xi, na.rm = TRUE),
  min_Delta = if ("Delta" %in% names(dat)) min(dat$Delta, na.rm = TRUE) else NA_real_,
  pct_Delta_le0 = if ("Delta" %in% names(dat)) mean(dat$Delta <= 0, na.rm = TRUE) else NA_real_,
  pct_kappa_clamped = if ("kappa_clamped" %in% names(dat)) mean(dat$kappa_clamped == 1, na.rm = TRUE) else NA_real_
)

# ---------- 3) plot pack (timeseries_end, tail_zoom, phase_end with direction) ----------
dir_create(out_dir)

add_tail_flag <- function(df, frac = 0.25) {
  n <- nrow(df)
  i0 <- max(1, floor((1 - frac) * n))
  df %>%
    mutate(idx = row_number(),
           is_tail = idx >= i0,
           t_tail_start = time[i0])
}

plot_timeseries_end <- function(df, vars = c("e","omega","d","u","xi"),
                                tag, out_dir, tail_frac = 0.25) {
  
  df2 <- add_tail_flag(df, tail_frac)
  
  df_long <- df2 %>%
    select(time, all_of(vars), is_tail, t_tail_start) %>%
    pivot_longer(cols = all_of(vars), names_to = "var", values_to = "value")
  
  df_end <- df_long %>% group_by(var) %>% slice_tail(n = 1) %>% ungroup()
  
  p <- ggplot(df_long, aes(x = time, y = value)) +
    geom_line(aes(alpha = is_tail)) +
    geom_vline(xintercept = unique(df_long$t_tail_start), linetype = "dashed") +
    geom_point(data = df_end, size = 2) +
    geom_text(
      data = df_end,
      aes(label = sprintf("end=%.4g", value)),
      hjust = -0.05, vjust = 0.5,
      check_overlap = TRUE
    ) +
    facet_wrap(~var, scales = "free_y", ncol = 1) +
    scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1.0), guide = "none") +
    labs(
      title = paste0("Time paths (tail emphasized) - ", tag),
      subtitle = paste0("Tail = last ", round(100*tail_frac), "%; dashed line = tail start"),
      x = "t", y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(5.5, 35, 5.5, 5.5))
  
  path <- file.path(out_dir, paste0("timeseries_end_", tag, ".png"))
  ggsave(path, p, width = 7.5, height = 9, dpi = 160)
  invisible(path)
}

plot_tail_zoom <- function(df, vars = c("e","omega","d","u","xi"),
                           tag, out_dir, tail_frac = 0.25) {
  
  df_tail <- tail_window(df, tail_frac) %>%
    select(time, all_of(vars)) %>%
    pivot_longer(-time, names_to = "var", values_to = "value")
  
  p <- ggplot(df_tail, aes(x = time, y = value)) +
    geom_line() +
    facet_wrap(~var, scales = "free_y", ncol = 1) +
    labs(
      title = paste0("Tail zoom (final dynamics) - ", tag),
      subtitle = paste0("Only last ", round(100*tail_frac), "% of time"),
      x = "t", y = NULL
    )
  
  path <- file.path(out_dir, paste0("tail_zoom_", tag, ".png"))
  ggsave(path, p, width = 7, height = 9, dpi = 160)
  invisible(path)
}

phase_end_plot <- function(df, x, y, tag, out_dir,
                           tail_frac = 0.25, arrow_every = 40) {
  
  df <- df %>% arrange(time)
  
  n <- nrow(df)
  i0 <- max(1, floor((1 - tail_frac) * n))
  df <- df %>% mutate(idx = row_number(), is_tail = idx >= i0)
  
  p_start <- df[1, ]
  p_end   <- df[n, ]
  
  tail_df <- df %>% filter(is_tail) %>%
    mutate(x0 = .data[[x]],
           y0 = .data[[y]],
           x1 = lead(.data[[x]]),
           y1 = lead(.data[[y]])) %>%
    filter(!is.na(x1), idx %% arrow_every == 0)
  
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_path(aes(alpha = is_tail)) +
    scale_alpha_manual(values = c(`FALSE` = 0.25, `TRUE` = 1.0), guide = "none") +
    geom_point(data = p_start, size = 2) +
    geom_point(data = p_end, size = 2) +
    geom_segment(
      data = tail_df,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      arrow = arrow(length = unit(0.12, "inches")),
      inherit.aes = FALSE
    ) +
    labs(
      title = paste0("Phase (tail emphasized + direction): ", y, " vs ", x, " - ", tag),
      subtitle = paste0("Tail = last ", round(100*tail_frac), "% | points mark start/end"),
      x = x, y = y
    )
  
  path <- file.path(out_dir, paste0("phase_end_", y, "_vs_", x, "_", tag, ".png"))
  ggsave(path, p, width = 6.5, height = 5.5, dpi = 160)
  invisible(path)
}

tag <- paste0(mu_name, "_", Case)

plot_files <- c(
  timeseries_end = plot_timeseries_end(dat, tag = tag, out_dir = out_dir, tail_frac = tail_frac),
  tail_zoom      = plot_tail_zoom(dat, tag = tag, out_dir = out_dir, tail_frac = tail_frac),
  phase_omega_e  = phase_end_plot(dat, x="e",     y="omega", tag = tag, out_dir = out_dir, tail_frac = tail_frac),
  phase_d_omega  = phase_end_plot(dat, x="omega", y="d",     tag = tag, out_dir = out_dir, tail_frac = tail_frac),
  phase_u_d      = phase_end_plot(dat, x="d",     y="u",     tag = tag, out_dir = out_dir, tail_frac = tail_frac),
  phase_xi_d     = phase_end_plot(dat, x="d",     y="xi",    tag = tag, out_dir = out_dir, tail_frac = tail_frac),
  phase_xi_u     = phase_end_plot(dat, x="u",     y="xi",    tag = tag, out_dir = out_dir, tail_frac = tail_frac)
)

stopifnot(all(file.exists(plot_files)))

# ---------- 4) write packet + analysis placeholders ----------
stage1_packet_text <- paste0(
  "========================\nSTAGE 1 RUN PACKET (v1)\n========================\n\n",
  "(1) HEADER\n",
  "RunID: ", RunID, "\n",
  "Date: ", Sys.Date(), "\n",
  "Model stamp: ", model_stamp, "\n",
  "Case: ", Case, "\n",
  "μ varied: ", mu_name, "\n",
  "μ value(s): ", mu_value, "\n\n",
  "(5) SIMULATION: TAIL SUMMARY (last ", round(100*tail_frac), "%)\n",
  sprintf("e mean/sd: %.6g / %.6g\n", tail_stats$tail_mean_e, tail_stats$tail_sd_e),
  sprintf("omega mean/sd: %.6g / %.6g\n", tail_stats$tail_mean_omega, tail_stats$tail_sd_omega),
  sprintf("d mean/sd: %.6g / %.6g\n", tail_stats$tail_mean_d, tail_stats$tail_sd_d),
  sprintf("u mean/sd: %.6g / %.6g\n", tail_stats$tail_mean_u, tail_stats$tail_sd_u),
  sprintf("xi mean/sd: %.6g / %.6g\n", tail_stats$tail_mean_xi, tail_stats$tail_sd_xi),
  sprintf("Period est (proxy): %.6g\n\n", tail_stats$period_est),
  "(6) ADMISSIBILITY\n",
  sprintf("min(e)=%.6g min(omega)=%.6g min(d)=%.6g min(u)=%.6g min(xi)=%.6g\n",
          audit_stats$min_e, audit_stats$min_omega, audit_stats$min_d, audit_stats$min_u, audit_stats$min_xi),
  sprintf("min(Delta)=%.6g | pct Delta<=0 = %.6g\n",
          audit_stats$min_Delta, audit_stats$pct_Delta_le0),
  sprintf("pct kappa clamped = %.6g\n\n", audit_stats$pct_kappa_clamped),
  "(7) PLOTS PROVIDED\n",
  paste(names(plot_files), collapse = ", "), "\n"
)

analysis_text_md <- read_if_exists("analysis_I_VI.md")
if (analysis_text_md == "") {
  analysis_text_md <- "I. Run summary\n\nII. Attractor\n\nIII. Regime classification\n\nIV. Tail diagnostics table\n\nV. Mechanism narrative\n\nVI. Next-step checklist\n"
}

writeLines(stage1_packet_text, "run_packet_filled.txt")
writeLines(analysis_text_md, "analysis_I_VI.md")

# ---------- 5) register_run (requires your previously defined registry module) ----------
if (!exists("register_run")) stop("register_run() not found. Source your registry module first.")

out <- register_run(
  root = registry_root,
  RunID = RunID,
  Case  = Case,
  mu_name  = mu_name,
  mu_value = mu_value,
  model_stamp = model_stamp,
  
  params = list(alpha = if (exists("alpha")) alpha else NA,
                beta  = if (exists("beta")) beta else NA,
                delta = if (exists("delta")) delta else NA,
                phi0  = if (exists("phi0")) phi0 else NA,
                phi1  = if (exists("phi1")) phi1 else NA,
                i     = if (exists("i")) i else NA,
                sigma = if (exists("sigma")) sigma else NA,
                u_star = if (exists("u_star")) u_star else NA,
                lambda_u = if (exists("lambda_u")) lambda_u else NA),
  
  reduced = list(e_star=NA, omega_star=NA, d_star=NA, H=NA, dH_dmu=NA),
  audit = audit_stats,
  tail  = tail_stats,
  
  stage1_packet_text = stage1_packet_text,
  analysis_text_md   = analysis_text_md,
  plot_paths = plot_files
)

# ---------- 6) print what you should paste back to me ----------
manifest <- read.csv(file.path(registry_root, "run_manifest.csv"), stringsAsFactors = FALSE)
manifest_row <- manifest[manifest$RunID == RunID, ]
print(manifest_row)

cat("\nSaved run folder:\n", out$run_dir, "\n")
cat("\nSaved plots:\n"); print(out$saved_plots)
cat("\nDelta audit: min_Delta=", audit_stats$min_Delta, " | pct_Delta_le0=", audit_stats$pct_Delta_le0, "\n")
