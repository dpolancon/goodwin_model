# ============================================================
# run_manager.R — minimal run-based folder wiring
# - Stage 1–4 stay in grid_search/
# - Stage 5 + analysis go into system_analysis/runs/<run_id>/
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(tibble)
})

make_run_id <- function(prefix = "run") {
  paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
}

init_run_dirs <- function(
    project_root = "outputs/wealth_goodwin",
    grid_subdir  = file.path("grid_search"),
    runs_subdir  = file.path("system_analysis", "runs"),
    run_id       = make_run_id()
) {
  grid_base <- file.path(project_root, grid_subdir)
  run_base  <- file.path(project_root, runs_subdir, run_id)
  
  # Stage 5 expects these buckets (as in your patched run_stage5)
  dirs <- list(
    # archive/provenance
    grid_base = grid_base,
    
    # run workspace
    run_base  = run_base,
    inputs    = file.path(run_base, "inputs"),
    reports   = file.path(run_base, "reports"),
    
    # if you want to mirror the stage naming
    stage5A = file.path(run_base, "stage5A_candidates"),
    stage5B = file.path(run_base, "stage5B_sims"),
    stage5C = file.path(run_base, "stage5C_robustness"),
    stage5D = file.path(run_base, "stage5D_ranking"),
    
    # keep a pointer to Stage 4 outputs (read-only)
    stage4 = file.path(grid_base, "stage4_scoring")
  )
  
  # create folders
  for (d in dirs[c("run_base","inputs","reports","stage5A","stage5B","stage5C","stage5D")]) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }
  
  dirs
}

write_run_manifest <- function(run_base) {
  files <- list.files(run_base, recursive = TRUE, full.names = TRUE)
  man <- tibble(path = files) %>%
    transform(
      rel_path = sub(paste0("^", run_base, "/?"), "", path),
      bytes    = suppressWarnings(file.size(path)),
      ext      = tools::file_ext(path)
    ) %>%
    man[order(man$rel_path), ]
  
  readr::write_csv(man, file.path(run_base, "manifest.csv"))
  man
}

# optional: snapshot cfg for perfect reproducibility
snapshot_cfg <- function(cfg, run_base) {
  saveRDS(cfg, file.path(run_base, "cfg_snapshot.rds"))
}


# ---- create a new run workspace ----
dirs <- init_run_dirs(project_root = "outputs/wealth_goodwin")

# optional: keep the exact cfg used
snapshot_cfg(cfg, dirs$run_base)

# ---- run Stage 5 into the run folder ----
final_5D <- run_stage5(cfg, dirs)

# ---- manifest the run outputs ----
man <- write_run_manifest(dirs$run_base)
print(head(man, 30))
