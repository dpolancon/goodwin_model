# ============================================================
# 10) 3D Hopf visualizations (interactive HTML via plotly)
# ============================================================
if (!requireNamespace("plotly", quietly = TRUE) ||
    !requireNamespace("htmlwidgets", quietly = TRUE)) {
  message("Skipping 3D plots: install.packages(c('plotly','htmlwidgets'))")
} else {
  
  hopf_roots <- pp$hopf_roots
  hopf_df <- hopf_roots %>% filter(!is.na(hopf_rF))
  
  if (nrow(hopf_df) == 0) {
    message("No Hopf roots found (all NA). Skipping 3D Hopf plots.")
  } else {
    
    # Helper: make z-matrix for plotly surface from a regular grid
    make_surface_matrix <- function(df, xvar, yvar, zvar) {
      xs <- sort(unique(df[[xvar]]))
      ys <- sort(unique(df[[yvar]]))
      full <- tidyr::expand_grid(x = xs, y = ys) %>%
        left_join(df %>% transmute(x = .data[[xvar]], y = .data[[yvar]], z = .data[[zvar]]),
                  by = c("x","y"))
      zmat <- matrix(full$z, nrow = length(xs), ncol = length(ys), byrow = FALSE)
      list(x = xs, y = ys, z = zmat)
    }
    
    # --- Choose default slices (median of grid) ---
    sigma_slice <- if (is.na(GRID_EXT$sigma_slice)) median(unique(grid$sigma)) else GRID_EXT$sigma_slice
    psi_slice   <- if (is.na(GRID_EXT$psi_slice))   median(unique(grid$psi))   else GRID_EXT$psi_slice
    
    # ------------------------------------------------------------
    # PLOT 1) 3D parameter-space map: (sigma, psi, phi2), color = hopf_rF
    # Shows where a Hopf exists in parameter space, and the required rF level.
    # ------------------------------------------------------------
    p1 <- plotly::plot_ly(
      hopf_df,
      x = ~sigma, y = ~psi, z = ~phi2,
      color = ~hopf_rF,
      type = "scatter3d", mode = "markers",
      marker = list(size = 2, opacity = 0.75)
    ) %>%
      plotly::layout(
        title = "Hopf existence in parameter space (color = hopf_rF)",
        scene = list(
          xaxis = list(title = "sigma"),
          yaxis = list(title = "psi"),
          zaxis = list(title = "phi2")
        )
      )
    
    htmlwidgets::saveWidget(p1, file.path(plot_dir, "hopf_paramspace_3d.html"), selfcontained = TRUE)
    
    # ------------------------------------------------------------
    # PLOT 2) 3D Hopf surface view: (psi, phi2, hopf_rF), color = sigma
    # Shows how the Hopf boundary in rF moves with (psi,phi2) and shifts with sigma.
    # ------------------------------------------------------------
    p2 <- plotly::plot_ly(
      hopf_df,
      x = ~psi, y = ~phi2, z = ~hopf_rF,
      color = ~sigma,
      type = "scatter3d", mode = "markers",
      marker = list(size = 2, opacity = 0.75)
    ) %>%
      plotly::layout(
        title = "Hopf boundary surface cloud: z = hopf_rF (color = sigma)",
        scene = list(
          xaxis = list(title = "psi"),
          yaxis = list(title = "phi2"),
          zaxis = list(title = "hopf_rF")
        )
      )
    
    htmlwidgets::saveWidget(p2, file.path(plot_dir, "hopf_surface_cloud_3d.html"), selfcontained = TRUE)
    
    # ------------------------------------------------------------
    # PLOT 3) True 3D surface (slice): fix sigma = sigma_slice,
    # surface: z = hopf_rF over (psi x phi2)
    # ------------------------------------------------------------
    hopf_slice_sigma <- hopf_df %>% filter(abs(sigma - sigma_slice) < 1e-12)
    
    if (nrow(hopf_slice_sigma) >= 6) { # minimal sanity
      surf <- make_surface_matrix(hopf_slice_sigma, xvar = "psi", yvar = "phi2", zvar = "hopf_rF")
      
      p3 <- plotly::plot_ly(
        x = surf$x, y = surf$y, z = t(surf$z), # transpose for plotly orientation
        type = "surface"
      ) %>%
        plotly::layout(
          title = paste0("Hopf surface slice (sigma = ", sigma_slice, "): z = hopf_rF(psi,phi2)"),
          scene = list(
            xaxis = list(title = "psi"),
            yaxis = list(title = "phi2"),
            zaxis = list(title = "hopf_rF")
          )
        )
      
      htmlwidgets::saveWidget(p3, file.path(plot_dir, "hopf_surface_slice_sigma.html"), selfcontained = TRUE)
    } else {
      message("Not enough Hopf points for sigma slice surface. Skipping hopf_surface_slice_sigma.html")
    }
    
    # Bonus: surface slice fixing psi (often also useful)
    hopf_slice_psi <- hopf_df %>% filter(abs(psi - psi_slice) < 1e-12)
    if (nrow(hopf_slice_psi) >= 6) {
      surf2 <- make_surface_matrix(hopf_slice_psi, xvar = "sigma", yvar = "phi2", zvar = "hopf_rF")
      
      p4 <- plotly::plot_ly(
        x = surf2$x, y = surf2$y, z = t(surf2$z),
        type = "surface"
      ) %>%
        plotly::layout(
          title = paste0("Hopf surface slice (psi = ", psi_slice, "): z = hopf_rF(sigma,phi2)"),
          scene = list(
            xaxis = list(title = "sigma"),
            yaxis = list(title = "phi2"),
            zaxis = list(title = "hopf_rF")
          )
        )
      
      htmlwidgets::saveWidget(p4, file.path(plot_dir, "hopf_surface_slice_psi.html"), selfcontained = TRUE)
    }
  }
}
