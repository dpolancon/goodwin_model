############################################################
## Canonical Goodwin (1967) — Benchmark Cycle
############################################################

library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)

## ---------------------------------------------------------
## 1) Targets and exogenous growth
## ---------------------------------------------------------
e_star     <- 0.65
omega_star <- 0.50

alpha <- 0.02
beta  <- 0.01
g     <- alpha + beta

## ---------------------------------------------------------
## 2) Structural calibration
## ---------------------------------------------------------
sigma <- (1 - omega_star) / g

phi1 <- 0.10
phi0 <- alpha - phi1 * e_star
Phi  <- function(e) phi0 + phi1 * e

stopifnot(
  e_star     > 0, e_star     < 1,
  omega_star > 0, omega_star < 1
)

## ---------------------------------------------------------
## 3) Canonical Goodwin ODE
## ---------------------------------------------------------
goodwin <- function(t, state, parms) {
  
  e     <- state["e"]
  omega <- state["omega"]
  
  de     <- ((1 - omega)/sigma - g) * e
  domega <- (Phi(e) - alpha) * omega
  
  list(c(de, domega))
}

## ---------------------------------------------------------
## 4) Simulation design
## ---------------------------------------------------------
amplitudes <- c(0.015, 0.03, 0.05, 0.10)
times <- seq(0, 400, by = 0.1)

## ---------------------------------------------------------
## 5) Run simulations
## ---------------------------------------------------------
sims <- lapply(amplitudes, function(a) {
  
  state0 <- c(
    e     = e_star     + a,
    omega = omega_star - a
  )
  
  sol <- ode(
    y      = state0,
    times  = times,
    func   = goodwin,
    parms  = NULL,
    method = "lsoda",
    rtol   = 1e-8,
    atol   = 1e-10
  )
  
  as.data.frame(sol) |>
    mutate(amplitude = factor(a))
})

df <- bind_rows(sims)

## ---------------------------------------------------------
## 6) Phase diagram
## ---------------------------------------------------------
p_phase <- ggplot(df, aes(x = e, y = omega)) +
  geom_path(aes(group = amplitude),
            linewidth = 0.9,
            colour = "steelblue") +
  geom_point(
    aes(x = e_star, y = omega_star),
    inherit.aes = FALSE,
    shape = 21,
    fill = "black",
    size = 2
  ) +
  labs(
    x = "Employment rate (e)",
    y = "Wage share (ω)"
  ) +
  theme_minimal()

ggsave(
  filename = "fig_goodwin_phase.pdf",
  plot     = p_phase,
  width    = 6,
  height   = 6,
  units    = "in",
  dpi      = 300
)

## ---------------------------------------------------------
## 7) Mild cycle and derived variables
## ---------------------------------------------------------
a_mild <- "0.015"

df_mild <- df |>
  filter(amplitude == a_mild) |>
  mutate(
    a_t  = exp(alpha * time),
    y_pc = a_t * e,
    g_Y  = (1 - omega) / sigma
  )

## ---------------------------------------------------------
## 8) Tidy data for stacked time paths
## ---------------------------------------------------------
df_mild_long <- df_mild |>
  select(time, e, omega, g_Y, y_pc) |>
  pivot_longer(
    cols = c(e, omega, g_Y, y_pc),
    names_to = "variable",
    values_to = "value"
  )

df_mild_long$variable <- factor(
  df_mild_long$variable,
  levels = c("e", "omega", "g_Y", "y_pc"),
  labels = c(
    "Employment rate (e)",
    "Wage share (ω)",
    "Capital accumulation  g_Y = (1 − ω)/σ",
    "Output per capita  y = a · e"
  )
)

## ---------------------------------------------------------
## 9) Stacked time paths (mild cycle)
## ---------------------------------------------------------
p_time_mild <- ggplot(df_mild_long,
                      aes(x = time, y = value)) +
  geom_line(linewidth = 0.5, colour = "darkred") +
  facet_wrap(~ variable, ncol = 1, scales = "free_y") +
  labs(
    x = "Time",
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1.2, "lines")
  )

ggsave(
  filename = "fig_goodwin_timepaths_mild.pdf",
  plot     = p_time_mild,
  width    = 7,
  height   = 8,
  units    = "in",
  dpi      = 300
)

## ---------------------------------------------------------
## 10) Output objects
## ---------------------------------------------------------
# df             : wide data (all amplitudes)
# df_mild        : mild cycle data with derived variables
# df_mild_long   : tidy stacked data for plotting
############################################################
