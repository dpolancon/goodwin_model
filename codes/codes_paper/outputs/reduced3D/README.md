# Reduced 3D Goodwin–Minsky (quasi-steady financialization) outputs

## What this run does
- Computes admissible interior steady state (closed-form) + saves steady_state.csv
- Builds Jacobian at SS, RH/Hopf functional H, eigenvalues
- Scans rF_bar for Hopf crossings (RH-positivity gated) and refines roots
- Runs simulations for regimes around Hopf roots + baseline
- Uses a common initial condition across regimes if common_state0=TRUE
- Classifies tail behavior (fixed point / damped / limit cycle / runaway/clamp)

## Key files
- steady_state.csv, jacobian.csv, rh_hopf.csv, eigenvalues.csv
- hopf_scan.csv, hopf_roots.csv, hopf_H_curve.png
- regime_summary.csv
- sim_<tag>.csv and plots: states_*, finance_*, phase_*, traj3d_*
- sessionInfo.txt

## Notes
- common_state0 = TRUE (reference tag: between_roots)
