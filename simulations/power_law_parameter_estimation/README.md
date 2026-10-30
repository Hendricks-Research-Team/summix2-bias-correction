# Power-Law Parameter Estimation

This directory contains the workflow used to estimate the reference-group-specific parameters `a` and `b` for the nonlinear power-law bias correction.

This workflow is for parameter estimation and is not an additional simulation scenario.

## Files

- `run_dirichlet_simulations.R` — runs the Dirichlet parameter-estimation simulations.
- `submit_dirichlet_simulations.sh` — submits the simulation jobs on Slurm.
- `summarize_dirichlet_bias.R` — summarizes bias across simulation replicates using `mean(abs(bias))` for power-law fitting.
- `fit_power_law_model.R` — fits `a / N^b` for each Dirichlet mixing-proportion combination and reference group.
- `summarize_parameters.R` — summarizes fitted `a` and `b` values across Dirichlet combinations to produce the final reference-group-specific parameters.

## Workflow

```text
run_dirichlet_simulations.R
        ↓
summarize_dirichlet_bias.R
        ↓
fit_power_law_model.R
        ↓
summarize_parameters.R
        ↓
final power-law parameter table
```

The final parameter table is consumed by the downstream correction workflow. Do not hard-code the final `a` and `b` values in correction scripts.
