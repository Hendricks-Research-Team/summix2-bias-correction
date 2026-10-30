# Bias-Correction Evaluation

This directory contains the scripts used to apply and evaluate the two post-estimation bias-correction methods in simulation Scenarios 1–3.

The shared correction functions are stored separately in `scripts/functions/correction_functions.R` because the same functions are also used for the real-data analysis.

## Methods

Two correction approaches are evaluated:

1. `1/N` correction
2. Power-law `a/N^b` correction

The power-law parameters are estimated separately by the workflow in `scripts/simulation/power_law_parameter_estimation/`. The Dirichlet parameter-estimation simulations are not an additional evaluation scenario.

## Workflow

```text
apply_corrections.R
        ↓
evaluate_corrections.R
        ↓
results used for Figure 3 and supplementary figures
```

### `apply_corrections.R`

Reads the validated Scenario 1–3 simulation outputs and applies both correction methods to the original Summix2 estimated mixing proportions.

Correction targeting follows the simulation design:

- Scenario 1: all five reference groups are directly corrected.
- Scenario 2: only the varying reference group is directly corrected.
- Scenario 3: only the varying reference group is directly corrected.

After direct correction, the full five-component vector is processed using the shared post-correction rule in `correction_functions.R`: corrected values are constrained to the valid proportion range, the vector is renormalized to sum to one, and if post-processing produces an all-zero vector, each of the five reference-group proportions is set to `1/5`.

The original uncorrected Summix2 estimates are preserved. The script also recalculates reference-group-specific signed bias, absolute bias, and total absolute bias for each corrected method.

### `evaluate_corrections.R`

Combines the corrected Scenario 1–3 outputs and compares three methods:

- Uncorrected
- `1/N`
- `a/N^b`

Evaluation is based on total absolute bias across the five reference groups. Replicate-level results are retained, and summary statistics are calculated across simulation replicates for each simulation condition and method.

The evaluation output is used by the manuscript plotting workflow, including Figure 3.

## Inputs

The correction workflow requires:

- validated Scenario 1 simulation outputs;
- validated Scenario 2/3 simulation outputs;
- the final power-law parameter summary produced by `summarize_parameters.R`; and
- `scripts/functions/correction_functions.R`.

Do not rerun the original simulations when applying corrections.

## Outputs

Corrected simulation results are written under `results/correction/`, preserving separate Scenario 1 and Scenario 2/3 output directories.

The evaluation step produces the manuscript-ready comparison datasets, including:

```text
correction_method_comparison.rds
correction_method_comparison_summary.rds
```

These outputs compare uncorrected Summix2, `1/N`, and `a/N^b` corrections without modifying the original simulation result files.
