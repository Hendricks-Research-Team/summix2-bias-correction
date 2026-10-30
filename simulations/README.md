# Simulation Workflows

This directory contains the simulation workflows used to evaluate the effect of reference sample size on Summix2 mixing-proportion estimates and to estimate the parameters used by the nonlinear bias correction.

The five reference groups are:

- AFR
- EAS
- EUR
- IAM
- SAS

The simulated admixed populations are referred to as:

- Latinx-like
- African American-like
- Peruvian-like

## Simulation design

The manuscript evaluates three simulation scenarios.

### Scenario 1: Equal reference sample sizes

All five reference groups have the same reference sample size.

This scenario evaluates the effect of jointly reducing reference sample size across all reference groups.

See:

`scenario_1_equal_sizes/`

### Scenario 2: One varying reference group, fixed reference groups N = 500

One reference group varies in sample size while the other four reference groups are fixed at N = 500.

The varying reference group is:

- IAM for the Latinx-like population
- AFR for the African American-like population
- IAM for the Peruvian-like population

See:

`scenario_2_3_one_varying/`

### Scenario 3: One varying reference group, fixed reference groups N = 150

One reference group varies in sample size while the other four reference groups are fixed at N = 150.

The varying reference group is the same as in Scenario 2:

- IAM for the Latinx-like population
- AFR for the African American-like population
- IAM for the Peruvian-like population

See:

`scenario_2_3_one_varying/`

## Power-law parameter-estimation simulations

The `power_law_parameter_estimation/` workflow is separate from Scenarios 1–3.

It uses 100 Dirichlet mixing-proportion combinations and reference sample sizes N = 1,...,1000 to estimate the parameters of the nonlinear power-law correction.

These simulations are used for parameter estimation, not as an additional correction-evaluation scenario.

## Simulation outputs

Replicate-level simulation results preserve the original Summix2 estimates and include reference-group-specific signed bias, reference-group-specific absolute bias, and total absolute bias (`total_absolute_deviation`).

Correction methods are applied downstream. The original simulation results should not be overwritten.

## Downstream workflow

```text
Scenario 1–3 simulations
        |
        +--> manuscript plots of uncorrected performance
        |
        +--> correction application
                |
                +--> correction evaluation
                |
                +--> corrected manuscript figures

Power-law parameter-estimation simulations
        |
        +--> summarize mean(abs(bias))
        |
        +--> fit a/N^b
        |
        +--> summarize final a and b
        |
        +--> power-law correction
```

Simulation, parameter estimation, correction application, correction evaluation, and manuscript plotting are intentionally kept as separate stages.
