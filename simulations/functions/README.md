# Shared Functions

This directory contains reusable functions shared across the Summix2 small-reference-sample analysis workflows.

## Files

- `prepare_reference_data.R` — prepares the chromosome 22 reference allele-frequency data and constructs the five reference groups: AFR, EAS, EUR, IAM, and SAS.
- `simulation_functions.R` — contains the shared simulation functions used by Scenario 1, Scenarios 2/3, and the power-law parameter-estimation simulations. It generates finite-reference allele frequencies, observed allele frequencies, runs Summix2, and calculates signed bias, absolute bias, and total absolute deviation.
- `correction_functions.R` — contains the shared 1/N and power-law correction functions, parameter validation, direction rule, truncation of negative corrected proportions, and renormalization.

These files contain shared functionality only. Scenario-specific execution scripts remain under `scripts/simulation/`, and correction application/evaluation scripts remain under `scripts/correction/`.
