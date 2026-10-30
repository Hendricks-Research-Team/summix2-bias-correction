# Scenario 1: Equal Reference Sample Sizes

Scenario 1 evaluates Summix2 estimation when all five reference groups have the same reference sample size.

## Reference groups

The five reference groups are:

- AFR
- EAS
- EUR
- IAM
- SAS

## Simulated admixed populations

Scenario 1 is evaluated for:

- Latinx-like
- African American-like
- Peruvian-like

## Design

For each simulation setting, all five reference groups use the same reference sample size:

\[
N_{\mathrm{AFR}}
=
N_{\mathrm{EAS}}
=
N_{\mathrm{EUR}}
=
N_{\mathrm{IAM}}
=
N_{\mathrm{SAS}}.
\]

The validated simulation evaluates each reference sample size across 1,000 simulation replicates.

## Outcomes

For each replicate, the simulation retains the original Summix2 mixing-proportion estimates and calculates:

- signed bias for each reference group;
- absolute bias for each reference group; and
- total absolute bias across the five reference groups.

Total absolute bias is stored as:

`total_absolute_deviation`

and is calculated as the sum of the five reference-group-specific absolute biases.

## Purpose

Scenario 1 measures the effect of jointly reducing the amount of reference data available for all five reference groups.

It is also used downstream to evaluate the linear `1/N` correction and nonlinear `a/N^b` correction.

For Scenario 1 correction evaluation, all five reference groups are directly corrected before the shared post-correction normalization step.

## Outputs

Simulation outputs should be stored under:

`results/simulation/scenario_1_equal_sizes/`

The original simulation results should not be modified when correction methods are applied.
