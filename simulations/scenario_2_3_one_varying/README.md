# Scenarios 2 and 3: One Varying Reference Group

Scenarios 2 and 3 evaluate the effect of reducing the reference sample size of one reference group while holding the other four reference groups fixed.

## Reference groups

The five reference groups are:

- AFR
- EAS
- EUR
- IAM
- SAS

## Varying reference group

The varying reference group depends on the simulated admixed population:

| Simulated admixed population | Varying reference group |
| --- | --- |
| Latinx-like | IAM |
| African American-like | AFR |
| Peruvian-like | IAM |

## Scenario 2

In Scenario 2, one reference group varies in sample size while the other four reference groups are fixed at:

\[
N = 500.
\]

This evaluates the effect of a small reference group when the remaining reference groups are comparatively large.

## Scenario 3

In Scenario 3, one reference group varies in sample size while the other four reference groups are fixed at:

\[
N = 150.
\]

This evaluates the same varying-reference-group problem when the fixed reference groups are more moderately sized.

## Replicates

Each population, scenario, and varying reference sample-size setting is evaluated across 1,000 simulation replicates.

## Outcomes

For each replicate, the simulation retains the original Summix2 mixing-proportion estimates and calculates:

- signed bias for each reference group;
- absolute bias for each reference group; and
- total absolute bias across all five reference groups.

Total absolute bias is stored as:

`total_absolute_deviation`

## Manuscript evaluation

The Scenario 2/3 figures distinguish:

- absolute bias in the varying reference group;
- absolute bias in the other reference groups; and
- total absolute bias.

This allows the analysis to assess both direct bias in the reference group with reduced sample size and spillover to the other reference groups.

## Correction evaluation

For Scenarios 2 and 3, the varying reference group is directly corrected.

The remaining reference-group estimates are initially retained, after which the complete five-component mixing-proportion vector undergoes the shared post-correction handling and renormalization.

Both the linear `1/N` correction and nonlinear `a/N^b` correction are evaluated using the existing Scenario 2/3 simulation results; the simulations are not rerun for correction evaluation.

## Outputs

Simulation outputs should be stored under:

`results/simulation/scenario_2_3_one_varying/`

The original simulation results should not be overwritten by corrected results.
