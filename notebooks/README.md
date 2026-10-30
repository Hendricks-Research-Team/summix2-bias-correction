# Manuscript Figures

This directory contains the manuscript plotting notebook for the Summix2 small-reference-sample analysis.

## Main notebook

`manuscript_plots.Rmd` is the single source for publication-ready manuscript figures and supplementary figures.

The notebook reads finalized analysis outputs produced elsewhere in the repository. It should not rerun simulations, estimate power-law parameters, or apply bias corrections.

## Figures

The manuscript plotting workflow includes:

- Figure 1 — Scenario 1: total absolute bias across reference sample sizes.
- Figure 2 — Scenarios 2 and 3: absolute bias in the varying reference group, absolute bias in the other reference groups, and total absolute bias.
- Figure 3 — correction-method evaluation using the Scenario 1–3 simulation results.
- Figure 4 — real-data stress-test and correction evaluation.
- Figure S1 — Scenario 1 reference-group-specific absolute bias.
- Figure S2 — Scenario 2 reference-group-specific absolute bias.
- Figure S3 — Scenario 3 reference-group-specific absolute bias.

## Terminology

Use the following terminology consistently:

- `reference group`
- `reference sample size (N)`
- `varying reference group`
- `fixed reference groups`
- `absolute bias` for simulation reference-group-specific absolute deviation from known truth
- `total absolute bias` for simulation `total_absolute_deviation`
- `absolute deviation` for real-data deviation from the benchmark
- `total absolute deviation` for the corresponding real-data summary

Reference groups are labeled:

- AFR
- EAS
- EUR
- IAM
- SAS

Simulated admixed populations are labeled:

- Latinx-like
- African American-like
- Peruvian-like

## Inputs

The notebook reads outputs from:

```text
results/simulation/scenario_1_equal_sizes/
results/simulation/scenario_2_3_one_varying/
results/correction/
results/real_data/
```

The exact input files used by each figure should be defined in the corresponding notebook section.

## Outputs

Publication figures should be saved under:

```text
results/figures/
```

PDF should be used for manuscript-quality vector output, with PNG versions generated when useful for previewing or sharing.

## Reproducibility rule

Keep analysis and plotting separate.

The notebook may:

- read finalized `.rds` outputs;
- reshape data for plotting;
- validate expected columns and replicate counts;
- create figures;
- save figures.

The notebook should not:

- rerun Summix2 simulations;
- regenerate Dirichlet parameter-estimation simulations;
- fit final power-law parameters;
- apply the `1/N` or `a/N^b` corrections;
- overwrite analysis results.

This separation keeps simulation, parameter estimation, correction application, correction evaluation, and manuscript visualization reproducible and independently auditable.
