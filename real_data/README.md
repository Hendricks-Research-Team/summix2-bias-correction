# Real-Data Evaluation

This directory contains the real-data stress-test and correction-evaluation workflow used to assess whether the simulation-derived reference-size patterns generalize to realistic reference-data constraints.

Because the true mixing proportions are unknown in the real data, results are evaluated relative to a large-reference benchmark rather than against known truth. Accordingly, this workflow uses the terms **deviation**, **absolute deviation**, and **total absolute deviation**, not bias.

## Reference groups

The five reference groups are:

- AFR
- EAS
- EUR
- IAM
- SAS

## Benchmark

The common uncorrected benchmark is:

- AFR = 500
- EAS = 500
- EUR = 500
- IAM = 45
- SAS = 500

This benchmark represents the largest available reference configuration used in the analysis and is treated as a stability benchmark rather than ground truth.

## Stress tests

### Jointly limited reference groups

The reference sample-size grid is:

```text
5, 15, 25, 35, 45, 75, 100, 150, 250, 500
```

AFR, EAS, EUR, and SAS use the requested value of N.

IAM is limited by available data and therefore uses:

```text
min(N, 45)
```

For correction application in this stress test, all five reference-group estimates are directly corrected before post-correction normalization.

### AFR-varying stress test

AFR varies over the same sample-size grid, while the remaining reference groups are fixed at:

```text
EAS = 500
EUR = 500
IAM = 45
SAS = 500
```

For correction application in this stress test, only AFR is directly corrected. The complete five-component vector is then processed by the shared post-correction normalization rule.

## Workflow

```text
run_real_data_stress_tests.R
        ↓
apply_real_data_corrections.R
        ↓
evaluate_real_data_corrections.R
        ↓
summarize_real_data_results.R
        ↓
Figure 4 in notebooks/manuscript_plots.Rmd
```

## Scripts

### `run_real_data_stress_tests.R`

Runs the uncorrected Summix2 real-data stress tests.

It produces:

```text
real_data_jointly_limited_uncorrected.rds
real_data_vary_AFR_uncorrected.rds
```

The script also checks that the shared large-reference configuration yields consistent benchmark estimates across the two stress tests.

### `apply_real_data_corrections.R`

Applies the two correction methods using the shared functions in:

```text
scripts/functions/correction_functions.R
```

The methods are:

- `1/N`
- `a/N^b`

The power-law correction requires the finalized parameter table produced by the power-law parameter-estimation workflow.

The script preserves the original uncorrected estimates and appends corrected mixing-proportion columns.

### `evaluate_real_data_corrections.R`

Evaluates each method relative to the uncorrected large-reference benchmark.

For reference group \(k\):

\[
D_k =
\left|
\hat{\pi}_{k,\mathrm{configuration}}
-
\hat{\pi}_{k,\mathrm{benchmark}}
\right|.
\]

Total absolute deviation is:

\[
D_{\mathrm{total}}
=
\sum_k D_k.
\]

The script saves reference-group-specific and total deviation results in long format for downstream plotting and summarization.

### `summarize_real_data_results.R`

Creates manuscript-support summaries from the evaluation output.

The summary includes:

- estimated mixing proportions;
- benchmark estimates;
- reference-group-specific absolute deviations;
- total absolute deviation;
- absolute improvement relative to uncorrected Summix2; and
- relative improvement relative to uncorrected Summix2.

Both long- and wide-format summaries are saved.

## Post-correction handling

The correction workflow uses the shared post-correction handling defined in `correction_functions.R`.

After direct correction:

1. corrected values are constrained according to the validated correction implementation;
2. the five mixing proportions are renormalized to sum to one; and
3. if post-processing yields an all-zero vector, each reference group is assigned \(1/K\), where \(K=5\).

## Important

Do not interpret the real-data benchmark as ground truth.

The real-data analysis is intended to evaluate **stability under reference sample-size constraints** and whether correction reduces deviation from the largest available reference configuration.
