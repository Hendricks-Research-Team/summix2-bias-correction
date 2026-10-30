# ============================================================
# Summarize real-data correction results
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

evaluation_file <- Sys.getenv(
  "REAL_DATA_EVALUATION_FILE",
  "../results/real_data/real_data_correction_evaluation_long.rds"
)

summary_output_file <- Sys.getenv(
  "REAL_DATA_SUMMARY_FILE",
  "../results/real_data/real_data_correction_summary.rds"
)

wide_output_file <- Sys.getenv(
  "REAL_DATA_SUMMARY_WIDE_FILE",
  "../results/real_data/real_data_correction_summary_wide.rds"
)

reference_groups <- c("AFR", "EAS", "EUR", "IAM", "SAS")
method_levels <- c("Uncorrected", "1/N", "a/N^b")

if (!file.exists(evaluation_file)) {
  stop("Real-data evaluation file does not exist: ", evaluation_file)
}

dir.create(
  dirname(summary_output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

evaluation <- readRDS(evaluation_file)

required_columns <- c(
  "StressTest",
  "ReferenceSampleSize",
  "Method",
  "ReferenceGroup",
  "Estimate",
  "BenchmarkEstimate",
  "Deviation",
  "AbsoluteDeviation",
  "TotalAbsoluteDeviation",
  "N_AFR",
  "N_EAS",
  "N_EUR",
  "N_IAM",
  "N_SAS"
)

missing_columns <- setdiff(required_columns, names(evaluation))

if (length(missing_columns) > 0L) {
  stop(
    "Real-data evaluation is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

evaluation <- evaluation %>%
  mutate(
    Method = factor(Method, levels = method_levels),
    ReferenceGroup = factor(
      ReferenceGroup,
      levels = reference_groups
    )
  )

reference_group_summary <- evaluation %>%
  select(
    StressTest,
    ReferenceSampleSize,
    Method,
    ReferenceGroup,
    Estimate,
    BenchmarkEstimate,
    Deviation,
    AbsoluteDeviation,
    N_AFR,
    N_EAS,
    N_EUR,
    N_IAM,
    N_SAS
  )

total_summary <- evaluation %>%
  distinct(
    StressTest,
    ReferenceSampleSize,
    Method,
    TotalAbsoluteDeviation,
    N_AFR,
    N_EAS,
    N_EUR,
    N_IAM,
    N_SAS
  )

uncorrected_total <- total_summary %>%
  filter(Method == "Uncorrected") %>%
  select(
    StressTest,
    ReferenceSampleSize,
    UncorrectedTotalAbsoluteDeviation =
      TotalAbsoluteDeviation
  )

total_summary <- total_summary %>%
  left_join(
    uncorrected_total,
    by = c(
      "StressTest",
      "ReferenceSampleSize"
    )
  ) %>%
  mutate(
    AbsoluteImprovement =
      UncorrectedTotalAbsoluteDeviation -
      TotalAbsoluteDeviation,
    RelativeImprovement = ifelse(
      UncorrectedTotalAbsoluteDeviation > 0,
      AbsoluteImprovement /
        UncorrectedTotalAbsoluteDeviation,
      NA_real_
    )
  )

summary_long <- reference_group_summary %>%
  left_join(
    total_summary %>%
      select(
        StressTest,
        ReferenceSampleSize,
        Method,
        TotalAbsoluteDeviation,
        UncorrectedTotalAbsoluteDeviation,
        AbsoluteImprovement,
        RelativeImprovement
      ),
    by = c(
      "StressTest",
      "ReferenceSampleSize",
      "Method"
    )
  ) %>%
  arrange(
    StressTest,
    ReferenceSampleSize,
    Method,
    ReferenceGroup
  )

summary_wide <- summary_long %>%
  select(
    StressTest,
    ReferenceSampleSize,
    Method,
    ReferenceGroup,
    Estimate,
    BenchmarkEstimate,
    AbsoluteDeviation,
    TotalAbsoluteDeviation,
    UncorrectedTotalAbsoluteDeviation,
    AbsoluteImprovement,
    RelativeImprovement,
    N_AFR,
    N_EAS,
    N_EUR,
    N_IAM,
    N_SAS
  ) %>%
  pivot_wider(
    names_from = ReferenceGroup,
    values_from = c(
      Estimate,
      BenchmarkEstimate,
      AbsoluteDeviation
    ),
    names_sep = "_"
  ) %>%
  arrange(
    StressTest,
    ReferenceSampleSize,
    Method
  )

benchmark_check <- total_summary %>%
  filter(
    Method == "Uncorrected",
    N_AFR == 500,
    N_EAS == 500,
    N_EUR == 500,
    N_IAM == 45,
    N_SAS == 500
  )

if (nrow(benchmark_check) == 0L) {
  stop(
    "The uncorrected 500/500/500/45/500 benchmark is missing."
  )
}

if (any(
  abs(
    benchmark_check$TotalAbsoluteDeviation
  ) > 1e-12
)) {
  stop(
    "The uncorrected benchmark does not have zero total absolute deviation."
  )
}

saveRDS(
  summary_long,
  summary_output_file,
  compress = TRUE
)

saveRDS(
  summary_wide,
  wide_output_file,
  compress = TRUE
)

message(
  "Saved long-format real-data summary to: ",
  normalizePath(
    summary_output_file,
    mustWork = FALSE
  ),
  "\nSaved wide-format real-data summary to: ",
  normalizePath(
    wide_output_file,
    mustWork = FALSE
  )
)
