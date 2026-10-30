# ============================================================
# Evaluate real-data bias corrections against benchmark
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ============================================================
# Configuration
# ============================================================

real_data_dir <- Sys.getenv(
  "REAL_DATA_DIR",
  "../results/real_data"
)

joint_input_file <- Sys.getenv(
  "JOINT_CORRECTED_FILE",
  file.path(
    real_data_dir,
    "real_data_jointly_limited_corrected.rds"
  )
)

vary_afr_input_file <- Sys.getenv(
  "VARY_AFR_CORRECTED_FILE",
  file.path(
    real_data_dir,
    "real_data_vary_AFR_corrected.rds"
  )
)

joint_output_file <- Sys.getenv(
  "JOINT_EVALUATION_FILE",
  file.path(
    real_data_dir,
    "real_data_jointly_limited_evaluation.rds"
  )
)

vary_afr_output_file <- Sys.getenv(
  "VARY_AFR_EVALUATION_FILE",
  file.path(
    real_data_dir,
    "real_data_vary_AFR_evaluation.rds"
  )
)

combined_output_file <- Sys.getenv(
  "COMBINED_EVALUATION_FILE",
  file.path(
    real_data_dir,
    "real_data_correction_evaluation_long.rds"
  )
)

reference_groups <- c(
  "AFR",
  "EAS",
  "EUR",
  "IAM",
  "SAS"
)

method_levels <- c(
  "Uncorrected",
  "1/N",
  "a/N^b"
)

for (
  input_file in c(
    joint_input_file,
    vary_afr_input_file
  )
) {

  if (!file.exists(input_file)) {
    stop(
      "Required corrected real-data file does not exist: ",
      input_file
    )
  }
}

# ============================================================
# Read corrected real-data results
# ============================================================

joint_results <- readRDS(
  joint_input_file
)

vary_afr_results <- readRDS(
  vary_afr_input_file
)

# ============================================================
# Validate required columns
# ============================================================

required_columns <- c(
  reference_groups,
  paste0(
    "linear_corrected_",
    reference_groups
  ),
  paste0(
    "powerlaw_corrected_",
    reference_groups
  ),
  "N_AFR",
  "N_EAS",
  "N_EUR",
  "N_IAM",
  "N_SAS"
)

validate_real_data_columns <- function(
    data,
    data_name
) {

  missing_columns <- setdiff(
    required_columns,
    names(data)
  )

  if (
    length(missing_columns) > 0L
  ) {
    stop(
      data_name,
      " is missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  invisible(TRUE)
}

validate_real_data_columns(
  joint_results,
  "Jointly limited real-data results"
)

validate_real_data_columns(
  vary_afr_results,
  "AFR-varying real-data results"
)

# ============================================================
# Extract common uncorrected benchmark
#
# AFR/EAS/EUR/SAS = 500
# IAM = 45
# ============================================================

get_benchmark_row <- function(
    data,
    data_name
) {

  benchmark_row <- data %>%
    filter(
      N_AFR == 500,
      N_EAS == 500,
      N_EUR == 500,
      N_IAM == 45,
      N_SAS == 500
    )

  if (
    nrow(benchmark_row) != 1L
  ) {
    stop(
      data_name,
      " must contain exactly one ",
      "500/500/500/45/500 benchmark row."
    )
  }

  benchmark_row
}

joint_benchmark_row <- get_benchmark_row(
  joint_results,
  "Jointly limited real-data results"
)

vary_afr_benchmark_row <- get_benchmark_row(
  vary_afr_results,
  "AFR-varying real-data results"
)

joint_benchmark <- as.numeric(
  joint_benchmark_row[
    ,
    reference_groups,
    drop = TRUE
  ]
)

names(
  joint_benchmark
) <- reference_groups

vary_afr_benchmark <- as.numeric(
  vary_afr_benchmark_row[
    ,
    reference_groups,
    drop = TRUE
  ]
)

names(
  vary_afr_benchmark
) <- reference_groups

benchmark_difference <- max(
  abs(
    joint_benchmark -
      vary_afr_benchmark
  )
)

if (
  benchmark_difference >
    1e-10
) {
  stop(
    "The two real-data stress tests do not share ",
    "the same uncorrected benchmark estimates."
  )
}

benchmark <- joint_benchmark

# ============================================================
# Convert one stress test to long method format
# ============================================================

prepare_method_estimates <- function(
    data,
    stress_test_label,
    x_variable
) {

  if (
    !x_variable %in%
      names(data)
  ) {
    stop(
      "Missing x-axis variable: ",
      x_variable
    )
  }

  uncorrected_long <- data %>%
    select(
      all_of(
        x_variable
      ),
      N_AFR,
      N_EAS,
      N_EUR,
      N_IAM,
      N_SAS,
      all_of(
        reference_groups
      )
    ) %>%
    pivot_longer(
      cols = all_of(
        reference_groups
      ),
      names_to =
        "ReferenceGroup",
      values_to =
        "Estimate"
    ) %>%
    mutate(
      Method =
        "Uncorrected"
    )

  linear_long <- data %>%
    select(
      all_of(
        x_variable
      ),
      N_AFR,
      N_EAS,
      N_EUR,
      N_IAM,
      N_SAS,
      all_of(
        paste0(
          "linear_corrected_",
          reference_groups
        )
      )
    ) %>%
    pivot_longer(
      cols = all_of(
        paste0(
          "linear_corrected_",
          reference_groups
        )
      ),
      names_to =
        "ReferenceGroup",
      names_prefix =
        "linear_corrected_",
      values_to =
        "Estimate"
    ) %>%
    mutate(
      Method =
        "1/N"
    )

  powerlaw_long <- data %>%
    select(
      all_of(
        x_variable
      ),
      N_AFR,
      N_EAS,
      N_EUR,
      N_IAM,
      N_SAS,
      all_of(
        paste0(
          "powerlaw_corrected_",
          reference_groups
        )
      )
    ) %>%
    pivot_longer(
      cols = all_of(
        paste0(
          "powerlaw_corrected_",
          reference_groups
        )
      ),
      names_to =
        "ReferenceGroup",
      names_prefix =
        "powerlaw_corrected_",
      values_to =
        "Estimate"
    ) %>%
    mutate(
      Method =
        "a/N^b"
    )

  bind_rows(
    uncorrected_long,
    linear_long,
    powerlaw_long
  ) %>%
    mutate(
      StressTest =
        stress_test_label,
      Method = factor(
        Method,
        levels =
          method_levels
      ),
      ReferenceGroup = factor(
        ReferenceGroup,
        levels =
          reference_groups
      ),
      ReferenceSampleSize =
        .data[
          [x_variable]
        ]
    )
}

# ============================================================
# Prepare jointly limited stress test
# ============================================================

if (
  !"N" %in%
    names(joint_results)
) {
  stop(
    "Jointly limited results are missing N."
  )
}

joint_long <- prepare_method_estimates(
  data =
    joint_results,
  stress_test_label =
    "jointly_limited",
  x_variable =
    "N"
)

# ============================================================
# Prepare AFR-varying stress test
# ============================================================

if (
  !"AFR_N" %in%
    names(vary_afr_results)
) {
  stop(
    "AFR-varying results are missing AFR_N."
  )
}

vary_afr_long <- prepare_method_estimates(
  data =
    vary_afr_results,
  stress_test_label =
    "vary_AFR",
  x_variable =
    "AFR_N"
)

# ============================================================
# Add benchmark and calculate absolute deviation
# ============================================================

benchmark_table <- tibble(
  ReferenceGroup = factor(
    reference_groups,
    levels =
      reference_groups
  ),
  BenchmarkEstimate =
    unname(
      benchmark[
        reference_groups
      ]
    )
)

calculate_deviation <- function(
    long_data
) {

  long_data %>%
    left_join(
      benchmark_table,
      by =
        "ReferenceGroup"
    ) %>%
    mutate(
      Deviation =
        Estimate -
        BenchmarkEstimate,
      AbsoluteDeviation =
        abs(
          Deviation
        )
    )
}

joint_evaluation_long <- calculate_deviation(
  joint_long
)

vary_afr_evaluation_long <- calculate_deviation(
  vary_afr_long
)

# ============================================================
# Calculate total absolute deviation
# ============================================================

summarize_total_deviation <- function(
    long_data
) {

  long_data %>%
    group_by(
      StressTest,
      ReferenceSampleSize,
      Method,
      N_AFR,
      N_EAS,
      N_EUR,
      N_IAM,
      N_SAS
    ) %>%
    summarise(
      TotalAbsoluteDeviation =
        sum(
          AbsoluteDeviation
        ),
      .groups =
        "drop"
    )
}

joint_total <- summarize_total_deviation(
  joint_evaluation_long
)

vary_afr_total <- summarize_total_deviation(
  vary_afr_evaluation_long
)

# ============================================================
# Add total deviation back to reference-group-level results
# ============================================================

join_keys <- c(
  "StressTest",
  "ReferenceSampleSize",
  "Method",
  "N_AFR",
  "N_EAS",
  "N_EUR",
  "N_IAM",
  "N_SAS"
)

joint_evaluation_long <- joint_evaluation_long %>%
  left_join(
    joint_total,
    by =
      join_keys
  )

vary_afr_evaluation_long <- vary_afr_evaluation_long %>%
  left_join(
    vary_afr_total,
    by =
      join_keys
  )

# ============================================================
# Validate benchmark behavior
#
# Uncorrected benchmark must have zero deviation.
# ============================================================

benchmark_check <- bind_rows(
  joint_evaluation_long,
  vary_afr_evaluation_long
) %>%
  filter(
    N_AFR == 500,
    N_EAS == 500,
    N_EUR == 500,
    N_IAM == 45,
    N_SAS == 500,
    Method ==
      "Uncorrected"
  )

if (
  any(
    abs(
      benchmark_check$Deviation
    ) >
      1e-12
  )
) {
  stop(
    "Uncorrected benchmark configuration ",
    "does not have zero deviation from itself."
  )
}

# ============================================================
# Create compact evaluation tables
# ============================================================

joint_evaluation <- joint_evaluation_long %>%
  select(
    StressTest,
    ReferenceSampleSize,
    Method,
    ReferenceGroup,
    Estimate,
    BenchmarkEstimate,
    Deviation,
    AbsoluteDeviation,
    TotalAbsoluteDeviation,
    N_AFR,
    N_EAS,
    N_EUR,
    N_IAM,
    N_SAS
  ) %>%
  arrange(
    ReferenceSampleSize,
    Method,
    ReferenceGroup
  )

vary_afr_evaluation <- vary_afr_evaluation_long %>%
  select(
    StressTest,
    ReferenceSampleSize,
    Method,
    ReferenceGroup,
    Estimate,
    BenchmarkEstimate,
    Deviation,
    AbsoluteDeviation,
    TotalAbsoluteDeviation,
    N_AFR,
    N_EAS,
    N_EUR,
    N_IAM,
    N_SAS
  ) %>%
  arrange(
    ReferenceSampleSize,
    Method,
    ReferenceGroup
  )

combined_evaluation <- bind_rows(
  joint_evaluation,
  vary_afr_evaluation
)

# ============================================================
# Save evaluation results
# ============================================================

saveRDS(
  joint_evaluation,
  joint_output_file,
  compress = TRUE
)

saveRDS(
  vary_afr_evaluation,
  vary_afr_output_file,
  compress = TRUE
)

saveRDS(
  combined_evaluation,
  combined_output_file,
  compress = TRUE
)

message(
  "Saved jointly limited evaluation to: ",
  normalizePath(
    joint_output_file,
    mustWork = FALSE
  ),
  "\nSaved AFR-varying evaluation to: ",
  normalizePath(
    vary_afr_output_file,
    mustWork = FALSE
  ),
  "\nSaved combined long-format evaluation to: ",
  normalizePath(
    combined_output_file,
    mustWork = FALSE
  )
)
