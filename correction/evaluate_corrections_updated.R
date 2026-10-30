# ============================================================
# Evaluate bias-correction methods across simulation scenarios
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ============================================================
# Configuration
# ============================================================

corrected_dir <- Sys.getenv(
  "CORRECTED_DIR",
  "../results/correction"
)

output_file <- Sys.getenv(
  "OUTPUT_FILE",
  "../results/correction/correction_method_comparison.rds"
)

summary_file <- Sys.getenv(
  "SUMMARY_FILE",
  "../results/correction/correction_method_comparison_summary.rds"
)

varying_summary_file <- Sys.getenv(
  "VARYING_SUMMARY_FILE",
  "../results/correction/varying_ancestry_correction_summary.rds"
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

if (!dir.exists(corrected_dir)) {
  stop(
    "Corrected simulation directory does not exist: ",
    corrected_dir
  )
}

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(summary_file),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(varying_summary_file),
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Read corrected simulation outputs
# ============================================================

corrected_files <- list.files(
  corrected_dir,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

# Exclude evaluation outputs if the script is rerun.
corrected_files <- corrected_files[
  !basename(corrected_files) %in% c(
    basename(output_file),
    basename(summary_file),
    basename(varying_summary_file)
  )
]

if (length(corrected_files) == 0L) {
  stop(
    "No corrected simulation .rds files were found in: ",
    corrected_dir
  )
}

simulation_results <- bind_rows(
  lapply(
    corrected_files,
    function(input_file) {
      result <- readRDS(input_file)
      result$source_file <- input_file
      result
    }
  )
)

# ============================================================
# Validate required columns
# ============================================================

required_metadata <- c(
  "Scenario",
  "Group",
  "Replicate",
  "N_AFR",
  "N_EAS",
  "N_EUR",
  "N_IAM",
  "N_SAS"
)

required_bias_columns <- c(
  "total_absolute_deviation",
  "linear_total_absolute_deviation",
  "powerlaw_total_absolute_deviation"
)

required_ancestry_bias_columns <- c(
  paste0("absolute_bias_", reference_groups),
  paste0("linear_absolute_bias_", reference_groups),
  paste0("powerlaw_absolute_bias_", reference_groups)
)

required_columns <- c(
  required_metadata,
  required_bias_columns,
  required_ancestry_bias_columns
)

missing_columns <- setdiff(
  required_columns,
  names(simulation_results)
)

if (length(missing_columns) > 0L) {
  stop(
    "Corrected simulation results are missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

if (anyNA(
  simulation_results[
    ,
    required_columns,
    drop = FALSE
  ]
)) {
  stop(
    "Required correction-evaluation columns contain missing values."
  )
}

# ============================================================
# Preserve optional Scenario 2/3 metadata
# ============================================================

optional_metadata <- intersect(
  c(
    "VaryingAncestry",
    "N_var",
    "N_fixed"
  ),
  names(simulation_results)
)

# ============================================================
# Define plotted reference sample size
# ============================================================

if ("N_var" %in% names(simulation_results)) {

  simulation_results <- simulation_results %>%
    mutate(
      reference_sample_size = case_when(
        Scenario == "scenario_1_equal_sizes" ~ N_AFR,
        !is.na(N_var) ~ N_var,
        Group == "AfricanAmerican" ~ N_AFR,
        Group %in% c("Latinx", "Peruvian") ~ N_IAM,
        TRUE ~ NA_real_
      )
    )

} else {

  simulation_results <- simulation_results %>%
    mutate(
      reference_sample_size = case_when(
        Scenario == "scenario_1_equal_sizes" ~ N_AFR,
        Group == "AfricanAmerican" ~ N_AFR,
        Group %in% c("Latinx", "Peruvian") ~ N_IAM,
        TRUE ~ NA_real_
      )
    )
}

if (anyNA(simulation_results$reference_sample_size)) {
  stop(
    "Unable to determine reference sample size for one or more rows."
  )
}

# ============================================================
# Extract ancestry-specific absolute bias
#
# Scenario 1:
#   retain all five ancestry-specific absolute-bias values.
#
# Scenarios 2/3:
#   additionally identify the absolute bias for VaryingAncestry,
#   which is the primary targeted outcome for those scenarios.
# ============================================================

for (reference_group in reference_groups) {

  simulation_results[[paste0(
    "uncorrected_absolute_bias_",
    reference_group
  )]] <- simulation_results[[paste0(
    "absolute_bias_",
    reference_group
  )]]

  simulation_results[[paste0(
    "linear_absolute_bias_eval_",
    reference_group
  )]] <- simulation_results[[paste0(
    "linear_absolute_bias_",
    reference_group
  )]]

  simulation_results[[paste0(
    "powerlaw_absolute_bias_eval_",
    reference_group
  )]] <- simulation_results[[paste0(
    "powerlaw_absolute_bias_",
    reference_group
  )]]
}

simulation_results$uncorrected_varying_absolute_bias <- NA_real_
simulation_results$linear_varying_absolute_bias <- NA_real_
simulation_results$powerlaw_varying_absolute_bias <- NA_real_

if ("VaryingAncestry" %in% names(simulation_results)) {

  varying_rows <- which(
    !is.na(simulation_results$VaryingAncestry)
  )

  for (row_index in varying_rows) {

    varying_group <-
      simulation_results$VaryingAncestry[row_index]

    if (!varying_group %in% reference_groups) {
      stop(
        "Unknown VaryingAncestry: ",
        varying_group
      )
    }

    simulation_results$uncorrected_varying_absolute_bias[row_index] <-
      simulation_results[[paste0(
        "absolute_bias_",
        varying_group
      )]][row_index]

    simulation_results$linear_varying_absolute_bias[row_index] <-
      simulation_results[[paste0(
        "linear_absolute_bias_",
        varying_group
      )]][row_index]

    simulation_results$powerlaw_varying_absolute_bias[row_index] <-
      simulation_results[[paste0(
        "powerlaw_absolute_bias_",
        varying_group
      )]][row_index]
  }
}

# ============================================================
# Convert correction methods to long format
# ============================================================

comparison_results <- simulation_results %>%
  select(
    all_of(required_metadata),
    all_of(optional_metadata),
    reference_sample_size,
    source_file,
    total_absolute_deviation,
    linear_total_absolute_deviation,
    powerlaw_total_absolute_deviation,
    uncorrected_varying_absolute_bias,
    linear_varying_absolute_bias,
    powerlaw_varying_absolute_bias
  ) %>%
  pivot_longer(
    cols = c(
      total_absolute_deviation,
      linear_total_absolute_deviation,
      powerlaw_total_absolute_deviation
    ),
    names_to = "Method",
    values_to = "total_absolute_bias"
  ) %>%
  mutate(
    Method = recode(
      Method,
      total_absolute_deviation = "Uncorrected",
      linear_total_absolute_deviation = "1/N",
      powerlaw_total_absolute_deviation = "a/N^b"
    ),
    Method = factor(
      Method,
      levels = method_levels
    ),
    varying_absolute_bias = case_when(
      Method == "Uncorrected" ~ uncorrected_varying_absolute_bias,
      Method == "1/N" ~ linear_varying_absolute_bias,
      Method == "a/N^b" ~ powerlaw_varying_absolute_bias,
      TRUE ~ NA_real_
    )
  )

# ============================================================
# Validate matched replicate structure
# ============================================================

condition_columns <- c(
  "Scenario",
  "Group",
  "reference_sample_size"
)

if ("VaryingAncestry" %in% names(comparison_results)) {
  condition_columns <- c(
    condition_columns,
    "VaryingAncestry"
  )
}

replicate_counts <- comparison_results %>%
  group_by(
    across(
      all_of(
        c(
          condition_columns,
          "Method"
        )
      )
    )
  ) %>%
  summarise(
    n_replicates = n(),
    .groups = "drop"
  )

replicate_count_check <- replicate_counts %>%
  group_by(
    across(
      all_of(condition_columns)
    )
  ) %>%
  summarise(
    n_methods = n_distinct(Method),
    n_unique_replicate_counts = n_distinct(n_replicates),
    .groups = "drop"
  )

if (any(
  replicate_count_check$n_methods != length(method_levels)
)) {
  stop(
    "One or more simulation conditions do not contain all three methods."
  )
}

if (any(
  replicate_count_check$n_unique_replicate_counts != 1L
)) {
  stop(
    "Correction methods have unequal replicate counts for one or more conditions."
  )
}

# ============================================================
# Validate one observation per replicate and method
# ============================================================

duplicate_check <- comparison_results %>%
  count(
    across(
      all_of(
        c(
          condition_columns,
          "Replicate",
          "Method"
        )
      )
    ),
    name = "n"
  ) %>%
  filter(
    n != 1L
  )

if (nrow(duplicate_check) > 0L) {
  stop(
    "Duplicate replicate × method observations were found."
  )
}

# ============================================================
# Summarize correction performance
# ============================================================

comparison_summary <- comparison_results %>%
  group_by(
    across(
      all_of(
        c(
          condition_columns,
          "Method"
        )
      )
    )
  ) %>%
  summarise(
    mean_total_absolute_bias = mean(
      total_absolute_bias
    ),
    median_total_absolute_bias = median(
      total_absolute_bias
    ),
    sd_total_absolute_bias = sd(
      total_absolute_bias
    ),
    q025_total_absolute_bias = quantile(
      total_absolute_bias,
      0.025
    ),
    q975_total_absolute_bias = quantile(
      total_absolute_bias,
      0.975
    ),
    n_replicates = n(),
    .groups = "drop"
  ) %>%
  arrange(
    Scenario,
    Group,
    reference_sample_size,
    Method
  )

# ============================================================
# Summarize the targeted varying-ancestry absolute bias
#
# This is the primary outcome for Scenarios 2/3.
# Scenario 1 has no VaryingAncestry, so it is excluded here.
# ============================================================

varying_ancestry_summary <- comparison_results %>%
  filter(
    !is.na(VaryingAncestry),
    !is.na(varying_absolute_bias)
  ) %>%
  group_by(
    Scenario,
    Group,
    VaryingAncestry,
    reference_sample_size,
    Method
  ) %>%
  summarise(
    mean_varying_absolute_bias = mean(
      varying_absolute_bias
    ),
    median_varying_absolute_bias = median(
      varying_absolute_bias
    ),
    sd_varying_absolute_bias = sd(
      varying_absolute_bias
    ),
    q025_varying_absolute_bias = quantile(
      varying_absolute_bias,
      0.025
    ),
    q975_varying_absolute_bias = quantile(
      varying_absolute_bias,
      0.975
    ),
    n_replicates = n(),
    .groups = "drop"
  ) %>%
  arrange(
    Scenario,
    Group,
    VaryingAncestry,
    reference_sample_size,
    Method
  )

uncorrected_varying_summary <- varying_ancestry_summary %>%
  filter(
    Method == "Uncorrected"
  ) %>%
  select(
    Scenario,
    Group,
    VaryingAncestry,
    reference_sample_size,
    uncorrected_mean_varying_absolute_bias =
      mean_varying_absolute_bias
  )

varying_ancestry_summary <- varying_ancestry_summary %>%
  left_join(
    uncorrected_varying_summary,
    by = c(
      "Scenario",
      "Group",
      "VaryingAncestry",
      "reference_sample_size"
    )
  ) %>%
  mutate(
    mean_varying_bias_reduction =
      uncorrected_mean_varying_absolute_bias -
      mean_varying_absolute_bias,
    relative_varying_bias_reduction = ifelse(
      uncorrected_mean_varying_absolute_bias > 0,
      mean_varying_bias_reduction /
        uncorrected_mean_varying_absolute_bias,
      NA_real_
    )
  )

# ============================================================
# Calculate improvement relative to uncorrected Summix2
# ============================================================

uncorrected_summary <- comparison_summary %>%
  filter(
    Method == "Uncorrected"
  ) %>%
  select(
    all_of(condition_columns),
    uncorrected_mean_total_absolute_bias =
      mean_total_absolute_bias
  )

comparison_summary <- comparison_summary %>%
  left_join(
    uncorrected_summary,
    by = condition_columns
  ) %>%
  mutate(
    mean_absolute_bias_reduction =
      uncorrected_mean_total_absolute_bias -
      mean_total_absolute_bias,
    relative_bias_reduction = ifelse(
      uncorrected_mean_total_absolute_bias > 0,
      mean_absolute_bias_reduction /
        uncorrected_mean_total_absolute_bias,
      NA_real_
    )
  )

# ============================================================
# Save evaluation datasets
# ============================================================

saveRDS(
  comparison_results,
  output_file,
  compress = TRUE
)

saveRDS(
  comparison_summary,
  summary_file,
  compress = TRUE
)

saveRDS(
  varying_ancestry_summary,
  varying_summary_file,
  compress = TRUE
)

message(
  "Saved replicate-level correction comparison to: ",
  normalizePath(
    output_file,
    mustWork = FALSE
  ),
  "\nSaved summarized correction comparison to: ",
  normalizePath(
    summary_file,
    mustWork = FALSE
  ),
  "\nSaved varying-ancestry correction summary to: ",
  normalizePath(
    varying_summary_file,
    mustWork = FALSE
  ),
  "\nMethods: ",
  paste(
    method_levels,
    collapse = ", "
  )
)
