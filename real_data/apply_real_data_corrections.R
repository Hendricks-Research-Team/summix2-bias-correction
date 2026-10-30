# ============================================================
# Apply bias corrections to real-data stress-test results
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ============================================================
# Configuration
# ============================================================

correction_functions_file <- Sys.getenv(
  "CORRECTION_FUNCTIONS",
  "../correction/correction_functions.R"
)

parameter_file <- Sys.getenv(
  "PARAMETER_FILE",
  "../results/power_law/power_law_parameter_summary.rds"
)

real_data_dir <- Sys.getenv(
  "REAL_DATA_DIR",
  "../results/real_data"
)

joint_input_file <- Sys.getenv(
  "JOINT_INPUT_FILE",
  file.path(
    real_data_dir,
    "real_data_jointly_limited_uncorrected.rds"
  )
)

vary_afr_input_file <- Sys.getenv(
  "VARY_AFR_INPUT_FILE",
  file.path(
    real_data_dir,
    "real_data_vary_AFR_uncorrected.rds"
  )
)

joint_output_file <- Sys.getenv(
  "JOINT_OUTPUT_FILE",
  file.path(
    real_data_dir,
    "real_data_jointly_limited_corrected.rds"
  )
)

vary_afr_output_file <- Sys.getenv(
  "VARY_AFR_OUTPUT_FILE",
  file.path(
    real_data_dir,
    "real_data_vary_AFR_corrected.rds"
  )
)

for (
  input_file in c(
    correction_functions_file,
    parameter_file,
    joint_input_file,
    vary_afr_input_file
  )
) {

  if (!file.exists(input_file)) {
    stop(
      "Required input file does not exist: ",
      input_file
    )
  }
}

source(
  correction_functions_file
)

parameter_table <- readRDS(
  parameter_file
)

parameter_table <- validate_power_law_parameters(
  parameter_table
)

# ============================================================
# Required columns
# ============================================================

estimate_columns <- reference_groups

reference_size_columns <- paste0(
  "N_",
  reference_groups
)

required_columns <- c(
  estimate_columns,
  reference_size_columns
)

# ============================================================
# Validate real-data result structure
# ============================================================

validate_real_data_results <- function(
    real_data_results,
    data_name
) {

  missing_columns <- setdiff(
    required_columns,
    names(real_data_results)
  )

  if (length(missing_columns) > 0L) {
    stop(
      data_name,
      " is missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  if (nrow(real_data_results) == 0L) {
    stop(
      data_name,
      " contains no rows."
    )
  }

  invisible(TRUE)
}

# ============================================================
# Apply corrections to one real-data result row
# ============================================================

apply_real_data_row_corrections <- function(
    result_row,
    groups_to_correct
) {

  estimates <- as.numeric(
    result_row[
      ,
      estimate_columns,
      drop = TRUE
    ]
  )

  names(estimates) <- reference_groups

  reference_sizes <- as.numeric(
    result_row[
      ,
      reference_size_columns,
      drop = TRUE
    ]
  )

  names(reference_sizes) <- reference_groups

  linear_corrected <- apply_linear_correction(
    estimates = estimates,
    reference_sizes = reference_sizes,
    groups_to_correct = groups_to_correct
  )

  powerlaw_corrected <- apply_power_law_correction(
    estimates = estimates,
    reference_sizes = reference_sizes,
    parameter_table = parameter_table,
    groups_to_correct = groups_to_correct
  )

  list(
    linear = linear_corrected,
    powerlaw = powerlaw_corrected
  )
}

# ============================================================
# Apply corrections to a full real-data result table
# ============================================================

apply_real_data_corrections <- function(
    real_data_results,
    groups_to_correct,
    data_name
) {

  validate_real_data_results(
    real_data_results,
    data_name
  )

  if (
    !all(
      groups_to_correct %in%
        reference_groups
    )
  ) {
    stop(
      "groups_to_correct contains an unknown reference group."
    )
  }

  n_rows <- nrow(
    real_data_results
  )

  linear_corrected_matrix <- t(
    vapply(
      seq_len(
        n_rows
      ),
      function(row_index) {

        corrected <- apply_real_data_row_corrections(
          result_row =
            real_data_results[
              row_index,
              ,
              drop = FALSE
            ],
          groups_to_correct =
            groups_to_correct
        )

        corrected$linear
      },
      numeric(
        length(
          reference_groups
        )
      )
    )
  )

  powerlaw_corrected_matrix <- t(
    vapply(
      seq_len(
        n_rows
      ),
      function(row_index) {

        corrected <- apply_real_data_row_corrections(
          result_row =
            real_data_results[
              row_index,
              ,
              drop = FALSE
            ],
          groups_to_correct =
            groups_to_correct
        )

        corrected$powerlaw
      },
      numeric(
        length(
          reference_groups
        )
      )
    )
  )

  colnames(
    linear_corrected_matrix
  ) <- paste0(
    "linear_corrected_",
    reference_groups
  )

  colnames(
    powerlaw_corrected_matrix
  ) <- paste0(
    "powerlaw_corrected_",
    reference_groups
  )

  corrected_results <- bind_cols(
    real_data_results,
    as.data.frame(
      linear_corrected_matrix
    ),
    as.data.frame(
      powerlaw_corrected_matrix
    )
  )

  # ----------------------------------------------------------
  # Validate corrected mixing proportions
  # ----------------------------------------------------------

  linear_columns <- paste0(
    "linear_corrected_",
    reference_groups
  )

  powerlaw_columns <- paste0(
    "powerlaw_corrected_",
    reference_groups
  )

  linear_sums <- rowSums(
    corrected_results[
      ,
      linear_columns,
      drop = FALSE
    ]
  )

  powerlaw_sums <- rowSums(
    corrected_results[
      ,
      powerlaw_columns,
      drop = FALSE
    ]
  )

  if (
    any(
      abs(
        linear_sums - 1
      ) > 1e-10
    )
  ) {
    stop(
      data_name,
      ": one or more 1/N corrected vectors do not sum to one."
    )
  }

  if (
    any(
      abs(
        powerlaw_sums - 1
      ) > 1e-10
    )
  ) {
    stop(
      data_name,
      ": one or more power-law corrected vectors do not sum to one."
    )
  }

  if (
    any(
      as.matrix(
        corrected_results[
          ,
          linear_columns,
          drop = FALSE
        ]
      ) < 0
    ) ||
      any(
        as.matrix(
          corrected_results[
            ,
            linear_columns,
            drop = FALSE
          ]
        ) > 1
      )
  ) {
    stop(
      data_name,
      ": invalid 1/N corrected mixing proportions detected."
    )
  }

  if (
    any(
      as.matrix(
        corrected_results[
          ,
          powerlaw_columns,
          drop = FALSE
        ]
      ) < 0
    ) ||
      any(
        as.matrix(
          corrected_results[
            ,
            powerlaw_columns,
            drop = FALSE
          ]
        ) > 1
      )
  ) {
    stop(
      data_name,
      ": invalid power-law corrected mixing proportions detected."
    )
  }

  corrected_results
}

# ============================================================
# Read uncorrected real-data stress-test results
# ============================================================

joint_results <- readRDS(
  joint_input_file
)

vary_afr_results <- readRDS(
  vary_afr_input_file
)

# ============================================================
# Stress test A
# Jointly limited reference groups
#
# All five reference groups are directly corrected.
# ============================================================

joint_corrected <- apply_real_data_corrections(
  real_data_results = joint_results,
  groups_to_correct = reference_groups,
  data_name = "Jointly limited real-data results"
)

# ============================================================
# Stress test B
# AFR varies while other reference groups remain fixed
#
# Only AFR is directly corrected.
# All five estimates are then clamped and renormalized by
# the shared correction functions.
# ============================================================

vary_afr_corrected <- apply_real_data_corrections(
  real_data_results = vary_afr_results,
  groups_to_correct = "AFR",
  data_name = "AFR-varying real-data results"
)

# ============================================================
# Record correction targeting
# ============================================================

joint_corrected$CorrectionTarget <-
  "all_reference_groups"

vary_afr_corrected$CorrectionTarget <-
  "AFR_only"

# ============================================================
# Save corrected real-data results
# ============================================================

saveRDS(
  joint_corrected,
  joint_output_file,
  compress = TRUE
)

saveRDS(
  vary_afr_corrected,
  vary_afr_output_file,
  compress = TRUE
)

message(
  "Saved jointly limited corrected results to: ",
  normalizePath(
    joint_output_file,
    mustWork = FALSE
  ),
  "\nSaved AFR-varying corrected results to: ",
  normalizePath(
    vary_afr_output_file,
    mustWork = FALSE
  )
)
