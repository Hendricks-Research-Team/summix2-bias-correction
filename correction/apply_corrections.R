# ============================================================
# Apply bias corrections to simulation results
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

correction_functions_file <- Sys.getenv(
  "CORRECTION_FUNCTIONS",
  "correction_functions.R"
)

parameter_file <- Sys.getenv(
  "PARAMETER_FILE",
  "../results/power_law/power_law_parameter_summary.rds"
)

scenario1_dir <- Sys.getenv(
  "SCENARIO1_DIR",
  "../results/simulation/scenario_1_equal_sizes"
)

scenario23_dir <- Sys.getenv(
  "SCENARIO23_DIR",
  "../results/simulation/scenario_2_3_one_varying"
)

output_dir <- Sys.getenv(
  "OUTPUT_DIR",
  "../results/correction"
)

if (!file.exists(correction_functions_file)) {
  stop("Correction functions file does not exist: ", correction_functions_file)
}

if (!file.exists(parameter_file)) {
  stop("Power-law parameter file does not exist: ", parameter_file)
}

for (input_dir in c(scenario1_dir, scenario23_dir)) {
  if (!dir.exists(input_dir)) {
    stop("Simulation input directory does not exist: ", input_dir)
  }
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source(correction_functions_file)

parameter_table <- readRDS(parameter_file)
parameter_table <- validate_power_law_parameters(parameter_table)

estimate_columns <- reference_groups
reference_size_columns <- paste0("N_", reference_groups)
truth_columns <- paste0("true_", reference_groups)

required_columns <- c(
  "Scenario",
  "Group",
  "Replicate",
  estimate_columns,
  reference_size_columns,
  truth_columns
)

get_groups_to_correct <- function(simulation_row) {
  scenario <- simulation_row$Scenario[[1]]

  if (scenario == "scenario_1_equal_sizes") {
    return(reference_groups)
  }

  if (scenario %in% c(
    "scenario_2_one_varying_fixed_500",
    "scenario_3_one_varying_fixed_150"
  )) {
    if (!"VaryingAncestry" %in% names(simulation_row)) {
      stop("Scenario 2/3 result is missing VaryingAncestry.")
    }

    varying_group <- simulation_row$VaryingAncestry[[1]]

    if (!varying_group %in% reference_groups) {
      stop("Unknown varying reference group: ", varying_group)
    }

    return(varying_group)
  }

  stop("Unknown simulation scenario: ", scenario)
}

apply_corrections_to_results <- function(simulation_results) {
  missing_columns <- setdiff(required_columns, names(simulation_results))

  if (length(missing_columns) > 0L) {
    stop(
      "Simulation results are missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (nrow(simulation_results) == 0L) {
    stop("Simulation result file contains no rows.")
  }

  linear_corrected_matrix <- t(
    vapply(
      seq_len(nrow(simulation_results)),
      function(row_index) {
        simulation_row <- simulation_results[
          row_index,
          ,
          drop = FALSE
        ]

        estimates <- as.numeric(
          simulation_row[
            ,
            estimate_columns,
            drop = TRUE
          ]
        )
        names(estimates) <- reference_groups

        reference_sizes <- as.numeric(
          simulation_row[
            ,
            reference_size_columns,
            drop = TRUE
          ]
        )
        names(reference_sizes) <- reference_groups

        groups_to_correct <- get_groups_to_correct(simulation_row)

        apply_linear_correction(
          estimates = estimates,
          reference_sizes = reference_sizes,
          groups_to_correct = groups_to_correct
        )
      },
      numeric(length(reference_groups))
    )
  )

  colnames(linear_corrected_matrix) <- paste0(
    "linear_corrected_",
    reference_groups
  )

  powerlaw_corrected_matrix <- t(
    vapply(
      seq_len(nrow(simulation_results)),
      function(row_index) {
        simulation_row <- simulation_results[
          row_index,
          ,
          drop = FALSE
        ]

        estimates <- as.numeric(
          simulation_row[
            ,
            estimate_columns,
            drop = TRUE
          ]
        )
        names(estimates) <- reference_groups

        reference_sizes <- as.numeric(
          simulation_row[
            ,
            reference_size_columns,
            drop = TRUE
          ]
        )
        names(reference_sizes) <- reference_groups

        groups_to_correct <- get_groups_to_correct(simulation_row)

        apply_power_law_correction(
          estimates = estimates,
          reference_sizes = reference_sizes,
          parameter_table = parameter_table,
          groups_to_correct = groups_to_correct
        )
      },
      numeric(length(reference_groups))
    )
  )

  colnames(powerlaw_corrected_matrix) <- paste0(
    "powerlaw_corrected_",
    reference_groups
  )

  corrected_results <- bind_cols(
    simulation_results,
    as.data.frame(linear_corrected_matrix),
    as.data.frame(powerlaw_corrected_matrix)
  )

  for (reference_group in reference_groups) {
    truth_column <- paste0("true_", reference_group)

    linear_corrected_column <- paste0(
      "linear_corrected_",
      reference_group
    )
    linear_bias_column <- paste0(
      "linear_bias_",
      reference_group
    )
    linear_absolute_bias_column <- paste0(
      "linear_absolute_bias_",
      reference_group
    )

    powerlaw_corrected_column <- paste0(
      "powerlaw_corrected_",
      reference_group
    )
    powerlaw_bias_column <- paste0(
      "powerlaw_bias_",
      reference_group
    )
    powerlaw_absolute_bias_column <- paste0(
      "powerlaw_absolute_bias_",
      reference_group
    )

    corrected_results[[linear_bias_column]] <-
      corrected_results[[linear_corrected_column]] -
      corrected_results[[truth_column]]

    corrected_results[[linear_absolute_bias_column]] <-
      abs(corrected_results[[linear_bias_column]])

    corrected_results[[powerlaw_bias_column]] <-
      corrected_results[[powerlaw_corrected_column]] -
      corrected_results[[truth_column]]

    corrected_results[[powerlaw_absolute_bias_column]] <-
      abs(corrected_results[[powerlaw_bias_column]])
  }

  corrected_results$linear_total_absolute_deviation <- rowSums(
    corrected_results[
      ,
      paste0("linear_absolute_bias_", reference_groups),
      drop = FALSE
    ]
  )

  corrected_results$powerlaw_total_absolute_deviation <- rowSums(
    corrected_results[
      ,
      paste0("powerlaw_absolute_bias_", reference_groups),
      drop = FALSE
    ]
  )

  corrected_results
}

process_simulation_directory <- function(input_dir, output_subdir) {
  input_files <- list.files(
    input_dir,
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(input_files) == 0L) {
    stop("No simulation .rds files were found in: ", input_dir)
  }

  output_root <- file.path(output_dir, output_subdir)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  normalized_input_dir <- normalizePath(input_dir, mustWork = TRUE)

  message(
    "Processing ",
    length(input_files),
    " files from: ",
    input_dir
  )

  for (input_file in input_files) {
    simulation_results <- readRDS(input_file)
    corrected_results <- apply_corrections_to_results(simulation_results)

    normalized_input_file <- normalizePath(
      input_file,
      mustWork = TRUE
    )

    relative_path <- substring(
      normalized_input_file,
      nchar(normalized_input_dir) + 2L
    )

    output_file <- file.path(
      output_root,
      relative_path
    )

    dir.create(
      dirname(output_file),
      recursive = TRUE,
      showWarnings = FALSE
    )

    saveRDS(
      corrected_results,
      output_file,
      compress = TRUE
    )

    message("Saved: ", output_file)
  }
}

process_simulation_directory(
  input_dir = scenario1_dir,
  output_subdir = "scenario_1_equal_sizes"
)

process_simulation_directory(
  input_dir = scenario23_dir,
  output_subdir = "scenario_2_3_one_varying"
)

message("Bias-correction application complete.")
