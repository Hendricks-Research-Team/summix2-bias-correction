options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ============================================================
# Summarize Dirichlet simulation bias
# ============================================================

# -----------------------------
# Configuration
# -----------------------------
input_dir <- Sys.getenv(
  "INPUT_DIR",
  "/home/math/tifours/summix_sim_outputs_100k/power_law_parameter_estimation"
)

output_file <- Sys.getenv(
  "OUTPUT_FILE",
  "dirichlet_bias_summary_long.rds"
)

ancestries <- c("AFR", "EAS", "EUR", "IAM", "SAS")

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir)
}

# -----------------------------
# Locate simulation outputs
# -----------------------------
input_files <- list.files(
  input_dir,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(input_files) == 0L) {
  stop("No .rds simulation output files were found in: ", input_dir)
}

message("Found ", length(input_files), " simulation output files.")

# -----------------------------
# Read simulation outputs
# -----------------------------
simulation_results <- lapply(
  input_files,
  function(input_file) {
    result <- readRDS(input_file)

    required_columns <- c(
      "ProportionID",
      "N_AFR",
      "N_EAS",
      "N_EUR",
      "N_IAM",
      "N_SAS",
      paste0("bias_", ancestries)
    )

    missing_columns <- setdiff(
      required_columns,
      names(result)
    )

    if (length(missing_columns) > 0L) {
      stop(
        "File is missing required columns: ",
        input_file,
        "\nMissing columns: ",
        paste(missing_columns, collapse = ", ")
      )
    }

    result
  }
)

simulation_results <- bind_rows(simulation_results)

# -----------------------------
# Validate sample-size structure
# -----------------------------
reference_size_columns <- paste0("N_", ancestries)

equal_reference_sizes <- apply(
  simulation_results[, reference_size_columns, drop = FALSE],
  1,
  function(reference_sizes) {
    length(unique(reference_sizes)) == 1L
  }
)

if (!all(equal_reference_sizes)) {
  stop(
    "One or more Dirichlet simulation rows do not have equal ",
    "reference sample sizes across all five reference groups."
  )
}

simulation_results$sample_size <- simulation_results$N_AFR

# -----------------------------
# Summarize signed and absolute bias
# -----------------------------
bias_summary <- simulation_results %>%
  group_by(
    ProportionID,
    sample_size
  ) %>%
  summarise(
    across(
      all_of(paste0("bias_", ancestries)),
      list(
        mean_signed_bias = ~ mean(.x, na.rm = TRUE),
        mean_absolute_bias = ~ mean(abs(.x), na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    n_replicates = n(),
    .groups = "drop"
  )

bias_summary_long <- bias_summary %>%
  pivot_longer(
    cols = matches("^bias_(AFR|EAS|EUR|IAM|SAS)_(mean_signed_bias|mean_absolute_bias)$"),
    names_to = c("ancestry", ".value"),
    names_pattern = "^bias_(AFR|EAS|EUR|IAM|SAS)_(mean_signed_bias|mean_absolute_bias)$"
  ) %>%
  arrange(
    ProportionID,
    ancestry,
    sample_size
  )

# -----------------------------
# Validate summary
# -----------------------------
expected_rows <- length(unique(simulation_results$ProportionID)) *
  length(unique(simulation_results$sample_size)) *
  length(ancestries)

if (nrow(bias_summary_long) != expected_rows) {
  stop(
    "Unexpected number of rows in the summarized bias dataset. ",
    "Expected ", expected_rows,
    " but found ", nrow(bias_summary_long), "."
  )
}

if (anyNA(bias_summary_long$mean_signed_bias)) {
  stop("The summarized dataset contains missing mean signed bias values.")
}

if (anyNA(bias_summary_long$mean_absolute_bias)) {
  stop("The summarized dataset contains missing mean absolute bias values.")
}

# -----------------------------
# Save results
# -----------------------------
saveRDS(
  bias_summary_long,
  output_file,
  compress = TRUE
)

message(
  "Saved ",
  nrow(bias_summary_long),
  " summarized bias rows to: ",
  normalizePath(output_file, mustWork = FALSE)
)
