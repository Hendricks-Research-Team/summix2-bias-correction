options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages({
  library(dplyr)
})

# ============================================================
# Summarize power-law correction parameters
# ============================================================

# -----------------------------
# Configuration
# -----------------------------
input_file <- Sys.getenv(
  "INPUT_FILE",
  "power_law_fit_parameters.rds"
)

output_file <- Sys.getenv(
  "OUTPUT_FILE",
  "power_law_parameter_summary.rds"
)

ancestries <- c("AFR", "EAS", "EUR", "IAM", "SAS")

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

# -----------------------------
# Read fitted parameters
# -----------------------------
fit_summary <- readRDS(input_file)

required_columns <- c(
  "ProportionID",
  "ancestry",
  "a",
  "b",
  "converged"
)

missing_columns <- setdiff(
  required_columns,
  names(fit_summary)
)

if (length(missing_columns) > 0L) {
  stop(
    "Input file is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (any(!fit_summary$ancestry %in% ancestries)) {
  stop("The input contains an unexpected ancestry label.")
}

if (anyNA(fit_summary$converged)) {
  stop("The converged column contains missing values.")
}

if (anyDuplicated(fit_summary[c("ProportionID", "ancestry")])) {
  stop("Duplicate ProportionID × ancestry fits were found.")
}

# -----------------------------
# Validate fit coverage
# -----------------------------
fit_counts <- fit_summary %>%
  group_by(ancestry) %>%
  summarise(
    n_fits = n(),
    n_converged = sum(converged),
    n_failed = sum(!converged),
    .groups = "drop"
  )

missing_ancestries <- setdiff(
  ancestries,
  fit_counts$ancestry
)

if (length(missing_ancestries) > 0L) {
  stop(
    "No fitted models were found for: ",
    paste(missing_ancestries, collapse = ", ")
  )
}

# -----------------------------
# Summarize converged fits
# -----------------------------
parameter_summary <- fit_summary %>%
  filter(
    converged,
    !is.na(a),
    !is.na(b)
  ) %>%
  group_by(ancestry) %>%
  summarise(
    median_a = median(a),
    median_b = median(b),
    .groups = "drop"
  ) %>%
  right_join(
    fit_counts,
    by = "ancestry"
  ) %>%
  mutate(
    convergence_rate = n_converged / n_fits
  ) %>%
  arrange(
    match(ancestry, ancestries)
  )

if (anyNA(parameter_summary$median_a) ||
    anyNA(parameter_summary$median_b)) {
  stop(
    "At least one ancestry has no converged fits available ",
    "for parameter summarization."
  )
}

# -----------------------------
# Save results
# -----------------------------
saveRDS(
  parameter_summary,
  output_file,
  compress = TRUE
)

message(
  "Saved power-law parameter summary to: ",
  normalizePath(output_file, mustWork = FALSE)
)

print(parameter_summary)
