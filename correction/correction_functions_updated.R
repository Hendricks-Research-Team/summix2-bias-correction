# ============================================================
# Shared bias-correction functions
# ============================================================

reference_groups <- c("AFR", "EAS", "EUR", "IAM", "SAS")
K <- length(reference_groups)
default_threshold <- 1 / K

validate_reference_vector <- function(x, vector_name, positive = FALSE) {
  if (!is.numeric(x)) stop(vector_name, " must be numeric.")
  if (is.null(names(x))) stop(vector_name, " must be a named vector.")

  missing_groups <- setdiff(reference_groups, names(x))

  if (length(missing_groups) > 0L) {
    stop(
      vector_name,
      " is missing reference groups: ",
      paste(missing_groups, collapse = ", ")
    )
  }

  x <- x[reference_groups]

  if (anyNA(x)) stop(vector_name, " contains missing values.")

  if (positive && any(x <= 0)) {
    stop(vector_name, " must contain only positive values.")
  }

  x
}


# ============================================================
# Apply correction direction to selected reference groups only
# ============================================================

apply_direction_rule <- function(
    estimates,
    correction_magnitude,
    groups_to_correct = reference_groups,
    threshold = default_threshold
) {
  estimates <- validate_reference_vector(estimates, "estimates")

  correction_magnitude <- validate_reference_vector(
    correction_magnitude,
    "correction_magnitude"
  )

  if (!is.numeric(threshold) ||
      length(threshold) != 1L ||
      is.na(threshold)) {
    stop("threshold must be a single numeric value.")
  }

  if (!all(groups_to_correct %in% reference_groups)) {
    stop("groups_to_correct contains an unknown reference group.")
  }

  corrected <- estimates

  for (reference_group in groups_to_correct) {
    estimate <- estimates[reference_group]
    magnitude <- correction_magnitude[reference_group]

    if (estimate < threshold) {
      corrected[reference_group] <- estimate - magnitude
    } else if (estimate > threshold) {
      corrected[reference_group] <- estimate + magnitude
    } else {
      corrected[reference_group] <- estimate
    }
  }

  corrected
}


# ============================================================
# Post-correction handling
#
# 1. Truncate negative values to zero.
# 2. Renormalize the five estimates to sum to one.
# 3. If truncation produces an all-zero vector, assign 1/K
#    to every reference group.
#
# Positive values are NOT clamped to one before renormalization.
# ============================================================

renormalize_corrected_estimates <- function(corrected_estimates) {
  corrected_estimates <- validate_reference_vector(
    corrected_estimates,
    "corrected_estimates"
  )

  corrected_estimates <- pmax(corrected_estimates, 0)

  corrected_sum <- sum(corrected_estimates)

  if (corrected_sum == 0) {
    corrected_estimates <- rep(
      1 / length(reference_groups),
      length(reference_groups)
    )
  } else {
    corrected_estimates <- corrected_estimates / corrected_sum
  }

  names(corrected_estimates) <- reference_groups
  corrected_estimates
}


# ============================================================
# Linear 1/N correction
# ============================================================

apply_linear_correction <- function(
    estimates,
    reference_sizes,
    groups_to_correct = reference_groups,
    threshold = default_threshold
) {
  estimates <- validate_reference_vector(estimates, "estimates")

  reference_sizes <- validate_reference_vector(
    reference_sizes,
    "reference_sizes",
    positive = TRUE
  )

  correction_magnitude <- 1 / reference_sizes

  corrected <- apply_direction_rule(
    estimates = estimates,
    correction_magnitude = correction_magnitude,
    groups_to_correct = groups_to_correct,
    threshold = threshold
  )

  renormalize_corrected_estimates(corrected)
}


# ============================================================
# Validate power-law parameter table
# ============================================================

validate_power_law_parameters <- function(parameter_table) {
  required_columns <- c("ancestry", "median_a", "median_b")
  missing_columns <- setdiff(required_columns, names(parameter_table))

  if (length(missing_columns) > 0L) {
    stop(
      "Power-law parameter table is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (anyDuplicated(parameter_table$ancestry)) {
    stop(
      "Power-law parameter table contains duplicate reference-group rows."
    )
  }

  missing_groups <- setdiff(reference_groups, parameter_table$ancestry)

  if (length(missing_groups) > 0L) {
    stop(
      "Power-law parameter table is missing reference groups: ",
      paste(missing_groups, collapse = ", ")
    )
  }

  parameter_table <- parameter_table[
    match(reference_groups, parameter_table$ancestry),
    ,
    drop = FALSE
  ]

  if (anyNA(parameter_table$median_a) ||
      anyNA(parameter_table$median_b)) {
    stop(
      "Power-law parameter table contains missing median_a or median_b values."
    )
  }

  if (any(parameter_table$median_a < 0)) {
    stop("median_a must be non-negative.")
  }

  if (any(parameter_table$median_b < 0)) {
    stop("median_b must be non-negative.")
  }

  parameter_table
}


# ============================================================
# Nonlinear power-law correction
#
# correction magnitude = a / N^b
# ============================================================

apply_power_law_correction <- function(
    estimates,
    reference_sizes,
    parameter_table,
    groups_to_correct = reference_groups,
    threshold = default_threshold
) {
  estimates <- validate_reference_vector(estimates, "estimates")

  reference_sizes <- validate_reference_vector(
    reference_sizes,
    "reference_sizes",
    positive = TRUE
  )

  parameter_table <- validate_power_law_parameters(parameter_table)

  a_parameters <- setNames(
    parameter_table$median_a,
    parameter_table$ancestry
  )

  b_parameters <- setNames(
    parameter_table$median_b,
    parameter_table$ancestry
  )

  correction_magnitude <-
    a_parameters / (reference_sizes ^ b_parameters)

  corrected <- apply_direction_rule(
    estimates = estimates,
    correction_magnitude = correction_magnitude,
    groups_to_correct = groups_to_correct,
    threshold = threshold
  )

  renormalize_corrected_estimates(corrected)
}
