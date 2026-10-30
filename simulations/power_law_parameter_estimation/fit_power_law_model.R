options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages({
  library(dplyr)
})

# ============================================================
# Fit power-law models to summarized Dirichlet simulation bias
# ============================================================

input_file <- Sys.getenv(
  "INPUT_FILE",
  "dirichlet_bias_summary_long.rds"
)

output_file <- Sys.getenv(
  "OUTPUT_FILE",
  "power_law_fit_parameters.rds"
)

ancestries <- c("AFR", "EAS", "EUR", "IAM", "SAS")

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

bias_summary_long <- readRDS(input_file)

required_columns <- c(
  "ProportionID",
  "sample_size",
  "ancestry",
  "mean_signed_bias",
  "mean_absolute_bias"
)

missing_columns <- setdiff(
  required_columns,
  names(bias_summary_long)
)

if (length(missing_columns) > 0L) {
  stop(
    "Input file is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyNA(bias_summary_long[, required_columns])) {
  stop("Required input columns contain missing values.")
}

if (any(!bias_summary_long$ancestry %in% ancestries)) {
  stop("The input contains an unexpected ancestry label.")
}

if (any(
  bias_summary_long$sample_size < 1L |
    bias_summary_long$sample_size > 1000L
)) {
  stop("sample_size must be between 1 and 1000.")
}

if (any(bias_summary_long$mean_absolute_bias < 0)) {
  stop("mean_absolute_bias cannot contain negative values.")
}

trajectory_check <- bias_summary_long %>%
  group_by(ProportionID, ancestry) %>%
  summarise(
    n_sample_sizes = n_distinct(sample_size),
    min_sample_size = min(sample_size),
    max_sample_size = max(sample_size),
    .groups = "drop"
  )

if (any(
  trajectory_check$n_sample_sizes != 1000L |
    trajectory_check$min_sample_size != 1L |
    trajectory_check$max_sample_size != 1000L
)) {
  stop(
    "One or more ProportionID × ancestry combinations do not ",
    "contain the complete N = 1,...,1000 trajectory."
  )
}

# -----------------------------
# Power-law model
# -----------------------------
power_law_model <- function(N, a, b) {
  a / N^b
}

# -----------------------------
# Fit one model per proportion × ancestry
# -----------------------------
fit_groups <- split(
  bias_summary_long,
  interaction(
    bias_summary_long$ProportionID,
    bias_summary_long$ancestry,
    drop = TRUE
  )
)

fit_results <- lapply(
  fit_groups,
  function(subset_df) {
    subset_df <- subset_df[
      order(subset_df$sample_size),
      ,
      drop = FALSE
    ]

    proportion_id <- unique(subset_df$ProportionID)
    ancestry <- unique(subset_df$ancestry)

    warning_messages <- character(0)

    fit <- tryCatch(
      withCallingHandlers(
        nls(
          mean_absolute_bias ~ power_law_model(sample_size, a, b),
          data = subset_df,
          start = list(a = 1, b = 0.5),
          control = nls.control(warnOnly = TRUE)
        ),
        warning = function(w) {
          warning_messages <<- c(
            warning_messages,
            conditionMessage(w)
          )
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )

    warning_message <- if (length(warning_messages) > 0L) {
      paste(unique(warning_messages), collapse = " | ")
    } else {
      NA_character_
    }

    if (inherits(fit, "error")) {
      return(
        data.frame(
          ProportionID = proportion_id,
          ancestry = ancestry,
          a = NA_real_,
          b = NA_real_,
          converged = FALSE,
          iterations = NA_integer_,
          stop_code = NA_integer_,
          stop_message = conditionMessage(fit),
          warning_message = warning_message,
          n_sample_sizes = nrow(subset_df),
          stringsAsFactors = FALSE
        )
      )
    }

    coefficients <- coef(fit)
    convergence_info <- fit$convInfo

    data.frame(
      ProportionID = proportion_id,
      ancestry = ancestry,
      a = unname(coefficients["a"]),
      b = unname(coefficients["b"]),
      converged = isTRUE(convergence_info$isConv),
      iterations = convergence_info$finIter,
      stop_code = convergence_info$stopCode,
      stop_message = convergence_info$stopMessage,
      warning_message = warning_message,
      n_sample_sizes = nrow(subset_df),
      stringsAsFactors = FALSE
    )
  }
)

fit_summary <- bind_rows(fit_results) %>%
  arrange(ProportionID, match(ancestry, ancestries))

expected_fits <- length(unique(bias_summary_long$ProportionID)) *
  length(ancestries)

if (nrow(fit_summary) != expected_fits) {
  stop(
    "Unexpected number of fitted models. Expected ",
    expected_fits,
    " but found ",
    nrow(fit_summary),
    "."
  )
}

if (anyDuplicated(fit_summary[c("ProportionID", "ancestry")])) {
  stop("Duplicate ProportionID × ancestry fits were produced.")
}

saveRDS(
  fit_summary,
  output_file,
  compress = TRUE
)

message(
  "Saved ",
  nrow(fit_summary),
  " power-law model fits to: ",
  normalizePath(output_file, mustWork = FALSE),
  "\nConverged fits: ",
  sum(fit_summary$converged),
  "\nFailed or non-converged fits: ",
  sum(!fit_summary$converged)
)
