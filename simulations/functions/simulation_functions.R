# ============================================================
# Simulation functions
# ============================================================

simulate_finite_reference_af <- function(true_af, N) {
  vapply(
    true_af,
    function(p) {
      counts <- rmultinom(
        n = 1L,
        size = 2L * N,
        prob = c(p, 1 - p)
      )

      counts[1L, 1L] / (2L * N)
    },
    numeric(1)
  )
}


simulate_observed_af <- function(
    refdat,
    truth_named,
    ntot,
    ancestries
) {
  group_sizes <- ntot * truth_named[ancestries]

  if (any(abs(group_sizes - round(group_sizes)) > 1e-8)) {
    stop(
      "ntot multiplied by each mixing proportion must produce integers."
    )
  }

  group_sizes <- setNames(
    as.integer(round(group_sizes)),
    ancestries
  )

  if (sum(group_sizes) != ntot) {
    stop("Ancestry-specific observed sample sizes do not sum to ntot.")
  }

  genotype_counts <- matrix(
    0,
    nrow = nrow(refdat),
    ncol = 3L
  )

  for (ancestry in ancestries) {
    ancestry_n <- group_sizes[[ancestry]]

    if (ancestry_n == 0L) {
      next
    }

    ancestry_counts <- t(
      vapply(
        refdat[[ancestry]],
        function(p) {
          rmultinom(
            n = 1L,
            size = ancestry_n,
            prob = c(
              p^2,
              2 * p * (1 - p),
              (1 - p)^2
            )
          )[, 1L]
        },
        numeric(3)
      )
    )

    genotype_counts <- genotype_counts + ancestry_counts
  }

  (
    2 * genotype_counts[, 1L] +
      genotype_counts[, 2L]
  ) / (2 * ntot)
}


run_simulation <- function(
    N_per_ancestry,
    n_reps,
    ntot,
    truth_named,
    refdat,
    ancestries,
    scenario,
    group,
    simulation_seed,
    varying_ancestry = NA_character_,
    N_var = NA_integer_,
    N_fixed = NA_integer_
) {
  if (!setequal(names(N_per_ancestry), ancestries)) {
    stop(
      "N_per_ancestry must be a named vector containing: ",
      paste(ancestries, collapse = ", ")
    )
  }

  N_per_ancestry <- N_per_ancestry[ancestries]

  if (
    anyNA(N_per_ancestry) ||
    any(N_per_ancestry < 1L) ||
    any(N_per_ancestry != as.integer(N_per_ancestry))
  ) {
    stop("All reference sample sizes must be positive integers.")
  }

  truth_named <- truth_named[ancestries]

  if (anyNA(truth_named)) {
    stop("The truth vector is missing one or more reference groups.")
  }

  if (abs(sum(truth_named) - 1) > 1e-10) {
    stop(
      "The mixing proportions for ", group,
      " do not sum to one. Sum = ", sum(truth_named)
    )
  }

  set.seed(simulation_seed)

  sim_data <- data.frame(
    rsID = refdat$rsID,
    AFR = NA_real_,
    EAS = NA_real_,
    EUR = NA_real_,
    IAM = NA_real_,
    SAS = NA_real_,
    AF = NA_real_
  )

  result_list <- vector("list", n_reps)

  for (replicate_id in seq_len(n_reps)) {
    for (ancestry in ancestries) {
      sim_data[[ancestry]] <- simulate_finite_reference_af(
        true_af = refdat[[ancestry]],
        N = N_per_ancestry[[ancestry]]
      )
    }

    sim_data$AF <- simulate_observed_af(
      refdat = refdat,
      truth_named = truth_named,
      ntot = ntot,
      ancestries = ancestries
    )

    summix_result <- Summix::summix(
      data = sim_data,
      reference = ancestries,
      observed = "AF",
      goodness.of.fit = TRUE,
      override_removeSmallRef = TRUE
    )

    summix_result$Replicate <- replicate_id
    result_list[[replicate_id]] <- summix_result
  }

  estimates <- dplyr::bind_rows(result_list)

  expected_columns <- c(
    "goodness.of.fit",
    "iterations",
    "time",
    "filtered",
    ancestries,
    "Replicate"
  )

  missing_columns <- setdiff(
    expected_columns,
    names(estimates)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Summix output is missing expected columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  simulation_results <- estimates %>%
    dplyr::select(dplyr::all_of(expected_columns)) %>%
    dplyr::mutate(
      Scenario = scenario,
      Group = group,
      VaryingAncestry = varying_ancestry,
      N_var = N_var,
      N_fixed = N_fixed,
      N_AFR = N_per_ancestry[["AFR"]],
      N_EAS = N_per_ancestry[["EAS"]],
      N_EUR = N_per_ancestry[["EUR"]],
      N_IAM = N_per_ancestry[["IAM"]],
      N_SAS = N_per_ancestry[["SAS"]],
      true_AFR = truth_named[["AFR"]],
      true_EAS = truth_named[["EAS"]],
      true_EUR = truth_named[["EUR"]],
      true_IAM = truth_named[["IAM"]],
      true_SAS = truth_named[["SAS"]],
      bias_AFR = AFR - true_AFR,
      bias_EAS = EAS - true_EAS,
      bias_EUR = EUR - true_EUR,
      bias_IAM = IAM - true_IAM,
      bias_SAS = SAS - true_SAS,
      absolute_bias_AFR = abs(bias_AFR),
      absolute_bias_EAS = abs(bias_EAS),
      absolute_bias_EUR = abs(bias_EUR),
      absolute_bias_IAM = abs(bias_IAM),
      absolute_bias_SAS = abs(bias_SAS),
      total_absolute_deviation =
        absolute_bias_AFR +
        absolute_bias_EAS +
        absolute_bias_EUR +
        absolute_bias_IAM +
        absolute_bias_SAS,
      NumberSNPs = nrow(refdat),
      ObservedPopulationSize = ntot,
      SimulationSeed = simulation_seed
    )

  simulation_results
}
