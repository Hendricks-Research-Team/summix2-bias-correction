options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages({
  library(dplyr)
  library(Summix)
})

# ============================================================
# Dirichlet simulations for power-law parameter estimation
# ============================================================

# -----------------------------
# Configuration
# -----------------------------
sample_sizes <- 1:1000

n_reps <- as.integer(Sys.getenv("REPS", "1000"))
ntot <- as.integer(Sys.getenv("NTOT", "10000"))

input_file <- Sys.getenv(
  "REFERENCE_INPUT",
  "/home/math/tifours/simulations/small_chr22_final.gz"
)

proportions_file <- Sys.getenv(
  "PROPORTIONS_FILE",
  "/home/math/tifours/summix2-bias-correction/simulations/power_law_parameter_estimation/dirichlet_100_mixing_proportions.csv"
)

base_out_dir <- Sys.getenv(
  "OUTPUT_DIR",
  "/home/math/tifours/summix_sim_outputs_100k/power_law_parameter_estimation"
)

if (is.na(n_reps) || n_reps < 1L) {
  stop("REPS must be a positive integer.")
}

if (is.na(ntot) || ntot < 1L) {
  stop("NTOT must be a positive integer.")
}

if (!file.exists(input_file)) {
  stop("Reference input file does not exist: ", input_file)
}

if (!file.exists(proportions_file)) {
  stop("Mixing-proportion input file does not exist: ", proportions_file)
}

# -----------------------------
# Shared functions
# -----------------------------
command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)

if (length(file_arg) != 1L) {
  stop("Unable to determine the current script path.")
}

script_dir <- dirname(
  normalizePath(
    sub("^--file=", "", file_arg),
    mustWork = TRUE
  )
)

source(
  file.path(
    script_dir,
    "..",
    "functions",
    "prepare_reference_data.R"
  )
)

source(
  file.path(
    script_dir,
    "..",
    "functions",
    "simulation_functions.R"
  )
)

# -----------------------------
# Select N from Slurm array
# -----------------------------
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

if (
  is.na(task_id) ||
  task_id < 1L ||
  task_id > length(sample_sizes)
) {
  stop(
    "Invalid SLURM_ARRAY_TASK_ID: ", task_id,
    ". It must be between 1 and ", length(sample_sizes), "."
  )
}

N <- sample_sizes[task_id]

ancestries <- c("AFR", "EAS", "EUR", "IAM", "SAS")

N_per_ancestry <- setNames(
  rep(N, length(ancestries)),
  ancestries
)

# -----------------------------
# Read mixing proportions
# -----------------------------
mixing_proportions <- read.csv(
  proportions_file,
  check.names = FALSE
)

required_columns <- c(
  "proportion_id",
  ancestries
)

missing_columns <- setdiff(
  required_columns,
  names(mixing_proportions)
)

if (length(missing_columns) > 0L) {
  stop(
    "Mixing-proportion file is missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyDuplicated(mixing_proportions$proportion_id)) {
  stop("The proportion_id column contains duplicates.")
}

row_sums <- rowSums(
  mixing_proportions[, ancestries, drop = FALSE]
)

if (any(abs(row_sums - 1) > 1e-10)) {
  stop("One or more mixing-proportion rows do not sum to one.")
}

# -----------------------------
# Prepare reference data
# -----------------------------
refdat <- prepare_reference_data(input_file)

if (nrow(refdat) == 0L) {
  stop("No variants remained after reference filtering.")
}

# -----------------------------
# Output location
# -----------------------------
out_dir <- file.path(
  base_out_dir,
  sprintf("N_%04d", N)
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message(
  "Analysis: power_law_parameter_estimation",
  "\nN: ", N,
  "\nMixing-proportion combinations: ", nrow(mixing_proportions),
  "\nReplicates per combination: ", n_reps,
  "\nObserved population size: ", ntot,
  "\nVariants: ", nrow(refdat),
  "\nOutput directory: ", out_dir
)

# -----------------------------
# Run simulations
# -----------------------------
for (row_index in seq_len(nrow(mixing_proportions))) {
  proportion_id <- mixing_proportions$proportion_id[row_index]

  truth_named <- unlist(
    mixing_proportions[
      row_index,
      ancestries,
      drop = FALSE
    ],
    use.names = TRUE
  )

  # Convert Dirichlet proportions to integer-compatible population counts
  expected_counts <- ntot * truth_named
  ancestry_counts <- floor(expected_counts)

  remainder <- ntot - sum(ancestry_counts)

  if (remainder > 0L) {
    fractional_parts <- expected_counts - ancestry_counts

    add_to <- order(
      fractional_parts,
      decreasing = TRUE
    )[seq_len(remainder)]

    ancestry_counts[add_to] <-
      ancestry_counts[add_to] + 1L
  }

  # Use realizable proportions as the simulation truth
  truth_named <- ancestry_counts / ntot
  names(truth_named) <- ancestries

  simulation_seed <- 400000L +
    task_id * 1000L +
    proportion_id

  simulation_results <- run_simulation(
    N_per_ancestry = N_per_ancestry,
    n_reps = n_reps,
    ntot = ntot,
    truth_named = truth_named,
    refdat = refdat,
    ancestries = ancestries,
    scenario = "power_law_parameter_estimation",
    group = paste0("Dirichlet_", proportion_id),
    simulation_seed = simulation_seed
  )

  simulation_results$ProportionID <- proportion_id

  out_file <- file.path(
    out_dir,
    sprintf(
      "summix_power_law_N-%04d_proportion-%03d.rds",
      N,
      proportion_id
    )
  )

  saveRDS(
    simulation_results,
    out_file,
    compress = TRUE
  )

  message(
    "Saved proportion ",
    proportion_id,
    " of ",
    nrow(mixing_proportions),
    ": ",
    out_file
  )
}
