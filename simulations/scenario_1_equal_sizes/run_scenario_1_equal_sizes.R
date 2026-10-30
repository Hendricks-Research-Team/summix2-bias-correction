options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages({
  library(dplyr)
  library(Summix)
})

# ============================================================
# Scenario 1: all five reference groups have the same size N
# ============================================================

# -----------------------------
# Configuration
# -----------------------------
sample_sizes <- c(
  5L, 10L, 25L, 35L, 45L, 50L, 75L,
  100L, 150L, 250L, 500L, 750L, 1000L
)

n_reps <- as.integer(Sys.getenv("REPS", "1000"))
ntot <- as.integer(Sys.getenv("NTOT", "10000"))

input_file <- Sys.getenv(
  "REFERENCE_INPUT",
  "/home/math/tifours/simulations/small_chr22_final.gz"
)

base_out_dir <- Sys.getenv(
  "OUTPUT_DIR",
  "/home/math/tifours/summix_sim_outputs_100k/scenario_1_equal_sizes"
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
# Study group and truth vector
# -----------------------------
group <- Sys.getenv("GROUP", "Latinx")

valid_groups <- c("Latinx", "AfricanAmerican", "Peruvian")

if (!group %in% valid_groups) {
  stop(
    "Invalid GROUP: ", group,
    ". Valid groups are: ",
    paste(valid_groups, collapse = ", ")
  )
}

truth_profiles <- list(
  Latinx = c(
    AFR = 0.044,
    EAS = 0.049,
    EUR = 0.438,
    IAM = 0.469,
    SAS = 0.000
  ),
  AfricanAmerican = c(
    AFR = 0.852,
    EAS = 0.000,
    EUR = 0.148,
    IAM = 0.000,
    SAS = 0.000
  ),
  Peruvian = c(
    AFR = 0.027,
    EAS = 0.009,
    EUR = 0.196,
    IAM = 0.768,
    SAS = 0.000
  )
)

ancestries <- c("AFR", "EAS", "EUR", "IAM", "SAS")
truth_named <- truth_profiles[[group]][ancestries]

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

# -----------------------------
# Reference sample sizes
# -----------------------------
N_per_ancestry <- setNames(
  rep(N, length(ancestries)),
  ancestries
)

# -----------------------------
# Deterministic simulation seed
# -----------------------------
scenario_seed_offset <- 100000L

group_seed_offset <- switch(
  group,
  Latinx = 10000L,
  AfricanAmerican = 20000L,
  Peruvian = 30000L
)

simulation_seed <- scenario_seed_offset +
  group_seed_offset +
  task_id

# -----------------------------
# Output location
# -----------------------------
out_dir <- file.path(base_out_dir, group)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message(
  "Scenario: scenario_1_equal_sizes",
  "\nGroup: ", group,
  "\nN: ", N,
  "\nReplicates: ", n_reps,
  "\nObserved population size: ", ntot,
  "\nSimulation seed: ", simulation_seed,
  "\nOutput directory: ", out_dir
)

# -----------------------------
# Prepare reference data
# -----------------------------
refdat <- prepare_reference_data(input_file)

if (nrow(refdat) == 0L) {
  stop("No variants remained after reference filtering.")
}

message(
  "Using all ",
  nrow(refdat),
  " filtered chromosome 22 variants."
)

# -----------------------------
# Run simulation
# -----------------------------
simulation_results <- run_simulation(
  N_per_ancestry = N_per_ancestry,
  n_reps = n_reps,
  ntot = ntot,
  truth_named = truth_named,
  refdat = refdat,
  ancestries = ancestries,
  scenario = "scenario_1_equal_sizes",
  group = group,
  simulation_seed = simulation_seed
)

# -----------------------------
# Save results
# -----------------------------
out_file <- file.path(
  out_dir,
  sprintf(
    "summix_scenario1_group-%s_N-%04d.rds",
    group,
    N
  )
)

saveRDS(
  simulation_results,
  out_file,
  compress = TRUE
)

message(
  "Saved ",
  nrow(simulation_results),
  " replicates to: ",
  out_file
)
