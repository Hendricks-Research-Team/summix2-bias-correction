options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages({
  library(gtools)
})

# ============================================================
# Generate Dirichlet mixing-proportion combinations
# ============================================================

# -----------------------------
# Configuration
# -----------------------------
n_combinations <- 100L
dirichlet_seed <- 42L
ancestries <- c("AFR", "EAS", "EUR", "IAM", "SAS")

output_file <- Sys.getenv(
  "OUTPUT_FILE",
  "dirichlet_100_mixing_proportions.csv"
)

# -----------------------------
# Generate mixing proportions
# -----------------------------
set.seed(dirichlet_seed)

mixing_proportions <- gtools::rdirichlet(
  n = n_combinations,
  alpha = rep(1, length(ancestries))
)

colnames(mixing_proportions) <- ancestries

mixing_proportions <- data.frame(
  proportion_id = seq_len(n_combinations),
  mixing_proportions,
  check.names = FALSE
)

# -----------------------------
# Validate output
# -----------------------------
row_sums <- rowSums(
  mixing_proportions[, ancestries, drop = FALSE]
)

if (any(abs(row_sums - 1) > 1e-12)) {
  stop("One or more mixing-proportion rows do not sum to one.")
}

# -----------------------------
# Save results
# -----------------------------
write.csv(
  mixing_proportions,
  output_file,
  row.names = FALSE
)

message(
  "Saved ",
  n_combinations,
  " Dirichlet mixing-proportion combinations to: ",
  normalizePath(output_file, mustWork = FALSE)
)
