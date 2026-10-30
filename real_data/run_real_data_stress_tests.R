# ============================================================
# Run real-data Summix reference-size stress tests
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(Summix)
})

# ============================================================
# Configuration
# ============================================================

observed_file <- Sys.getenv(
  "OBSERVED_FILE",
  "small_chr22_final.gz"
)

reference_dir <- Sys.getenv(
  "REFERENCE_DIR",
  "/users/tifours/Downloads"
)

output_dir <- Sys.getenv(
  "OUTPUT_DIR",
  "../results/real_data"
)

observed_col <- "AFafr"

reference_groups <- c(
  "AFR",
  "EUR",
  "IAM",
  "SAS",
  "EAS"
)

reference_size_grid <- c(
  5,
  15,
  25,
  35,
  45,
  75,
  100,
  150,
  250,
  500
)

iam_max_n <- 45L

variant_keys <- c(
  "CHROM",
  "POS",
  "REF",
  "ALT"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(observed_file)) {
  stop(
    "Observed allele-frequency file does not exist: ",
    observed_file
  )
}

# ============================================================
# Load and filter observed chromosome 22 data
# ============================================================

chr22 <- read.table(
  observed_file,
  header = TRUE
)

required_observed_columns <- c(
  variant_keys,
  "AF_AFR",
  "AF_EAS",
  "AF_EUR",
  "AF_IAM",
  "AF_SAS",
  observed_col
)

missing_observed_columns <- setdiff(
  required_observed_columns,
  names(chr22)
)

if (length(missing_observed_columns) > 0L) {
  stop(
    "Observed data are missing required columns: ",
    paste(
      missing_observed_columns,
      collapse = ", "
    )
  )
}

chr22 <- chr22[
  (
    chr22$AF_AFR > 0.01 &
      chr22$AF_AFR < 0.99
  ) |
    (
      chr22$AF_IAM > 0.01 &
        chr22$AF_IAM < 0.99
    ) |
    (
      chr22$AF_EAS > 0.01 &
        chr22$AF_EAS < 0.99
    ) |
    (
      chr22$AF_EUR > 0.01 &
        chr22$AF_EUR < 0.99
    ) |
    (
      chr22$AF_SAS > 0.01 &
        chr22$AF_SAS < 0.99
    ),
]

chr22 <- na.omit(
  chr22
)

observed_data <- chr22 %>%
  select(
    all_of(variant_keys),
    all_of(observed_col)
  )

# ============================================================
# Reference-panel file helpers
# ============================================================

reference_file <- function(
    reference_group,
    sample_size
) {

  extension <- if (
    reference_group == "IAM"
  ) {
    ".tsv"
  } else {
    ".tsv.gz"
  }

  file.path(
    reference_dir,
    paste0(
      reference_group,
      "_data"
    ),
    paste0(
      "hgdp_1kg_subset_by_filtered_",
      reference_group,
      "_",
      sample_size,
      "_1_AC_AN_AF",
      extension
    )
  )
}

read_reference_af <- function(
    reference_group,
    sample_size
) {

  input_file <- reference_file(
    reference_group,
    sample_size
  )

  if (!file.exists(input_file)) {
    stop(
      "Reference file does not exist: ",
      input_file
    )
  }

  reference_data <- read_tsv(
    input_file,
    col_names = c(
      "CHROM",
      "POS",
      "REF",
      "ALT",
      "AC",
      "AN",
      "AF"
    ),
    show_col_types = FALSE
  ) %>%
    transmute(
      CHROM,
      POS,
      REF,
      ALT,
      AF = as.numeric(AF)
    )

  names(reference_data)[
    names(reference_data) == "AF"
  ] <- reference_group

  reference_data
}

# ============================================================
# Build one reference AF panel and run Summix
# ============================================================

run_real_data_summix <- function(
    reference_sizes
) {

  if (
    is.null(names(reference_sizes)) ||
      !all(
        reference_groups %in%
          names(reference_sizes)
      )
  ) {
    stop(
      "reference_sizes must be a named vector containing ",
      paste(
        reference_groups,
        collapse = ", "
      ),
      "."
    )
  }

  reference_sizes <- reference_sizes[
    reference_groups
  ]

  if (
    anyNA(reference_sizes) ||
      any(reference_sizes <= 0)
  ) {
    stop(
      "All reference sample sizes must be positive."
    )
  }

  AFR <- read_reference_af(
    "AFR",
    reference_sizes[["AFR"]]
  )

  EUR <- read_reference_af(
    "EUR",
    reference_sizes[["EUR"]]
  )

  IAM <- read_reference_af(
    "IAM",
    reference_sizes[["IAM"]]
  )

  SAS <- read_reference_af(
    "SAS",
    reference_sizes[["SAS"]]
  )

  EAS <- read_reference_af(
    "EAS",
    reference_sizes[["EAS"]]
  )

  panel <- AFR %>%
    inner_join(
      EUR,
      by = variant_keys
    ) %>%
    inner_join(
      IAM,
      by = variant_keys
    ) %>%
    inner_join(
      SAS,
      by = variant_keys
    ) %>%
    inner_join(
      EAS,
      by = variant_keys
    )

  panel <- panel %>%
    inner_join(
      observed_data,
      by = variant_keys
    )

  # Preserve original handling of missing reference AFs.
  panel[
    reference_groups
  ][
    is.na(
      panel[
        reference_groups
      ]
    )
  ] <- 0

  reference_matrix <- as.matrix(
    panel[
      ,
      reference_groups
    ]
  )

  panel <- panel[
    rowSums(
      reference_matrix
    ) > 0,
  ]

  fit <- Summix::summix(
    data = panel,
    reference = reference_groups,
    observed = observed_col,
    goodness.of.fit = TRUE,
    override_removeSmallRef = TRUE
  )

  output <- as.data.frame(
    fit
  )

  output$N_AFR <- as.integer(
    reference_sizes[["AFR"]]
  )
  output$N_EUR <- as.integer(
    reference_sizes[["EUR"]]
  )
  output$N_IAM <- as.integer(
    reference_sizes[["IAM"]]
  )
  output$N_SAS <- as.integer(
    reference_sizes[["SAS"]]
  )
  output$N_EAS <- as.integer(
    reference_sizes[["EAS"]]
  )

  output$filtered <- nrow(
    panel
  )

  output
}

# ============================================================
# Stress test A
# Jointly limited reference groups
#
# AFR/EAS/EUR/SAS use N.
# IAM uses min(N, 45).
# ============================================================

joint_results <- vector(
  "list",
  length(
    reference_size_grid
  )
)

for (
  i in seq_along(
    reference_size_grid
  )
) {

  N <- reference_size_grid[i]

  N_IAM <- min(
    N,
    iam_max_n
  )

  reference_sizes <- c(
    AFR = N,
    EUR = N,
    IAM = N_IAM,
    SAS = N,
    EAS = N
  )

  output <- run_real_data_summix(
    reference_sizes
  )

  output$StressTest <-
    "jointly_limited"

  output$N <- N
  output$IAM_N <- N_IAM

  joint_results[[i]] <- output

  message(
    "Joint stress test complete: N = ",
    N,
    ", IAM N = ",
    N_IAM,
    ", filtered variants = ",
    output$filtered[[1]]
  )
}

real_data_joint_limited <- bind_rows(
  joint_results
)

# ============================================================
# Stress test B
# AFR varies while other reference groups remain fixed
#
# EUR/EAS/SAS = 500
# IAM = 45
# ============================================================

vary_afr_results <- vector(
  "list",
  length(
    reference_size_grid
  )
)

for (
  i in seq_along(
    reference_size_grid
  )
) {

  N_AFR <- reference_size_grid[i]

  reference_sizes <- c(
    AFR = N_AFR,
    EUR = 500,
    IAM = iam_max_n,
    SAS = 500,
    EAS = 500
  )

  output <- run_real_data_summix(
    reference_sizes
  )

  output$StressTest <-
    "vary_AFR"

  output$AFR_N <- N_AFR
  output$EUR_N <- 500L
  output$IAM_N <- iam_max_n
  output$SAS_N <- 500L
  output$EAS_N <- 500L

  vary_afr_results[[i]] <- output

  message(
    "AFR stress test complete: AFR N = ",
    N_AFR,
    ", filtered variants = ",
    output$filtered[[1]]
  )
}

real_data_vary_afr <- bind_rows(
  vary_afr_results
)

# ============================================================
# Validate shared large-reference benchmark
#
# AFR/EAS/EUR/SAS = 500
# IAM = 45
# ============================================================

joint_benchmark <- real_data_joint_limited %>%
  filter(
    N == 500
  )

afr_benchmark <- real_data_vary_afr %>%
  filter(
    AFR_N == 500
  )

if (
  nrow(joint_benchmark) != 1L ||
    nrow(afr_benchmark) != 1L
) {
  stop(
    "Expected exactly one benchmark row in each stress test."
  )
}

benchmark_groups <- c(
  "AFR",
  "EUR",
  "IAM",
  "SAS",
  "EAS"
)

benchmark_difference <- max(
  abs(
    as.numeric(
      joint_benchmark[
        ,
        benchmark_groups,
        drop = TRUE
      ]
    ) -
      as.numeric(
        afr_benchmark[
          ,
          benchmark_groups,
          drop = TRUE
        ]
      )
  )
)

if (
  benchmark_difference >
    1e-10
) {
  warning(
    "The two stress tests produced different estimates ",
    "for the shared 500/500/500/45/500 benchmark."
  )
}

# ============================================================
# Save uncorrected real-data stress-test results
# ============================================================

joint_output_file <- file.path(
  output_dir,
  "real_data_jointly_limited_uncorrected.rds"
)

vary_afr_output_file <- file.path(
  output_dir,
  "real_data_vary_AFR_uncorrected.rds"
)

saveRDS(
  real_data_joint_limited,
  joint_output_file,
  compress = TRUE
)

saveRDS(
  real_data_vary_afr,
  vary_afr_output_file,
  compress = TRUE
)

message(
  "Saved jointly limited results to: ",
  normalizePath(
    joint_output_file,
    mustWork = FALSE
  ),
  "\nSaved AFR-varying results to: ",
  normalizePath(
    vary_afr_output_file,
    mustWork = FALSE
  )
)
