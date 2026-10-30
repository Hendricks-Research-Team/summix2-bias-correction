#!/bin/bash
#SBATCH --job-name=summix_powerlaw
#SBATCH --array=1-1000
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=logs/powerlaw_%x_%A_%a.out
#SBATCH --error=logs/powerlaw_%x_%A_%a.err

set -euo pipefail

module purge
module load R/4.4.0
export R_LIBS_USER=/home/tifours/Rlibs

export REPS="${REPS:-1000}"
export NTOT="${NTOT:-10000}"

export REFERENCE_INPUT="${REFERENCE_INPUT:-/home/math/tifours/simulations/small_chr22_final.gz}"

export PROPORTIONS_FILE="${PROPORTIONS_FILE:-/home/math/tifours/summix2-bias-correction/simulations/power_law_parameter_estimation/dirichlet_100_mixing_proportions.csv}"

export OUTPUT_DIR="${OUTPUT_DIR:-/home/math/tifours/summix2-bias-correction/results/simulation/power_law_parameter_estimation}"

Rscript /home/math/tifours/summix2-bias-correction/simulations/power_law_parameter_estimation/run_dirichlet_simulations.R