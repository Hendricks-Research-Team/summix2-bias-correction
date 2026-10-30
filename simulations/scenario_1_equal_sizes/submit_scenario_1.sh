#!/bin/bash
#SBATCH --job-name=summix_sc1
#SBATCH --array=1-13
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=logs/scenario1_%x_%A_%a.out
#SBATCH --error=logs/scenario1_%x_%A_%a.err

set -euo pipefail

module purge
module load R/4.4.0
export R_LIBS_USER=/home/tifours/Rlibs

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

: "${GROUP:?Set GROUP to Latinx, AfricanAmerican, or Peruvian}"

export GROUP
export REPS="${REPS:-1000}"
export NTOT="${NTOT:-10000}"
export REFERENCE_INPUT="${REFERENCE_INPUT:-/home/math/tifours/simulations/small_chr22_final.gz}"
export OUTPUT_DIR="${OUTPUT_DIR:-/home/math/tifours/summix2-bias-correction/results/simulation}"

Rscript /home/math/tifours/summix2-bias-correction/simulations/scenario_1_equal_sizes/run_scenario_1_equal_sizes.R
