#!/bin/bash
#SBATCH --job-name=rajive_bench2
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=36:00:00
#SBATCH --output=logs/slurm_bench2_%j.log
#SBATCH --error=logs/slurm_bench2_%j.log

set -euo pipefail

# Under sbatch, scripts are executed from a spool copy; use submit directory.
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

# Keep BLAS/OpenMP single-threaded to avoid nested over-subscription when
# benchmarks themselves use explicit parallelism.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

R_BIN="/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/R"
R_SCRIPT="/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript"

cleanup() {
  perl -0pi -e 's/run_heavy <- TRUE\s+# set to TRUE once to regenerate \.rds cache files/run_heavy <- FALSE         # set to TRUE once to regenerate .rds cache files/g' vignettes/benchmarking.Rmd || true
}
trap cleanup EXIT

perl -0pi -e 's/run_heavy <- FALSE\s+# set to TRUE once to regenerate \.rds cache files/run_heavy <- TRUE          # set to TRUE once to regenerate .rds cache files/g' vignettes/benchmarking.Rmd

"${R_BIN}" CMD INSTALL .
"${R_SCRIPT}" -e "rmarkdown::render('vignettes/benchmarking.Rmd')"

echo "EXIT CODE: $?"
