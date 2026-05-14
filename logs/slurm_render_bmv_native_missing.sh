#!/bin/bash
#SBATCH --job-name=rajive_bmv_native_missing
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm_render_bmv_native_missing_%j.log
#SBATCH --error=logs/slurm_render_bmv_native_missing_%j.log

set -euo pipefail

cd /exports/para-lipg-hpc/mdmanurung/RaJIVEutils

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

if ! command -v conda >/dev/null 2>&1; then
  echo "conda command not found on compute node"
  exit 127
fi

mkdir -p logs/vignette_renders

JOB_LIB="${SLURM_TMPDIR:-/tmp}/rajiveplus_bmv_native_missing_lib_${SLURM_JOB_ID:-manual}"
mkdir -p "${JOB_LIB}"

RMD="inst/analyses/bmv_native_missing_union.Rmd"

echo "====== BMV NATIVE MISSING RENDER START ======"
echo "Timestamp: $(date)"
echo "Artifact: ${RMD}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-8}"
echo "Output dir: logs/vignette_renders"
echo "Job library: ${JOB_LIB}"
echo "RUN_FULL_BMV_NATIVE=${RUN_FULL_BMV_NATIVE:-0}"

echo ""
echo "Installing current source tree into job-local library"
set +e
conda run -n R4_51 R CMD INSTALL . --library="${JOB_LIB}" --no-docs --no-demo --no-test-load 2>&1
INSTALL_EXIT=$?
set -e
echo "Install exit code: ${INSTALL_EXIT}"
if [ "${INSTALL_EXIT}" -ne 0 ]; then
  echo "SOURCE INSTALL FAILED"
  exit "${INSTALL_EXIT}"
fi

echo ""
echo "Rendering ${RMD}"
set +e
conda run -n R4_51 R --no-save -q -e ".libPaths(c('${JOB_LIB}', .libPaths())); Sys.setenv(RAJIVEPLUS_REPO_ROOT = getwd()); rmarkdown::render('${RMD}', output_format = 'rmarkdown::html_document', output_file = 'bmv_native_missing_union.html', output_dir = 'logs/vignette_renders', clean = TRUE)" 2>&1
RENDER_EXIT=$?
set -e

echo "Render exit code: ${RENDER_EXIT}"
echo "Timestamp: $(date)"
if [ "${RENDER_EXIT}" -ne 0 ]; then
  echo "BMV NATIVE MISSING RENDER FAILED: ${RMD}"
  exit "${RENDER_EXIT}"
fi

echo "BMV NATIVE MISSING RENDER PASSED: ${RMD}"
exit 0
