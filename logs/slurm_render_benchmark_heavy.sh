#!/bin/bash
#SBATCH --job-name=rajive_bench_heavy
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=36:00:00
#SBATCH --output=logs/slurm_render_benchmark_heavy_%j.log
#SBATCH --error=logs/slurm_render_benchmark_heavy_%j.log

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

JOB_LIB="${SLURM_TMPDIR:-/tmp}/rajiveplus_benchmark_heavy_lib_${SLURM_JOB_ID}"
mkdir -p "${JOB_LIB}"

RMD="inst/benchmarks/benchmarking_heavy.Rmd"

echo "====== HEAVY BENCHMARK RENDER START ======"
echo "Timestamp: $(date)"
echo "Vignette: ${RMD}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Output dir: logs/vignette_renders"
echo "Job library: ${JOB_LIB}"

echo ""
echo "Installing current source tree into job-local library"
conda run -n R4_51 R CMD INSTALL . --library="${JOB_LIB}" --no-docs --no-demo --no-test-load 2>&1
INSTALL_EXIT=$?
echo "Install exit code: ${INSTALL_EXIT}"
if [ $INSTALL_EXIT -ne 0 ]; then
  echo "SOURCE INSTALL FAILED"
  exit $INSTALL_EXIT
fi

echo ""
echo "Rendering ${RMD}"
conda run -n R4_51 R --no-save -q -e ".libPaths(c('${JOB_LIB}', .libPaths())); rmarkdown::render('${RMD}', output_format = 'rmarkdown::html_document', output_file = 'benchmarking_heavy.html', output_dir = 'logs/vignette_renders', clean = TRUE)" 2>&1
RENDER_EXIT=$?

echo "Render exit code: ${RENDER_EXIT}"
echo "Timestamp: $(date)"
if [ $RENDER_EXIT -ne 0 ]; then
  echo "VIGNETTE RENDER FAILED: ${RMD}"
  exit $RENDER_EXIT
fi

echo "VIGNETTE RENDER PASSED: ${RMD}"
exit 0
