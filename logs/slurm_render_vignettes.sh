#!/bin/bash
#SBATCH --job-name=rajive_vignettes
#SBATCH --partition=all
#SBATCH --array=1-6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm_render_vignettes_%A_%a.log
#SBATCH --error=logs/slurm_render_vignettes_%A_%a.log

cd /exports/para-lipg-hpc/mdmanurung/RaJIVEutils

if ! command -v conda >/dev/null 2>&1; then
  echo "conda command not found on compute node"
  exit 127
fi

mkdir -p logs/vignette_renders

JOB_LIB="${SLURM_TMPDIR:-/tmp}/rajiveplus_vignette_lib_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${JOB_LIB}"

VIGNETTES=(
  "vignettes/benchmarking.Rmd"
  "vignettes/cll_application.Rmd"
  "vignettes/function_gallery.Rmd"
  "vignettes/inference.Rmd"
  "vignettes/jackstraw_scaling.Rmd"
  "vignettes/microbiome_application.Rmd"
)

IDX=$((SLURM_ARRAY_TASK_ID - 1))
RMD="${VIGNETTES[$IDX]}"

if [ -z "$RMD" ]; then
  echo "No vignette mapped for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
  exit 2
fi

echo "====== VIGNETTE RENDER START ======"
echo "Timestamp: $(date)"
echo "Task: ${SLURM_ARRAY_TASK_ID}"
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
conda run -n R4_51 R --no-save -q -e ".libPaths(c('${JOB_LIB}', .libPaths())); rmarkdown::render('${RMD}', output_format = 'rmarkdown::html_vignette', output_dir = 'logs/vignette_renders', clean = TRUE)" 2>&1
RENDER_EXIT=$?

echo "Render exit code: ${RENDER_EXIT}"
echo "Timestamp: $(date)"
if [ $RENDER_EXIT -ne 0 ]; then
  echo "VIGNETTE RENDER FAILED: ${RMD}"
  exit $RENDER_EXIT
fi

echo "VIGNETTE RENDER PASSED: ${RMD}"
exit 0
