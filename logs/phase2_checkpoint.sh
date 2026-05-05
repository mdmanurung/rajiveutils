#!/bin/bash
set -euo pipefail

# Resolve project root relative to this script so renames/moves do not break checkpoints.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

INTERVAL_SECONDS=3600
LOOP_MODE=0
TRACK_REGEX='rajive_bench2|rajive_jack2|rajive_verify'
OUTFILE='logs/phase2_checkpoints.log'

required_files=(
  bench_results.rds
  bench_results_bench.rds
  peakram_single.rds
  scaling_results.rds
  parallel_results.rds
  jackstraw_time_vs_n.rds
  jackstraw_time_vs_p.rds
  jackstraw_time_vs_nnull.rds
  jackstraw_ram_vs_n.rds
  jackstraw_ram_vs_p.rds
  cll_preprocessed.rds
  cll_svd_list.rds
  cll_rajive_results.rds
  cll_jackstraw_results.rds
)

usage() {
  cat <<EOF
Usage: $0 [--loop] [--interval-seconds N] [--outfile PATH]

Options:
  --loop                 Run checkpoints repeatedly until tracked jobs finish.
  --interval-seconds N  Seconds between checkpoints in loop mode (default: 3600).
  --outfile PATH        Log file path (default: logs/phase2_checkpoints.log).
  --help                Show this help.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --loop)
      LOOP_MODE=1
      shift
      ;;
    --interval-seconds)
      INTERVAL_SECONDS="$2"
      shift 2
      ;;
    --outfile)
      OUTFILE="$2"
      shift 2
      ;;
    --help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1"
      usage
      exit 1
      ;;
  esac
done

run_checkpoint() {
  local ts
  ts="$(date '+%Y-%m-%d %H:%M:%S %Z')"

  local queue_lines
  queue_lines="$(squeue -h -u "$USER" -o '%A|%j|%T|%M|%l|%R' | grep -E "$TRACK_REGEX" || true)"

  local present_count=0
  local missing_count=0
  local missing_list=()

  for f in "${required_files[@]}"; do
    if [[ -f "vignettes/data/$f" ]]; then
      ((present_count+=1))
    else
      ((missing_count+=1))
      missing_list+=("$f")
    fi
  done

  {
    echo "=== Phase2 Checkpoint | $ts ==="
    if [[ -n "$queue_lines" ]]; then
      echo "Tracked jobs (jobid|name|state|elapsed|timelimit|reason):"
      echo "$queue_lines"
    else
      echo "Tracked jobs: none found in queue"
    fi

    echo "Cache completeness: $present_count/${#required_files[@]} present"
    if [[ ${#missing_list[@]} -gt 0 ]]; then
      echo "Missing cache files:"
      for f in "${missing_list[@]}"; do
        echo "- vignettes/data/$f"
      done
    else
      echo "Missing cache files: none"
    fi
    echo
  } | tee -a "$OUTFILE"

  if [[ -n "$queue_lines" ]]; then
    return 0
  fi
  return 1
}

if [[ "$LOOP_MODE" -eq 1 ]]; then
  while true; do
    if ! run_checkpoint; then
      echo "Tracked jobs finished; stopping loop." | tee -a "$OUTFILE"
      exit 0
    fi
    sleep "$INTERVAL_SECONDS"
  done
else
  run_checkpoint || true
fi
