#!/bin/bash
#SBATCH -J pannot_fill
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 5-00:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -uo pipefail
mkdir -p logs

ml WebProxy
ml Anaconda3
source activate base

# Optional: avoid CPU oversubscription if numpy/torch spawn threads
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

CHUNK_DIR="data/SwissProt/train_condensed_chunks"
CHUNK_PREFIX="train_condensed_"
MAX_PARALLEL=2
TOTAL_FILES=16   # 000..015

FAILED=()
pids=()
labels=()

run_one() {
  local i="$1"
  local idx id_file
  idx=$(printf "%03d" "${i}")
  id_file="${CHUNK_DIR}/${CHUNK_PREFIX}${idx}.txt"

  if [[ ! -f "${id_file}" ]]; then
    echo "Missing ${id_file}, skipping."
    return 0
  fi

  # If you still want to skip 008, keep this; otherwise remove it.
  if (( i == 8 )); then
    echo "Skipping index 8 (${CHUNK_PREFIX}008.txt)"
    return 0
  fi

  echo ">>> Processing ${id_file}"
  python -m data.fill_in_templates_parallel \
    --id_list_path "${id_file}" \
    --log_every 50
}

wait_one_slot() {
  # Wait for the oldest running PID, record failure if needed
  local pid="${pids[0]}"
  local lab="${labels[0]}"

  if ! wait "${pid}"; then
    echo "WARNING: file ${lab} failed (pid=${pid})"
    FAILED+=("${lab}")
  else
    echo "OK: file ${lab} finished"
  fi

  # Pop front
  pids=("${pids[@]:1}")
  labels=("${labels[@]:1}")
}

for i in $(seq 0 $((TOTAL_FILES - 1))); do
  lab=$(printf "%03d" "${i}")

  # If we already have MAX_PARALLEL running, wait for one to finish
  while (( ${#pids[@]} >= MAX_PARALLEL )); do
    wait_one_slot
  done

  # Launch next job
  ( run_one "${i}" ) &
  pids+=("$!")
  labels+=("${lab}")
done

# Wait for remaining jobs
while (( ${#pids[@]} > 0 )); do
  wait_one_slot
done

echo "======================================"
if (( ${#FAILED[@]} > 0 )); then
  echo "Completed with failures in files: ${FAILED[*]}"
  exit 1
else
  echo "All files completed successfully."
fi