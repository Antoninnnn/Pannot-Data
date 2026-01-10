#!/bin/bash
#SBATCH -J pannot_fill
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 5-00:00:00
#SBATCH --array=0-15
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

# Skip array task 8
if [[ "${SLURM_ARRAY_TASK_ID}" -eq 8 ]]; then
  echo "Skipping train_condensed_008.txt (array task 8)"
  exit 0
fi

ml WebProxy
ml Anaconda3
source activate base

# zero-pad array index to match train_condensed_000.txt naming
IDX=$(printf "%03d" "${SLURM_ARRAY_TASK_ID}")
ID_FILE="data/SwissProt/train_condensed_chunks/train_condensed_${IDX}.txt"

test -f "${ID_FILE}" || { echo "Missing ${ID_FILE}"; ls -lh data/SwissProt/train_condensed_chunks | head; exit 1; }

python -m data.fill_in_templates_parallel \
  --id_list_path "${ID_FILE}" \
  --log_every 50
