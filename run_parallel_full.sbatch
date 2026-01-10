#!/bin/bash
#SBATCH -J pannot_fill
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 5-00:00:00
#SBATCH --array=1-41
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail

mkdir -p logs

ml WebProxy
ml Anaconda3
source activate base

ID_FILE="data/SwissProt/train_ids/train_ids_${SLURM_ARRAY_TASK_ID}.txt"
test -f "${ID_FILE}" || { echo "Missing ${ID_FILE}"; exit 1; }

python -m data.fill_in_templates_parallel \
  --id_list_path "${ID_FILE}" \
  --log_every 50
