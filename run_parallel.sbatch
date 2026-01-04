#!/bin/bash
#SBATCH -J pannot_fill
#SBATCH -p cpu              # adjust if needed
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4                # UniProt parsing is CPU-bound
#SBATCH --mem=8G
#SBATCH -t 06:00:00
#SBATCH --array=1-40
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err


set -euo pipefail

ml WebProxy
ml Anaconda3
source activate base

cd /scratch/user/yining_yang/TAMU/PhD/Pannot/Pannot-Data
mkdir -p logs

ID_FILE="data/train_ids/train_id_${SLURM_ARRAY_TASK_ID}.txt"

python data/fill_in_templates_parallel.py \
  --id_list_path "${ID_FILE}"
