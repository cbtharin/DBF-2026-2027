#!/bin/bash
# One full-scale search on one node - the exp8 settings (9.9 h on the PC).
#     sbatch slurm/search.sbatch                      # as below
#     sbatch slurm/search.sbatch --seed 7             # extra arguments are passed on
# Results: results/<job name>_<job id>/ (ranked_designs.json, top25/, search_log.txt)
#SBATCH --job-name=dbf-search
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
# memory: 32 GB. seff reported 1.7 GB for a full search, but that was one
# process; with 8 GB, two of five searches lost their workers (2026-10-14)
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=YOUR_ID@uga.edu
cd "$SLURM_SUBMIT_DIR" || exit 1
source ./sapelo2_env.sh || exit 1
out="results/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
mkdir -p "$out"
"$PY" -u run_cluster.py --thorough --top 25 --finalists 32 --motor-counts 1 2 \
    --keep-propulsion 600 --refine 100 --swap-parents 12 \
    --seeds seeds --out "$out" "$@" 2>&1 | tee "$out/search_log.txt"
