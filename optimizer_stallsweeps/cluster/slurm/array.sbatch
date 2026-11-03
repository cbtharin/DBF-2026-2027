#!/bin/bash
# Many searches at once, one per line of a jobs file, each on its own node:
#     sbatch --array=1-4 slurm/array.sbatch slurm/jobs_example.txt
#     sbatch --array=1-8%4 slurm/array.sbatch slurm/jobs_battery.txt   # at most 4 running
# Each line is   name | run_cluster.py arguments   (# starts a comment).
# Results: results/<name>_<array job id>/
#SBATCH --job-name=dbf-array
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
# memory: 32 GB. seff reported 1.7 GB for a full search, but that was one
# process; with 8 GB, two of five searches lost their workers (2026-10-14)
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
cd "$SLURM_SUBMIT_DIR" || exit 1
jobs="${1:-slurm/jobs_example.txt}"
line=$(grep -v -E '^\s*(#|$)' "$jobs" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$line" ] || { echo "no line $SLURM_ARRAY_TASK_ID in $jobs"; exit 1; }
name=$(echo "${line%%|*}" | xargs)
read -r -a args <<< "${line#*|}"
source ./sapelo2_env.sh || exit 1
out="results/${name}_${SLURM_ARRAY_JOB_ID}"
mkdir -p "$out"
echo "job $SLURM_ARRAY_TASK_ID: $name: ${args[*]}"
"$PY" -u run_cluster.py "${args[@]}" --out "$out" 2>&1 | tee "$out/search_log.txt"
