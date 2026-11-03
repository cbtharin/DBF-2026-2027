#!/bin/bash
# Same answers as the PC?  How fast is this node?  ~5 minutes.
# Ask for the core count you will search with: the timing loads them all.
#SBATCH --job-name=dbf-check
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
cd "$SLURM_SUBMIT_DIR" || exit 1
source ./sapelo2_env.sh || exit 1
"$PY" -u check_cluster.py
