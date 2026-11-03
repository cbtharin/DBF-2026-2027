#!/bin/bash
# One search spread over several nodes with MPI.  Worth it for the screen,
# coarse sweep and swap scoring; the refinement (one Nelder-Mead per
# configuration, ~100 of them) does not get faster past ~100 workers.
# Needs mpi4py: check the module name with   module spider mpi4py
#SBATCH --job-name=dbf-mpi
#SBATCH --partition=batch
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
cd "$SLURM_SUBMIT_DIR" || exit 1
module purge
mod=$(module -t avail mpi4py 2>&1 | grep -E '^mpi4py/' | sort -V | tail -1)
[ -n "$mod" ] || { echo "no mpi4py module (module spider mpi4py)"; exit 1; }
module load "$mod"          # brings its own Python and MPI
out="results/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
mkdir -p "$out"
srun python3 -u -m mpi4py.futures run_cluster.py --mpi \
    --thorough --top 25 --finalists 32 --motor-counts 1 2 \
    --keep-propulsion 600 --refine 100 --swap-parents 12 \
    --seeds seeds --out "$out" "$@" 2>&1 | tee "$out/search_log.txt"
