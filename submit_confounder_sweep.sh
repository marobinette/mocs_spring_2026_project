#!/bin/bash
# Confound check: constant-omega sweep at omega=0.125.
# Matches the peak of the diversity-tension kernel at alpha=0.5 (alpha/4 = 0.125),
# the most bistable kernel case. Uses the paper's original model (S and I as scalars).
#
# Usage:
#   sbatch submit_confounder_sweep.sh
#
# Results land in Files/vacc/sweep_omega_0.1250.npz

#SBATCH --job-name=confound_0125
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=logs/confound_0125.out
#SBATCH --error=logs/confound_0125.err

mkdir -p logs Files/vacc

echo "Confound sweep omega=0.125 starting on $(hostname) at $(date)"

python vacc_sweep.py \
    --omega 0.125 \
    --workers $SLURM_CPUS_PER_TASK

echo "Confound sweep finished at $(date)"
