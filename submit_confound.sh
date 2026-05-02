#!/bin/bash
# Confound check: constant-omega sweep at omega=0.125
# This matches the peak switching rate of the diversity-tension kernel at alpha=0.5 (alpha/4 = 0.125).
# Compare the output to sweep_kernel_alpha_0.5.npz to test whether the kernel's
# functional form matters or whether amplitude alone explains the results.
#
# Usage:
#   sbatch submit_confound.sh
#
# Result lands in Files/vacc/sweep_omega_0.1250.npz

#SBATCH --job-name=omega_0125
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=logs/omega_0125.out
#SBATCH --error=logs/omega_0125.err

mkdir -p logs Files/vacc

echo "omega=0.125 sweep starting on $(hostname) at $(date)"

python vacc_sweep.py \
    --omega 0.125 \
    --workers $SLURM_CPUS_PER_TASK

echo "omega=0.125 sweep finished at $(date)"
