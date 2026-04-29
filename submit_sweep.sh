#!/bin/bash
# SLURM submission script for the diversity-tension kernel parameter sweep.
#
# Submits two jobs:
#   Job A — array of kernel sweeps, one per alpha value (9 jobs)
#   Job B — single baseline (constant omega) sweep
#
# Usage:
#   sbatch submit_sweep.sh
#
# Results land in Files/vacc/

# ---------------------------------------------------------------------------
# Job A: kernel sweeps (array, one task per alpha)
# ALPHA_VALUES = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0]
#                  0    1    2    3    4     5     6     7     8
# ---------------------------------------------------------------------------
sbatch <<'KERNEL_JOB'
#!/bin/bash
#SBATCH --job-name=kernel_sweep
#SBATCH --array=0-8
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=logs/kernel_%a.out
#SBATCH --error=logs/kernel_%a.err

mkdir -p logs Files/vacc

echo "Task $SLURM_ARRAY_TASK_ID starting on $(hostname) at $(date)"

python vacc_sweep.py \
    --alpha-index $SLURM_ARRAY_TASK_ID \
    --workers $SLURM_CPUS_PER_TASK

echo "Task $SLURM_ARRAY_TASK_ID finished at $(date)"
KERNEL_JOB

# ---------------------------------------------------------------------------
# Job B: constant-omega baseline
# ---------------------------------------------------------------------------
sbatch <<'BASELINE_JOB'
#!/bin/bash
#SBATCH --job-name=baseline_sweep
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=logs/baseline.out
#SBATCH --error=logs/baseline.err

mkdir -p logs Files/vacc

echo "Baseline starting on $(hostname) at $(date)"

python vacc_sweep.py \
    --baseline \
    --workers $SLURM_CPUS_PER_TASK

echo "Baseline finished at $(date)"
BASELINE_JOB
