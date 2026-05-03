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
# ALPHA_VALUES = np.logspace(-1, 2, 13)
# ~ [0.10, 0.18, 0.32, 0.56, 1.0, 1.78, 3.16, 5.62, 10.0, 17.8, 31.6, 56.2, 100.0]
#      0     1     2     3    4     5     6     7     8     9    10    11    12
# ---------------------------------------------------------------------------
sbatch <<'KERNEL_JOB'
#!/bin/bash
#SBATCH --job-name=kernel_sweep
#SBATCH --array=0-12
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=2:00:00
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
# Job B: constant-omega baseline (omega=5.0)
# ---------------------------------------------------------------------------
sbatch <<'BASELINE_JOB'
#!/bin/bash
#SBATCH --job-name=baseline_sweep
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=logs/baseline.out
#SBATCH --error=logs/baseline.err

mkdir -p logs Files/vacc

echo "Baseline starting on $(hostname) at $(date)"

python vacc_sweep.py \
    --baseline \
    --workers $SLURM_CPUS_PER_TASK

echo "Baseline finished at $(date)"
BASELINE_JOB

