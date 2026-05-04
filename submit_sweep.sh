#!/bin/bash
# SLURM submission script — diversity-tension kernel, inverted kernel, baseline.
#
# Submits three jobs for the selected network:
#   Job A — diversity-tension kernel array (13 alpha values)
#   Job B — constant-omega baseline
#   Job C — inverted kernel array (13 alpha values)
#
# Usage:
#   NETWORK=Thiers13 sbatch submit_sweep.sh
#   NETWORK=Synthetic_poisson_k5 sbatch submit_sweep.sh
#
# Or edit the default below and run:
#   sbatch submit_sweep.sh
#
# Results land in Files/vacc/

NETWORK="${NETWORK:-Thiers13}"   # override with:  NETWORK=Synthetic_poisson_k5 sbatch submit_sweep.sh

echo "Submitting all three sweeps for NETWORK=$NETWORK"

# ---------------------------------------------------------------------------
# Job A: diversity-tension kernel sweeps (array, one task per alpha)
# ALPHA_VALUES = np.logspace(-1, 2, 13)
# ~ [0.10, 0.18, 0.32, 0.56, 1.0, 1.78, 3.16, 5.62, 10.0, 17.8, 31.6, 56.2, 100.0]
#      0     1     2     3    4     5     6     7     8     9    10    11    12
# ---------------------------------------------------------------------------
sbatch <<KERNEL_JOB
#!/bin/bash
#SBATCH --job-name=kernel_${NETWORK}
#SBATCH --array=0-12
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --output=logs/kernel_${NETWORK}_%a.out
#SBATCH --error=logs/kernel_${NETWORK}_%a.err

mkdir -p logs Files/vacc

echo "Task \$SLURM_ARRAY_TASK_ID starting on \$(hostname) at \$(date)"

python vacc_sweep.py \\
    --alpha-index \$SLURM_ARRAY_TASK_ID \\
    --network ${NETWORK} \\
    --workers \$SLURM_CPUS_PER_TASK

echo "Task \$SLURM_ARRAY_TASK_ID finished at \$(date)"
KERNEL_JOB

# ---------------------------------------------------------------------------
# Job B: constant-omega baseline (omega=5.0)
# ---------------------------------------------------------------------------
sbatch <<BASELINE_JOB
#!/bin/bash
#SBATCH --job-name=baseline_${NETWORK}
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=5:00:00
#SBATCH --output=logs/baseline_${NETWORK}.out
#SBATCH --error=logs/baseline_${NETWORK}.err

mkdir -p logs Files/vacc

echo "Baseline starting on \$(hostname) at \$(date)"

python vacc_sweep.py \\
    --baseline \\
    --network ${NETWORK} \\
    --workers \$SLURM_CPUS_PER_TASK

echo "Baseline finished at \$(date)"
BASELINE_JOB

# ---------------------------------------------------------------------------
# Job C: inverted kernel sweeps (array, one task per alpha)
# Same ALPHA_VALUES as Job A — direct like-for-like comparison.
# ---------------------------------------------------------------------------
sbatch <<INVERTED_JOB
#!/bin/bash
#SBATCH --job-name=inverted_${NETWORK}
#SBATCH --array=0-12
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --output=logs/inverted_${NETWORK}_%a.out
#SBATCH --error=logs/inverted_${NETWORK}_%a.err

mkdir -p logs Files/vacc

echo "Task \$SLURM_ARRAY_TASK_ID starting on \$(hostname) at \$(date)"

python vacc_sweep.py \\
    --inverted-index \$SLURM_ARRAY_TASK_ID \\
    --network ${NETWORK} \\
    --workers \$SLURM_CPUS_PER_TASK

echo "Task \$SLURM_ARRAY_TASK_ID finished at \$(date)"
INVERTED_JOB
