#!/bin/bash
# SLURM submission script — diversity-tension kernel or baseline.
#
# Submits a 13-task array job (one per alpha value) for the kernel sweep, or
# a single job for the baseline (constant ω=5) sweep.
# Results land in Files/vacc/
#
# Usage:
#   sbatch submit_sweep.sh                                  # kernel, Thiers13
#   NETWORK=Synthetic_poisson_k5 sbatch submit_sweep.sh    # kernel, other network
#   BASELINE=1 sbatch submit_sweep.sh                      # baseline, Thiers13
#   BASELINE=1 NETWORK=Synthetic_poisson_k5 sbatch submit_sweep.sh

NETWORK="${NETWORK:-Thiers13}"
BASELINE="${BASELINE:-0}"

if [ "$BASELINE" = "1" ]; then

echo "Submitting baseline sweep for NETWORK=$NETWORK"

sbatch <<BASELINE_JOB
#!/bin/bash
#SBATCH --job-name=baseline_${NETWORK}
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --output=logs/baseline_${NETWORK}.out
#SBATCH --error=logs/baseline_${NETWORK}.err

mkdir -p logs Files/vacc

echo "Baseline job starting on \$(hostname) at \$(date)"

python vacc_sweep.py \\
    --baseline \\
    --network ${NETWORK} \\
    --workers \$SLURM_CPUS_PER_TASK

echo "Baseline job finished at \$(date)"
BASELINE_JOB

else

echo "Submitting kernel sweep for NETWORK=$NETWORK"

# ALPHA_VALUES = np.logspace(-1, 2, 13)
# ~ [0.10, 0.18, 0.32, 0.56, 1.0, 1.78, 3.16, 5.62, 10.0, 17.8, 31.6, 56.2, 100.0]
#      0     1     2     3    4     5     6     7     8     9    10    11    12

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

fi
