#!/bin/bash
# SLURM submission script — diversity-tension kernel + baseline, one array job.
#
# Submits a 12-task array job: tasks 0–10 run the kernel at each alpha value,
# task 11 runs the baseline (constant ω=5).  Results land in Files/vacc/
#
# Usage:
#   sbatch submit_sweep.sh                                  # Thiers13
#   NETWORK=Synthetic_poisson_k5 sbatch submit_sweep.sh    # other network

NETWORK="${NETWORK:-Thiers13}"

echo "Submitting combined sweep (11 alpha + baseline) for NETWORK=$NETWORK"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=sweep_${NETWORK}
#SBATCH --array=0-11
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --output=logs/sweep_${NETWORK}_%a.out
#SBATCH --error=logs/sweep_${NETWORK}_%a.err

mkdir -p logs Files/vacc

echo "Task \$SLURM_ARRAY_TASK_ID starting on \$(hostname) at \$(date)"

python vacc_sweep.py \\
    --task-index \$SLURM_ARRAY_TASK_ID \\
    --network ${NETWORK} \\
    --workers \$SLURM_CPUS_PER_TASK

echo "Task \$SLURM_ARRAY_TASK_ID finished at \$(date)"
EOF
