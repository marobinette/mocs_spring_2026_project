#!/bin/bash
# SLURM submission script — kernel sweep + baseline, one array job.
#
# Submits a 12-task array job: tasks 0–10 run the chosen kernel at each alpha
# value, task 11 runs the baseline (constant ω=5).  Results land in Files/vacc/
#
# Usage:
#   sbatch submit_sweep.sh
#   NETWORK=Synthetic_poisson_k5 sbatch submit_sweep.sh
#   KERNEL=tension_shifted sbatch submit_sweep.sh
#   KERNEL=tension_shifted NETWORK=Synthetic_poisson_k2 sbatch submit_sweep.sh
#   KERNEL=mixture NETWORK=Thiers13 sbatch submit_sweep.sh  
#   choices=["Thiers13", "Synthetic_poisson_k5", "Synthetic_poisson_k3", "LyonSchool", "Synthetic_poisson_k2"],

NETWORK="${NETWORK:-Thiers13}"
KERNEL="${KERNEL:-diversity_tension}"   # diversity_tension | tension_shifted

echo "Submitting sweep: kernel=${KERNEL}  network=${NETWORK}"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${KERNEL}_${NETWORK}
#SBATCH --array=0-11
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --output=logs/${KERNEL}_${NETWORK}_%a.out
#SBATCH --error=logs/${KERNEL}_${NETWORK}_%a.err

mkdir -p logs Files/vacc

echo "Task \$SLURM_ARRAY_TASK_ID starting on \$(hostname) at \$(date)"

python vacc_sweep.py \\
    --task-index \$SLURM_ARRAY_TASK_ID \\
    --network ${NETWORK} \\
    --kernel ${KERNEL} \\
    --workers \$SLURM_CPUS_PER_TASK

echo "Task \$SLURM_ARRAY_TASK_ID finished at \$(date)"
EOF
