#!/bin/bash
#SBATCH --job-name=mpi_passthrough
#SBATCH --nodes=3
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --time=00:10:00
#SBATCH --output=src/logs/slurm_passthrough_log.txt

#DEP:mpi.img

if [ -z "$IN_CONDATAINER" ] && command -v condatainer >/dev/null 2>&1; then
    if [ -n "$SLURM_JOB_ID" ]; then
        FULL_COMMAND=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}' | head -n 1)
        ORIGINAL_SCRIPT_PATH=$(echo "$FULL_COMMAND" | awk '{print $1}')
    else
        ORIGINAL_SCRIPT_PATH=$(realpath "$0")
    fi
    module purge && module load openmpi/4.1.5 && \
        mpiexec condatainer run "$ORIGINAL_SCRIPT_PATH"
    exit $?
fi

# Run under test folder not the src folder.
python src/mpi_file_handler.py
