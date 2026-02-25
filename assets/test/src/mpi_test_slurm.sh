#!/bin/bash
#SBATCH --job-name=mpi_test
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --time=00:10:00
#SBATCH --output=src/logs/slurm_log.txt

#DEP:mpi.img

echo ======In Job Testing========
echo In Job Testing
echo NNODES $NNODES
echo NTASKS $NTASKS
echo NTASKS_PER_NODE $NTASKS_PER_NODE
echo NCPUS $NCPUS
echo NCPUS_PER_TASK $NCPUS_PER_TASK
echo MEM $MEM
echo MEM MB $MEM_MB
echo MEM GB $MEM_GB
echo ============================

# Run under test folder not the src folder.
python src/mpi_file_handler.py
