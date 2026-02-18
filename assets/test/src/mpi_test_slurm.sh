#!/bin/bash
#SBATCH --job-name=python_mpi_file
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --time=00:10:00
#SBATCH --output=src/logs/slurm_log.txt

#DEP:mpi.img

# Run under test folder not the src folder.
python src/mpi_file_handler.py
