#!/bin/bash

condatainer scheduler | head -7
# --debug will keep the sbatch script (by default, it will be removed after submission)

# Slurm test
condatainer create grcm39/salmon/1.10.2/gencodeM6
# Test
# 1. chain dep submit (Slurm)
# 2. output path (current dir)
condatainer run --debug src/run_chain_dep.sh --auto-install
condatainer run --debug src/run_chain_external.sh --auto-install # should fail (external img cannot auto-install)

# LSF test (Cross-scheduler test)
condatainer run --debug src/run_chain_lsf.sh
# PBS test (Cross-scheduler test)
condatainer run --debug src/run_chain_pbs.sh

# GPU jobs test
condatainer create -p src/pytorch -f src/pytorch.def
condatainer run --debug src/gpu_test_slurm.sh

# MPI jobs test
condatainer overlay create --sparse mpi.img
ml av openmpi # openmpi/4.1.5
# Make sure install the same major.minor version of openmpi
condatainer e mpi.img -- mm-install mpi4py openmpi=4.1 -y
condatainer run --debug src/mpi_test_slurm.sh
cat src/logs/output.txt
rm src/logs/output.txt

# Test Passthrough workaround
sbatch src/mpi_test_slurm_passthrough.sh
cat src/logs/output.txt
rm src/logs/output.txt

# Pure slurm test
# salloc --ntasks=4 --time=00:30:00
# ml purge
# ml openmpi mpi4py
# mpiexec python src/mpi_file_handler.py
# condatainer run src/mpi_test_slurm.sh

# run args test
condatainer run --debug src/run_arg.sh test 1 3 'asdfsd asds'
