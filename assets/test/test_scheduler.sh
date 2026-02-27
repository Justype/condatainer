#!/bin/bash

condatainer scheduler | head -7
# --debug will keep the sbatch script (by default, it will be removed after submission)

# Slurm test
condatainer create grcm39/salmon/1.10.2/gencodeM6
# Test
# 1. chain dep submit (Slurm)
# 2. output path (current dir)
condatainer run --dry-run src/run_chain_slurm.sh
condatainer run --debug --auto-install src/run_chain_dep.sh
condatainer run --debug --auto-install src/run_chain_external.sh # should fail (external img cannot auto-install)

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
condatainer run --dry-run src/mpi_test_slurm.sh
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
condatainer run --dry-run src/run_arg.sh test 1 3 'asdfsd asds'
condatainer run --debug src/run_arg.sh test 1 3 'asdfsd asds'

condatainer run --dry-run -o src/logs/run_arg_test_12.out src/run_arg.sh test 12 'another arg'
condatainer run --debug -o src/logs/run_arg_test_13.out src/run_arg.sh test 13 'another arg'

# run dep test
JOB=$(condatainer run --debug src/run_long.sh)
condatainer run --debug --afterok "$JOB" src/run_dep.sh

JOB=
condatainer run --afterok "$JOB" src/run_dep.sh # should fail, ID is not valid
JOB=abc123
condatainer run --afterok "$JOB" src/run_dep.sh # should fail, ID is not valid
JOB=12345.67890
condatainer run --afterok "$JOB" --local src/run_dep.sh # should work (just testing the HTCondor like job ID)
JOB=12345.node12.school.edu
condatainer run --afterok "$JOB" --local src/run_dep.sh # should work (just testing the PBS like job ID)

rm src/logs/long_job_output.txt

# run array test
JOB=$(condatainer run --array src/array.txt --array-limit 2 src/run_array.sh --testing this)
condatainer run --afterok "$JOB" --array src/array.txt src/run_array.sh another test


