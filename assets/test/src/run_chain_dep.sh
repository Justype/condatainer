#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --output=src/logs/run_chain_dep.out
#DEP: salmon/1.10.2
#DEP: grcm39/salmon/1.10.2/gencodeM6

echo IN_CONDATAINER "$IN_CONDATAINER"
echo SALMON_INDEX_DIR $SALMON_INDEX_DIR
ls -l $SALMON_INDEX_DIR
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
salmon --version
