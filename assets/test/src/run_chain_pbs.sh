#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=1gb
#DEP: salmon/1.10.2
#DEP: grcm39/salmon/1.10.2/gencodeM6

echo IN_CONDATAINER "$IN_CONDATAINER"
echo SALMON_INDEX_DIR $SALMON_INDEX_DIR
ls -l $SALMON_INDEX_DIR
echo NNODES $NNODES
echo NTASKS $NTASKS
echo NCPUS $NCPUS
echo MEM $MEM
echo MEM MB $MEM_MB
echo MEM GB $MEM_GB
salmon --version
