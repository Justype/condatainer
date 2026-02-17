#!/bin/bash
#BSUB -W 2:00
#BSUB -n 1
#BSUB -M 1G
#BSUB -R "span[hosts=1]"
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
