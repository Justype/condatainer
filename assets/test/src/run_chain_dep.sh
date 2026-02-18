#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --output=log/run_chain_dep.out
#DEP: salmon/1.10.2
#DEP: grcm39/salmon/1.10.2/gencodeM6

echo IN_CONDATAINER "$IN_CONDATAINER"
echo SALMON_INDEX_DIR $SALMON_INDEX_DIR
ls -l $SALMON_INDEX_DIR
salmon --version
