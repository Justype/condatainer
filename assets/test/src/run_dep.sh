#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=100
#SBATCH --output=src/logs/run_dep.out

# This script will run after the long job (src/run_long.sh) finishes, testing the --afterok dependency feature.

if [ -s src/logs/long_job_output.txt ]; then
    echo "Dependency test passed: Long job output is present."
    echo "Contents of long job output:"
    cat src/logs/long_job_output.txt
else
    echo "Dependency test failed: Long job output is missing or empty."
fi
