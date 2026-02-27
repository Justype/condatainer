#!/bin/bash
#SBATCH --job-name=run_array
#SBATCH --output=src/logs/run_array.out
##SBATCH --error=src/logs/run_array.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=100M
#SBATCH --time=00:10:00

# Run the array job
echo n args $#
echo args $@
for arg in "$@"; do
    sleep 10
    echo "Processing argument: $arg"
done
