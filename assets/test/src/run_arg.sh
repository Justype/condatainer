#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=2GB
#SBATCH --output=src/logs/run_arg_%j.out

echo run_arg.sh: "$@"
echo number of args: $#
