#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=100
#SBATCH --output=src/logs/run_long.out

echo "This is a long-running job. Sleeping for 100 seconds..."
echo "Job started at: $(date)" > src/logs/long_job_output.txt
sleep 100
echo "Job finished at: $(date)" >> src/logs/long_job_output.txt
echo "Job completed."
