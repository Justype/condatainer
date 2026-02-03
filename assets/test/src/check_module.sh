#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB

module purge
module load samtools/1.22.1 bcftools/1.22 seqtk/1.5

bcftools --version
