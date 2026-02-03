#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#DEP:samtools/1.22.1
#DEP:bcftools/1.22
#DEP:seqtk/1.5

samtools --version
