#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB

ml purge
ml samtools/1.22.1 bcftools/1.22
ml seqtk/1.5

seqtk
