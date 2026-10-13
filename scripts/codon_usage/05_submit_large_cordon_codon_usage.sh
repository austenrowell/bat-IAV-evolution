#!/bin/bash
#SBATCH --job-name=large_cordon
#SBATCH --output=logs/large_cordon_%j.out
#SBATCH --error=logs/large_cordon_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL
#SBATCH --time=12:00:00

set -euo pipefail

cd /lustre/work/aurowell/influenza/cds

mkdir -p logs

# Rscript must be available in the job environment.
# Uncomment/adapt these lines on HPCC if R is not already available:
# module purge
# module load ansys/v241
# module load gcc/14.2.0
# module load R/4.4.1

Rscript 04_run_large_cordon_codon_usage.R
