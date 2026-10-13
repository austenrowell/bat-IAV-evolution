#!/bin/bash
#SBATCH --job-name=iqtree_subsets
#SBATCH --output=logs/iqtree_%A_%a.out
#SBATCH --error=logs/iqtree_%A_%a.err
#SBATCH --array=1-8
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --export=ALL
#SBATCH --mail-user=austen.rowell@ttu.edu
#SBATCH --mail-type=ALL

SEG_NUM=$SLURM_ARRAY_TASK_ID
SEG="segment_${SEG_NUM}"

ALIGN="${SEG}/${SEG}_subset_trim.fasta"
PREFIX="${SEG}/${SEG}"

echo "[$SEG] Running IQ-TREE on $ALIGN"

iqtree \
  -s "$ALIGN" \
  -pre "$PREFIX" \
  -bb 10000 \
  -bnni \
  -m GTR+F+I \
  -nt AUTO -redo

echo "[$SEG] Done."

