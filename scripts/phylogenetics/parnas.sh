#!/bin/bash
#SBATCH --job-name=parnas_thresh98
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -p nocona
#SBATCH --output=logs/parnas_thresh98_%A_%a.out
#SBATCH --error=logs/parnas_thresh98_%A_%a.err
#SBATCH --array=1-8
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G


# Load or activate Parnas environment if needed:
# module load parnas
# source activate parnas_env

# Determine segment
SEG_NUM=$SLURM_ARRAY_TASK_ID
SEG="segment_${SEG_NUM}"

# Input files
TREE="${SEG}/${SEG}_tree.nwk"
ALIGN="${SEG}/${SEG}_align.fasta"

# Output file (list of representatives)
OUT="${SEG}/${SEG}_threshold98_reps.txt"

echo "[$SEG] Running PARNAS with threshold=98"

# Parnas invocation
parnas \
  -t "$TREE" \
  --nt "$ALIGN" \
  --threshold 98 \
  --cover > "$OUT"

echo "[$SEG] Representatives written to $OUT"
