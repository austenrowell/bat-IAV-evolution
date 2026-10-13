#!/bin/bash
#SBATCH --job-name=iqtree_time_scaled
#SBATCH --output=logs/iqtree_time_scaled_%A_%a.out
#SBATCH --error=logs/iqtree_time_scaled_%A_%a.err
#SBATCH --array=0-10
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL

set -euo pipefail
IFS=$'\n\t'

# list of segment dirs
SEGMENTS=(
  segment_1
  segment_2
  segment_3_PA
  segment_4
  segment_5
  segment_6
  segment_7_M1
  segment_7_M2
  segment_8_NEP
  segment_8_NS1
  segment_8_NS2
)

SEG="${SEGMENTS[$SLURM_ARRAY_TASK_ID]}"

ALIGN="${SEG}/${SEG}_clean.fasta"
DATES="${SEG}/${SEG}_dates_subset.txt"

# OPTIONAL: set OUTGROUP (comma-separated taxa exactly matching alignment names)
# OUTGROUP=""                # e.g. OUTGROUP="bat_taxon,pangolin_taxon"
OUTGROUP=""

# number of CI replicates (set lower for test runs)
CI_REPS=1000

echo "[${SEG}] START `date`"
echo "[${SEG}] ALIGN=${ALIGN}"
echo "[${SEG}] DATES=${DATES}"

# sanity checks
if [ ! -f "$ALIGN" ]; then
  echo "[${SEG}] ERROR: alignment not found: $ALIGN"
  exit 1
fi
if [ ! -f "$DATES" ]; then
  echo "[${SEG}] ERROR: dates file not found: $DATES"
  exit 1
fi

# build command: infer ML tree (ModelFinder + ML) and then run LSD2 dating in one step
# This follows the documented workflow: iqtree -s ALN --date DATE_FILE
CMD=( iqtree -s "$ALIGN" \
         --date "$DATES" \
         --date-ci "$CI_REPS" \
         -m GTR+F+I \
         -alrt 1000 \
         -bb 10000 \
         -bnni \
         -nt AUTO \
         -redo )

# add outgroup option if provided
if [ -n "$OUTGROUP" ]; then
  CMD+=( -o "$OUTGROUP" )
fi

echo "[${SEG}] Running IQ-TREE (inference + LSD2 dating)"
echo "[${SEG}] Command: ${CMD[*]}"

"${CMD[@]}" || { echo "[${SEG}] ERROR: iqtree failed"; exit 1; }

# outputs:
# - ${ALIGN}.timetree.nwk (Newick time tree)
# - ${ALIGN}.timetree.nex (Nexus time tree)
# - ${ALIGN}.timetree.lsd (LSD2 report)
echo "[${SEG}] DONE `date`"