#!/bin/bash
#SBATCH --job-name=beast75_50M
#SBATCH --output=logs/beast75_50M_%A_%a.out
#SBATCH --error=logs/beast75_50M_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL

set -euo pipefail
IFS=$'\n\t'

cd /lustre/work/aurowell/influenza/cds/beast_ctmc_75_subset

BEAST_EXE="/lustre/work/aurowell/influenza/cds/beast_ctmc_75_subset/BEASTv10.5.0/bin/beast"

XMLS=(
  "xml_50M/segment_3_PA_75subset_lsd2.timetree_50M.xml"
  "xml_50M/segment_7_M2_75subset_lsd2.timetree_50M.xml"
)

XML="${XMLS[$SLURM_ARRAY_TASK_ID]}"
XML_BASE="$(basename "$XML" .xml)"

# Remove the _50M suffix to recover the original tree basename.
TREE_BASE="${XML_BASE%_50M}"

RUN_DIR="runs_50M/${XML_BASE}"
mkdir -p "$RUN_DIR"

cp "$XML" "$RUN_DIR/$(basename "$XML")"
cp "input_trees/${TREE_BASE}.nex" "$RUN_DIR/${TREE_BASE}.nex"

cd "$RUN_DIR"

echo "Started: $(date)"
echo "Node: $(hostname)"
echo "XML: $(basename "$XML")"
echo "Tree: ${TREE_BASE}.nex"
echo "Run directory: $RUN_DIR"

export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

"$BEAST_EXE" \
  -overwrite \
  -beagle_CPU \
  -beagle_double \
  -beagle_threading_off \
  -threads "$SLURM_CPUS_PER_TASK" \
  "$(basename "$XML")"

echo "Finished: $(date)"
