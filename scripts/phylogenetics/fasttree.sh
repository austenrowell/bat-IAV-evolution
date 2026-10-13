#!/bin/bash
#SBATCH -J fasttree
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-7
#SBATCH -p nocona
#SBATCH -o logs/fasttree.%A_%a.out
#SBATCH -e logs/fasttree.%A_%a.err
#SBATCH --export=ALL
#SBATCH --mail-user=austen.rowell@ttu.edu
#SBATCH --mail-type=ALL

# Load FastTree module or activate environment
#module load fasttree  # Or use conda activate if needed

# Define the segment folders
segments=(segment_1 segment_2 segment_3 segment_4 segment_5 segment_6 segment_7 segment_8)

# Select the segment based on the SLURM array index
segment=${segments[$SLURM_ARRAY_TASK_ID]}

# Define input and output files
input_file="${segment}/${segment}_trim.fasta"
output_file="${segment}/${segment}_tree.nwk"

# Check if aligned input file exists
if [[ -f "$input_file" ]]; then
    echo "[$segment] Building tree with FastTree..."
    FastTreeMP -nt -gtr -fastest "$input_file" > "$output_file"
    echo "[$segment] Tree written to $output_file"
else
    echo "[$segment] Aligned file $input_file not found. Skipping."
fi

