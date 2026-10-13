#!/bin/bash
#SBATCH -J align
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --array=0-7
#SBATCH -p nocona
#SBATCH -o logs/align.%A_%a.out
#SBATCH -e logs/align.%A_%a.err
#SBATCH --export=ALL
#SBATCH --mail-user=austen.rowell@ttu.edu
#SBATCH --mail-type=ALL

#module load mafft  # or activate conda env if needed

segments=(segment_1 segment_2 segment_3 segment_4 segment_5 segment_6 segment_7 segment_8)

segment=${segments[$SLURM_ARRAY_TASK_ID]}
input_file="${segment}/${segment}_subset_threshold98.fasta"
output_file="${segment}/${segment}_subset_align.fasta"

if [[ -f "$input_file" ]]; then
    echo "[$segment] Aligning $input_file..."
    mafft --thread 2 --auto "$input_file" > "$output_file"
    echo "[$segment] Done."
else
    echo "[$segment] Input file $input_file not found. Skipping."
fi
