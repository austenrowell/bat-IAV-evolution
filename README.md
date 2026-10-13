# Bat Influenza A Virus Analyses

Scripts for phylogenetic and sequence analyses of bat-associated
influenza A viruses.

## Overview

This repository contains workflows for:

- Selection analyses using PAML and HyPhy.
- Codon-usage analyses using coRdon.
- Preparing analysis inputs and summarizing results.

## Repository contents

- `scripts/` — Input preparation, analysis, and plotting scripts.
- `data/` — Input data or instructions for obtaining the data.
- `results/` — Summary tables and selected analysis outputs.

## Requirements

The analyses use Python, R, PAML, HyPhy, and coRdon.
Individual scripts specify additional package requirements.

Shell submission scripts are configured for a SLURM cluster.
Update working directories, environment names, and cluster settings
before running them.

## Input data

Analyses require coding-sequence FASTA files, corresponding
phylogenetic trees, and sequence metadata.

The metadata file `bat_sequences.csv` identifies bat-associated
sequences through the `Accession` column. Codon-usage scripts can
also use a `Host` column for group comparisons.

Input paths and expected file structures are documented in the
individual workflows.

## Running the analyses

1. Install the required software and packages.
2. Obtain the input sequences, trees, and metadata.
3. Update the input and output paths.
4. Run the relevant input-preparation script.
5. Run the analysis locally or submit the cluster jobs.
6. Run the summary and plotting scripts.

Reduced-dataset and large-dataset analyses use separate workflows.
See the README included with each workflow for details.

## Outputs

Outputs include selection-test summaries, codon-usage metrics,
statistical comparisons, and figures. Output filenames and locations
are specified in the scripts.

## Reproducibility

Record software versions, analysis settings, and the input dataset
used for each run. Retain accession lists and metadata needed to
identify the analyzed sequences.

## Citation

TBA

## License

See `LICENSE` for the terms of reuse.

## Contact

For questions or problems, open a GitHub issue.
