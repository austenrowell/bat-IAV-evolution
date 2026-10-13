Bat IAV codon-usage analysis


Reduced 75-taxon datasets (June 2, 2026):
  run_cordon_codon_usage_validation_with_metric_tests_FIXED.R
  Input root: realign_inputs_combined_75
  Output: codon_bias_results_75_cordon

Large dataset (June 3, 2026):
  01_prepare_large_hyphy_codon_inputs.py
  04_run_large_cordon_codon_usage.R
  05_submit_large_cordon_codon_usage.sh
  Input root for R analysis: large_consistent_inputs
  Output: codon_bias_results_large_cordon

Both R analyses use bat_sequences.csv and coRdon-derived codon counts, ENC,
RSCU, GC3 and total GC, with PERMANOVA, dispersion and metric tests.
Review input paths and the installed R environment before execution.
The optional comparison with earlier Python results uses existing output files.
