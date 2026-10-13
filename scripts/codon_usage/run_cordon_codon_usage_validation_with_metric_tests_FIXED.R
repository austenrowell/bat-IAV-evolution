#!/usr/bin/env Rscript

############################################################
# coRdon-based codon-usage validation + statistical analysis
#
# Purpose:
#   Uses Bioconductor coRdon to read FASTA files, build codonTable
#   objects, extract codon counts, and calculate coRdon ENC values.
#   RSCU, GC3, and GC_total are calculated from the coRdon-derived
#   codon counts.
#
#   Then performs:
#     1. adonis2 PERMANOVA on scaled RSCU profiles by gene
#     2. betadisper tests for homogeneity of RSCU dispersion
#     3. Wilcoxon/Kruskal-Wallis tests for GC3, GC_total, ENC_cordon,
#        and valid_codons
#     4. Per-codon RSCU Wilcoxon/Kruskal-Wallis tests
#     5. BH correction for multiple testing
#
# Input:
#   realign_inputs_combined_75/*/segment_*/segment_*.retained_original_ungapped.fasta
#   bat_sequences.csv
#
# Output:
#   codon_bias_results_75_cordon/
############################################################

options(stringsAsFactors = FALSE)

############################################################
# 0. User settings
############################################################

fasta_root <- "realign_inputs_combined_75"
fasta_pattern <- "*/segment_*/segment_*.retained_original_ungapped.fasta"

metadata_file <- "bat_sequences.csv"
bat_csv <- "bat_sequences.csv"

accession_col <- "Accession"
host_col <- "Host"

outdir <- "codon_bias_results_75_cordon"

min_codons <- 30
permutations <- 999
set.seed(1)

############################################################
# 1. Install/load packages
############################################################

cran_packages <- c("tidyverse", "readr", "stringr", "vegan")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_packages <- c("Biostrings", "coRdon")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(vegan)
  library(Biostrings)
  library(coRdon)
})

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

############################################################
# 2. Genetic-code definitions for RSCU/GC calculations
############################################################

codon_to_aa <- c(
  TTT="F", TTC="F", TTA="L", TTG="L",
  TCT="S", TCC="S", TCA="S", TCG="S",
  TAT="Y", TAC="Y", TAA="*", TAG="*",
  TGT="C", TGC="C", TGA="*", TGG="W",
  CTT="L", CTC="L", CTA="L", CTG="L",
  CCT="P", CCC="P", CCA="P", CCG="P",
  CAT="H", CAC="H", CAA="Q", CAG="Q",
  CGT="R", CGC="R", CGA="R", CGG="R",
  ATT="I", ATC="I", ATA="I", ATG="M",
  ACT="T", ACC="T", ACA="T", ACG="T",
  AAT="N", AAC="N", AAA="K", AAG="K",
  AGT="S", AGC="S", AGA="R", AGG="R",
  GTT="V", GTC="V", GTA="V", GTG="V",
  GCT="A", GCC="A", GCA="A", GCG="A",
  GAT="D", GAC="D", GAA="E", GAG="E",
  GGT="G", GGC="G", GGA="G", GGG="G"
)

sense_codons <- sort(names(codon_to_aa)[codon_to_aa != "*"])
stop_codons <- sort(names(codon_to_aa)[codon_to_aa == "*"])

aa_to_codons <- split(names(codon_to_aa)[codon_to_aa != "*"], codon_to_aa[codon_to_aa != "*"])
aa_to_codons <- lapply(aa_to_codons, sort)

segment_to_gene <- c(
  segment_1 = "PB2",
  segment_2 = "PB1",
  segment_3_PA = "PA",
  segment_4 = "HA",
  segment_5 = "NP",
  segment_6 = "NA",
  segment_7_M1 = "M1",
  segment_7_M2 = "M2",
  segment_8_NS1 = "NS1",
  segment_8_NS2 = "NS2/NEP",
  segment_8_NEP = "NS2/NEP"
)

############################################################
# 3. Helper functions
############################################################

normalize_accession <- function(x) {
  if (is.na(x) || is.null(x)) return("")
  s <- str_trim(as.character(x))
  if (s == "") return("")
  s <- str_split(s, "\\s+", simplify = TRUE)[1]
  if (str_detect(s, "\\|")) {
    parts <- unlist(str_split(s, "\\|"))
    parts <- parts[parts != ""]
    acc_like <- parts[str_detect(parts, "^[A-Z]{1,3}_?\\d+\\.\\d+$|^[A-Z]{1,3}\\d+\\.\\d+$")]
    if (length(acc_like) > 0) {
      s <- acc_like[1]
    } else if (length(parts) > 0) {
      s <- parts[1]
    }
  }
  s
}

accession_variants <- function(x) {
  s <- normalize_accession(x)
  vars <- unique(c(
    s,
    str_replace_all(s, "\\.", "_"),
    str_replace_all(s, "_", "."),
    str_replace(s, "[._]\\d+$", "")
  ))
  vars[!is.na(vars) & vars != ""]
}

infer_segment <- function(path, header = "") {
  text <- paste(path, header)
  ordered_segments <- names(segment_to_gene)[order(nchar(names(segment_to_gene)), decreasing = TRUE)]
  for (seg in ordered_segments) {
    if (str_detect(text, fixed(seg))) return(seg)
  }
  m <- str_extract(text, "segment_\\d+(?:_[A-Za-z0-9]+)?")
  if (!is.na(m)) return(m)
  basename(dirname(path))
}

load_table <- function(path) {
  if (is.null(path) || is.na(path) || !file.exists(path)) return(tibble())
  if (str_detect(tolower(path), "\\.(tsv|txt)$")) {
    read_tsv(path, show_col_types = FALSE, na = c("", "NA", "NaN"))
  } else {
    read_csv(path, show_col_types = FALSE, na = c("", "NA", "NaN"))
  }
}

make_metadata_lookup <- function(metadata, accession_col) {
  lookup <- new.env(parent = emptyenv())
  if (nrow(metadata) == 0 || !(accession_col %in% colnames(metadata))) return(lookup)
  for (i in seq_len(nrow(metadata))) {
    acc <- metadata[[accession_col]][i]
    for (v in accession_variants(acc)) {
      assign(v, as.list(metadata[i, , drop = FALSE]), envir = lookup)
    }
  }
  lookup
}

lookup_metadata_row <- function(variants, lookup) {
  for (v in variants) {
    if (exists(v, envir = lookup, inherits = FALSE)) return(get(v, envir = lookup, inherits = FALSE))
  }
  NULL
}

safe_host_group <- function(md_row, host_col, is_bat) {
  if (!is.null(md_row) && !is.null(host_col) && host_col %in% names(md_row)) {
    raw <- as.character(md_row[[host_col]])
    if (!is.na(raw) && str_trim(raw) != "") {
      if (str_detect(tolower(raw), "bat")) return("Bat")
      return(raw)
    }
  }
  ifelse(is_bat, "Bat", "Non_bat")
}

calc_rscu_from_counts <- function(counts_named) {
  out <- setNames(rep(NA_real_, length(sense_codons)), sense_codons)
  for (aa in names(aa_to_codons)) {
    codons <- aa_to_codons[[aa]]
    total <- sum(counts_named[codons], na.rm = TRUE)
    n <- length(codons)
    if (total > 0) {
      expected <- total / n
      out[codons] <- counts_named[codons] / expected
    }
  }
  out
}

calc_gc_metrics_from_counts <- function(counts_named) {
  total_codons <- sum(counts_named[sense_codons], na.rm = TRUE)
  if (total_codons == 0) return(c(GC_total = NA_real_, GC3 = NA_real_))
  gc_total_count <- 0
  gc3_count <- 0
  for (codon in sense_codons) {
    n <- counts_named[[codon]]
    if (is.na(n) || n == 0) next
    nts <- unlist(strsplit(codon, ""))
    gc_total_count <- gc_total_count + n * sum(nts %in% c("G", "C"))
    gc3_count <- gc3_count + n * as.integer(nts[3] %in% c("G", "C"))
  }
  c(GC_total = gc_total_count / (3 * total_codons), GC3 = gc3_count / total_codons)
}

safe_enc <- function(ct) {
  out <- tryCatch(
    as.numeric(coRdon::ENC(ct, filtering = "soft")),
    error = function(e1) {
      tryCatch(
        as.numeric(coRdon::ENC(ct)),
        error = function(e2) rep(NA_real_, length(ct))
      )
    }
  )
  out
}

prepare_rscu_matrix <- function(df, codon_cols) {
  mat <- df %>% select(all_of(codon_cols)) %>% mutate(across(everything(), as.numeric)) %>% as.data.frame()
  for (col in colnames(mat)) {
    col_mean <- mean(mat[[col]], na.rm = TRUE)
    if (is.nan(col_mean) || !is.finite(col_mean)) col_mean <- 0
    mat[[col]][is.na(mat[[col]])] <- col_mean
  }
  keep_cols <- sapply(mat, function(x) sd(x, na.rm = TRUE) > 0)
  mat <- mat[, keep_cols, drop = FALSE]
  if (ncol(mat) < 2) return(NULL)
  scale(mat)
}

run_adonis2_for_gene <- function(gene_name, df, codon_cols, permutations) {
  sub <- df %>% filter(gene == gene_name) %>% filter(!is.na(host_group)) %>% mutate(host_group = factor(host_group))
  n_sequences <- nrow(sub)
  n_groups <- n_distinct(sub$host_group)
  groups <- paste(sort(unique(as.character(sub$host_group))), collapse = ";")
  if (n_sequences < 4 || n_groups < 2) {
    return(tibble(gene = gene_name, n_sequences = n_sequences, n_groups = n_groups, groups = groups,
                  method = "adonis2", distance = "euclidean_on_scaled_RSCU", permutations = permutations,
                  df_model = NA_real_, df_residual = NA_real_, sumsq_model = NA_real_, sumsq_residual = NA_real_,
                  r2 = NA_real_, pseudo_F = NA_real_, p_value = NA_real_, status = "not_run_insufficient_groups_or_sequences"))
  }
  mat <- prepare_rscu_matrix(sub, codon_cols)
  if (is.null(mat)) {
    return(tibble(gene = gene_name, n_sequences = n_sequences, n_groups = n_groups, groups = groups,
                  method = "adonis2", distance = "euclidean_on_scaled_RSCU", permutations = permutations,
                  df_model = NA_real_, df_residual = NA_real_, sumsq_model = NA_real_, sumsq_residual = NA_real_,
                  r2 = NA_real_, pseudo_F = NA_real_, p_value = NA_real_, status = "not_run_too_few_variable_codons"))
  }
  d <- dist(mat, method = "euclidean")
  ad <- vegan::adonis2(d ~ host_group, data = sub, permutations = permutations, by = "terms")
  ad_tbl <- as.data.frame(ad)
  tibble(gene = gene_name, n_sequences = n_sequences, n_groups = n_groups, groups = groups,
         method = "adonis2", distance = "euclidean_on_scaled_RSCU", permutations = permutations,
         df_model = ad_tbl$Df[1], df_residual = ad_tbl$Df[2], sumsq_model = ad_tbl$SumOfSqs[1],
         sumsq_residual = ad_tbl$SumOfSqs[2], r2 = ad_tbl$R2[1], pseudo_F = ad_tbl$F[1],
         p_value = ad_tbl$`Pr(>F)`[1], status = "OK")
}

run_betadisper_for_gene <- function(gene_name, df, codon_cols, permutations) {
  sub <- df %>% filter(gene == gene_name) %>% filter(!is.na(host_group)) %>% mutate(host_group = factor(host_group))
  if (nrow(sub) < 4 || n_distinct(sub$host_group) < 2) {
    return(tibble(gene = gene_name, dispersion_F = NA_real_, dispersion_p_value = NA_real_,
                  dispersion_status = "not_run_insufficient_groups_or_sequences"))
  }
  mat <- prepare_rscu_matrix(sub, codon_cols)
  if (is.null(mat)) {
    return(tibble(gene = gene_name, dispersion_F = NA_real_, dispersion_p_value = NA_real_,
                  dispersion_status = "not_run_too_few_variable_codons"))
  }
  d <- dist(mat, method = "euclidean")
  bd <- vegan::betadisper(d, sub$host_group)
  bd_perm <- vegan::permutest(bd, permutations = permutations)
  bd_tab <- as.data.frame(bd_perm$tab)
  tibble(gene = gene_name, dispersion_F = bd_tab$F[1], dispersion_p_value = bd_tab$`Pr(>F)`[1], dispersion_status = "OK")
}

run_group_metric_tests <- function(seq_df) {
  metrics_to_test <- c("GC3", "GC_total", "ENC_cordon", "valid_codons")
  missing_metrics <- setdiff(metrics_to_test, colnames(seq_df))
  if (length(missing_metrics) > 0) warning("Missing metrics skipped: ", paste(missing_metrics, collapse = ", "))
  metrics_to_test <- intersect(metrics_to_test, colnames(seq_df))
  results <- list()
  for (gene_name in sort(unique(seq_df$gene))) {
    sub_gene <- seq_df %>% filter(gene == gene_name) %>% filter(!is.na(host_group))
    for (metric in metrics_to_test) {
      sub <- sub_gene %>%
        transmute(gene = gene, host_group = host_group, value = .data[[metric]]) %>%
        filter(!is.na(value)) %>%
        mutate(host_group = factor(host_group))
      n_groups <- n_distinct(sub$host_group)
      groups <- paste(sort(unique(as.character(sub$host_group))), collapse = ";")
      if (nrow(sub) < 3 || n_groups < 2) {
        results[[length(results) + 1]] <- tibble(gene = gene_name, metric = metric, n_groups = n_groups,
                                                 groups = groups, test = "not_run_insufficient_groups_or_sequences",
                                                 statistic = NA_real_, p_value = NA_real_)
        next
      }
      row <- tibble(gene = gene_name, metric = metric, n_groups = n_groups, groups = groups)
      group_summary <- sub %>% group_by(host_group) %>%
        summarize(n = n(), mean = mean(value, na.rm = TRUE), median = median(value, na.rm = TRUE), .groups = "drop")
      for (i in seq_len(nrow(group_summary))) {
        grp <- as.character(group_summary$host_group[i])
        row[[paste0("n_", grp)]] <- group_summary$n[i]
        row[[paste0("mean_", grp)]] <- group_summary$mean[i]
        row[[paste0("median_", grp)]] <- group_summary$median[i]
      }
      if (n_groups == 2) {
        wt <- wilcox.test(value ~ host_group, data = sub, exact = FALSE)
        row$test <- "Wilcoxon rank-sum"
        row$statistic <- unname(wt$statistic)
        row$p_value <- wt$p.value
      } else {
        kt <- kruskal.test(value ~ host_group, data = sub)
        row$test <- "Kruskal-Wallis"
        row$statistic <- unname(kt$statistic)
        row$p_value <- kt$p.value
      }
      results[[length(results) + 1]] <- row
    }
  }
  out <- bind_rows(results)
  if (nrow(out) > 0 && "p_value" %in% colnames(out)) {
    out <- out %>% group_by(metric) %>% mutate(q_value_BH_within_metric = p.adjust(p_value, method = "BH")) %>%
      ungroup() %>% mutate(q_value_BH_all_metric_tests = p.adjust(p_value, method = "BH"),
                           significant_q_0.05_all_metric_tests = case_when(
                             is.na(q_value_BH_all_metric_tests) ~ NA,
                             q_value_BH_all_metric_tests < 0.05 ~ TRUE,
                             TRUE ~ FALSE))
  }
  out
}

run_per_codon_tests <- function(rscu_df, codon_cols) {
  results <- list()
  for (gene_name in sort(unique(rscu_df$gene))) {
    sub_gene <- rscu_df %>% filter(gene == gene_name) %>% filter(!is.na(host_group)) %>% mutate(host_group = factor(host_group))
    if (nrow(sub_gene) < 3 || n_distinct(sub_gene$host_group) < 2) next
    for (codon in codon_cols) {
      sub <- sub_gene %>%
        transmute(gene = gene, host_group = host_group, value = .data[[codon]]) %>%
        filter(!is.na(value))
      n_groups <- n_distinct(sub$host_group)
      if (nrow(sub) < 3 || n_groups < 2) next
      row <- tibble(gene = gene_name, codon = codon, amino_acid = codon_to_aa[[codon]], n_groups = n_groups,
                    groups = paste(sort(unique(as.character(sub$host_group))), collapse = ";"))
      group_summary <- sub %>% group_by(host_group) %>%
        summarize(n = n(), mean_RSCU = mean(value, na.rm = TRUE), median_RSCU = median(value, na.rm = TRUE), .groups = "drop")
      for (i in seq_len(nrow(group_summary))) {
        grp <- as.character(group_summary$host_group[i])
        row[[paste0("n_", grp)]] <- group_summary$n[i]
        row[[paste0("mean_RSCU_", grp)]] <- group_summary$mean_RSCU[i]
        row[[paste0("median_RSCU_", grp)]] <- group_summary$median_RSCU[i]
      }
      if (n_groups == 2) {
        wt <- wilcox.test(value ~ host_group, data = sub, exact = FALSE)
        row$test <- "Wilcoxon rank-sum"
        row$statistic <- unname(wt$statistic)
        row$p_value <- wt$p.value
      } else {
        kt <- kruskal.test(value ~ host_group, data = sub)
        row$test <- "Kruskal-Wallis"
        row$statistic <- unname(kt$statistic)
        row$p_value <- kt$p.value
      }
      results[[length(results) + 1]] <- row
    }
  }
  out <- bind_rows(results)
  if (nrow(out) > 0 && "p_value" %in% colnames(out)) {
    out <- out %>% group_by(gene) %>% mutate(q_value_BH_within_gene = p.adjust(p_value, method = "BH")) %>%
      ungroup() %>% mutate(q_value_BH_all_codon_tests = p.adjust(p_value, method = "BH"),
                           significant_q_0.05_all_codon_tests = case_when(
                             is.na(q_value_BH_all_codon_tests) ~ NA,
                             q_value_BH_all_codon_tests < 0.05 ~ TRUE,
                             TRUE ~ FALSE))
  }
  out
}

############################################################
# 4. Load metadata and collect FASTA files
############################################################

metadata <- load_table(metadata_file)
bat_df <- load_table(bat_csv)
metadata_lookup <- make_metadata_lookup(metadata, accession_col)

bat_accessions <- character()
if (nrow(bat_df) > 0 && accession_col %in% colnames(bat_df)) {
  bat_accessions <- unique(unlist(lapply(bat_df[[accession_col]], accession_variants)))
}

fasta_files <- Sys.glob(file.path(fasta_root, fasta_pattern))
fasta_files <- fasta_files[file.exists(fasta_files)]
if (length(fasta_files) == 0) stop("No FASTA files found with: ", file.path(fasta_root, fasta_pattern))

############################################################
# 5. Process FASTA files with coRdon
############################################################

seq_rows <- list()
codon_rows <- list()
rscu_rows <- list()

for (fasta in fasta_files) {
  message("Processing: ", fasta)
  dna <- coRdon::readSet(file = fasta)
  ct <- coRdon::codonTable(dna)
  counts <- as.data.frame(coRdon::codonCounts(ct))
  for (codon in sense_codons) if (!(codon %in% colnames(counts))) counts[[codon]] <- 0
  counts <- counts[, unique(c(sense_codons, intersect(colnames(counts), stop_codons))), drop = FALSE]

  seq_ids <- names(dna)
  if (is.null(seq_ids) || length(seq_ids) != nrow(counts)) seq_ids <- rownames(counts)
  if (is.null(seq_ids) || length(seq_ids) != nrow(counts)) seq_ids <- paste0("seq_", seq_len(nrow(counts)))

  enc_vals <- safe_enc(ct)
  if (length(enc_vals) != nrow(counts)) enc_vals <- rep(NA_real_, nrow(counts))

  for (i in seq_len(nrow(counts))) {
    header <- seq_ids[i]
    accession <- normalize_accession(header)
    variants <- accession_variants(accession)
    md_row <- lookup_metadata_row(variants, metadata_lookup)
    is_bat <- any(variants %in% bat_accessions)
    host_group <- safe_host_group(md_row, host_col, is_bat)
    segment <- infer_segment(fasta, header)
    gene <- ifelse(segment %in% names(segment_to_gene), segment_to_gene[[segment]], segment)

    counts_named <- as.numeric(counts[i, sense_codons, drop = TRUE])
    names(counts_named) <- sense_codons
    valid_codons <- sum(counts_named, na.rm = TRUE)
    if (valid_codons < min_codons) next

    rscu <- calc_rscu_from_counts(counts_named)
    gc_metrics <- calc_gc_metrics_from_counts(counts_named)

    base <- tibble(sequence_id = header, accession = accession, segment = segment, gene = gene,
                   host_group = host_group, is_bat_accession = is_bat, fasta_file = fasta)

    seq_rows[[length(seq_rows) + 1]] <- bind_cols(base, tibble(valid_codons = valid_codons,
                                                               GC_total = unname(gc_metrics["GC_total"]),
                                                               GC3 = unname(gc_metrics["GC3"]),
                                                               ENC_cordon = enc_vals[i]))
    codon_rows[[length(codon_rows) + 1]] <- bind_cols(base, as_tibble(as.list(counts_named)))
    rscu_rows[[length(rscu_rows) + 1]] <- bind_cols(base, as_tibble(as.list(rscu)))
  }
}

seq_df <- bind_rows(seq_rows)
codon_df <- bind_rows(codon_rows)
rscu_df <- bind_rows(rscu_rows)
if (nrow(seq_df) == 0) stop("No sequences retained after min_codons filtering.")

############################################################
# 6. Write coRdon-derived feature outputs
############################################################

write_csv(seq_df, file.path(outdir, "cordon_sequence_codon_bias_metrics.csv"))
write_csv(codon_df, file.path(outdir, "cordon_codon_counts.csv"))
write_csv(rscu_df, file.path(outdir, "cordon_rscu_matrix.csv"))

############################################################
# 7. Run statistical analyses
############################################################

codon_cols <- intersect(sense_codons, colnames(rscu_df))
genes <- rscu_df %>% filter(!is.na(gene)) %>% pull(gene) %>% unique() %>% sort()

adonis2_results <- purrr::map_dfr(genes, ~ run_adonis2_for_gene(.x, rscu_df, codon_cols, permutations)) %>%
  mutate(q_value_BH = p.adjust(p_value, method = "BH"),
         significant_q_0.05 = case_when(is.na(q_value_BH) ~ NA, q_value_BH < 0.05 ~ TRUE, TRUE ~ FALSE)) %>%
  arrange(q_value_BH)

dispersion_results <- purrr::map_dfr(genes, ~ run_betadisper_for_gene(.x, rscu_df, codon_cols, permutations)) %>%
  mutate(dispersion_q_value_BH = p.adjust(dispersion_p_value, method = "BH"),
         dispersion_significant_q_0.05 = case_when(is.na(dispersion_q_value_BH) ~ NA,
                                                   dispersion_q_value_BH < 0.05 ~ TRUE, TRUE ~ FALSE))

group_metric_tests <- run_group_metric_tests(seq_df)
per_codon_tests <- run_per_codon_tests(rscu_df, codon_cols)

group_counts <- rscu_df %>% filter(!is.na(gene), !is.na(host_group)) %>%
  count(gene, host_group, name = "n") %>%
  pivot_wider(names_from = host_group, values_from = n, values_fill = 0) %>% arrange(gene)

############################################################
# 8. Build manuscript-style summaries
############################################################

composition_metric_summary <- group_metric_tests %>%
  filter(metric %in% c("GC3", "GC_total", "ENC_cordon")) %>%
  mutate(significant_all_tests = !is.na(q_value_BH_all_metric_tests) & q_value_BH_all_metric_tests < 0.05,
         significant_within_metric = !is.na(q_value_BH_within_metric) & q_value_BH_within_metric < 0.05,
         label_all_tests = if_else(significant_all_tests, metric, NA_character_),
         label_within_metric = if_else(significant_within_metric, metric, NA_character_)) %>%
  group_by(gene) %>%
  summarize(significant_metrics_BH_all_metric_tests = paste(na.omit(unique(label_all_tests)), collapse = "; "),
            significant_metrics_BH_within_metric = paste(na.omit(unique(label_within_metric)), collapse = "; "),
            .groups = "drop") %>%
  mutate(significant_metrics_BH_all_metric_tests = if_else(significant_metrics_BH_all_metric_tests == "", "None", significant_metrics_BH_all_metric_tests),
         significant_metrics_BH_within_metric = if_else(significant_metrics_BH_within_metric == "", "None", significant_metrics_BH_within_metric))

valid_codon_summary <- group_metric_tests %>% filter(metric == "valid_codons") %>%
  transmute(gene, valid_codons_p_value = p_value,
            valid_codons_q_value_BH_all_metric_tests = q_value_BH_all_metric_tests)

per_codon_sig_summary <- per_codon_tests %>% group_by(gene) %>%
  summarize(n_codons_q_lt_0.05_within_gene = sum(q_value_BH_within_gene < 0.05, na.rm = TRUE),
            n_codons_q_lt_0.05_all_codon_tests = sum(q_value_BH_all_codon_tests < 0.05, na.rm = TRUE),
            .groups = "drop")

final_summary <- adonis2_results %>%
  left_join(dispersion_results, by = "gene") %>%
  left_join(group_counts, by = "gene") %>%
  left_join(composition_metric_summary, by = "gene") %>%
  left_join(valid_codon_summary, by = "gene") %>%
  left_join(per_codon_sig_summary, by = "gene") %>%
  mutate(significant_metrics_BH_all_metric_tests = replace_na(significant_metrics_BH_all_metric_tests, "None"),
         significant_metrics_BH_within_metric = replace_na(significant_metrics_BH_within_metric, "None"),
         n_codons_q_lt_0.05_within_gene = replace_na(n_codons_q_lt_0.05_within_gene, 0L),
         n_codons_q_lt_0.05_all_codon_tests = replace_na(n_codons_q_lt_0.05_all_codon_tests, 0L),
         rscu_profile_interpretation = case_when(significant_q_0.05 == TRUE ~ "Significant host-associated RSCU profile difference",
                                                 !is.na(p_value) & p_value < 0.05 & significant_q_0.05 == FALSE ~ "Nominally significant; not significant after BH correction",
                                                 TRUE ~ "Not significant"),
         dispersion_interpretation = case_when(dispersion_significant_q_0.05 == TRUE ~ "Significant dispersion difference; interpret PERMANOVA cautiously",
                                               dispersion_significant_q_0.05 == FALSE ~ "No significant dispersion difference",
                                               TRUE ~ "Dispersion test not available")) %>%
  arrange(q_value_BH)

############################################################
# 9. Optional comparison with original Python-derived outputs
############################################################

python_metrics_file <- "codon_bias_results_75/sequence_codon_bias_metrics.csv"
if (file.exists(python_metrics_file)) {
  py_metrics <- read_csv(python_metrics_file, show_col_types = FALSE, na = c("", "NA", "NaN")) %>%
    mutate(gene = case_when(is.na(gene) & segment == "segment_6" ~ "NA", TRUE ~ as.character(gene)))
  comparison <- seq_df %>%
    select(accession, segment, gene, host_group, valid_codons_cordon = valid_codons,
           GC_total_cordon = GC_total, GC3_cordon = GC3, ENC_cordon) %>%
    left_join(py_metrics %>% select(accession, segment, valid_codons_python = valid_codons,
                                    GC_total_python = GC_total, GC3_python = GC3, ENC_like_python = ENC_like),
              by = c("accession", "segment")) %>%
    mutate(valid_codons_diff = valid_codons_cordon - valid_codons_python,
           GC_total_diff = GC_total_cordon - GC_total_python,
           GC3_diff = GC3_cordon - GC3_python)
  write_csv(comparison, file.path(outdir, "cordon_vs_python_metric_comparison.csv"))
}

############################################################
# 10. Write outputs
############################################################

write_csv(adonis2_results, file.path(outdir, "cordon_adonis2_permanova_rscu_by_gene.csv"))
write_csv(dispersion_results, file.path(outdir, "cordon_betadisper_rscu_by_gene.csv"))
write_csv(group_counts, file.path(outdir, "cordon_host_group_counts_by_gene.csv"))
write_csv(group_metric_tests, file.path(outdir, "cordon_group_metric_tests_by_gene.csv"))
write_csv(per_codon_tests, file.path(outdir, "cordon_per_codon_rscu_tests_by_gene.csv"))
write_csv(final_summary, file.path(outdir, "cordon_final_summary.csv"))

readme_text <- c(
  "coRdon codon-usage validation outputs",
  "=====================================",
  "",
  "This workflow uses Bioconductor coRdon to read FASTA files, build codonTable objects,",
  "extract codon counts, and calculate coRdon ENC values.",
  "",
  "RSCU, GC3, and GC_total are calculated from the coRdon-derived codon counts.",
  "",
  "Statistical tests:",
  "- adonis2 PERMANOVA on scaled RSCU profiles by gene.",
  "- betadisper tests for homogeneity of multivariate RSCU dispersion.",
  "- Wilcoxon rank-sum tests for GC3, GC_total, ENC_cordon, and valid_codons when two host groups are present.",
  "- Kruskal-Wallis tests are used instead if more than two host groups are present.",
  "- Per-codon RSCU group comparisons with BH correction.",
  "",
  "Files:",
  "- cordon_sequence_codon_bias_metrics.csv: per-sequence GC3, GC_total, valid codons, and coRdon ENC.",
  "- cordon_codon_counts.csv: coRdon-derived sense codon counts per sequence.",
  "- cordon_rscu_matrix.csv: RSCU values calculated from coRdon-derived codon counts.",
  "- cordon_adonis2_permanova_rscu_by_gene.csv: adonis2 PERMANOVA on coRdon-derived RSCU profiles.",
  "- cordon_betadisper_rscu_by_gene.csv: dispersion tests for RSCU profiles.",
  "- cordon_group_metric_tests_by_gene.csv: univariate tests for GC3, GC_total, ENC_cordon, and valid_codons.",
  "- cordon_per_codon_rscu_tests_by_gene.csv: per-codon RSCU tests.",
  "- cordon_host_group_counts_by_gene.csv: group counts by gene.",
  "- cordon_final_summary.csv: merged coRdon/adonis2/metric-test summary.",
  "- cordon_vs_python_metric_comparison.csv: optional comparison to original Python outputs if available.",
  "",
  "CAI is not calculated here because no host/reference codon-usage table is supplied."
)
writeLines(readme_text, con = file.path(outdir, "README_cordon_outputs.txt"))

############################################################
# 11. Print summary
############################################################

cat("\nFinished coRdon codon-usage validation and statistical workflow.\n")
cat("Input FASTA files: ", length(fasta_files), "\n", sep = "")
cat("Sequences retained: ", nrow(seq_df), "\n", sep = "")
cat("Output directory: ", outdir, "\n\n", sep = "")

cat("Main outputs:\n")
cat("  cordon_sequence_codon_bias_metrics.csv\n")
cat("  cordon_codon_counts.csv\n")
cat("  cordon_rscu_matrix.csv\n")
cat("  cordon_adonis2_permanova_rscu_by_gene.csv\n")
cat("  cordon_betadisper_rscu_by_gene.csv\n")
cat("  cordon_group_metric_tests_by_gene.csv\n")
cat("  cordon_per_codon_rscu_tests_by_gene.csv\n")
cat("  cordon_final_summary.csv\n\n")

cat("Significant coRdon-derived adonis2 results at BH q < 0.05:\n")
print(
  final_summary %>%
    filter(significant_q_0.05 == TRUE) %>%
    select(gene, n_sequences, any_of(c("Bat", "Non_bat")), r2, pseudo_F, p_value, q_value_BH,
           dispersion_q_value_BH, significant_metrics_BH_all_metric_tests,
           n_codons_q_lt_0.05_all_codon_tests)
)

cat("\nSignificant univariate metric tests at BH q < 0.05 across all metric tests:\n")
print(
  group_metric_tests %>%
    filter(significant_q_0.05_all_metric_tests == TRUE) %>%
    select(gene, metric, any_of(c("n_Bat", "mean_Bat", "median_Bat", "n_Non_bat", "mean_Non_bat", "median_Non_bat")),
           test, statistic, p_value, q_value_BH_all_metric_tests)
)

cat("\nDone.\n")
