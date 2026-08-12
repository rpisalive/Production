# Paired-reference diagnostics and phosphoproteomic validation NanoSPLITS C10 RO-3306 experiment

# Purpose
# -------
# This script characterises the native paired-58 reference before RNA and
# protein correspondence is hidden for SCOT+ alignment. It does not retrain
# MOFA+ and it does not alter the preprocessing matrices.
#
# The script performs four analyses:
#   1. tests native MOFA+ factors for associations with sample QC and block;
#   2. measures recoverable RNA-protein geometry in the true paired cells,
#      including condition-stratified permutation controls;
#   3. measures paired RNA-protein gene and pathway-module concordance;
#   4. tests phosphopeptide condition effects in all 68 protein-QC cells,
#      before and after adjustment for the matched parent-protein abundance.
#
# Important interpretation
# ------------------------
# - The condition-stratified null preserves Control/Treated labels while
#   destroying cell correspondence. Outperforming this null therefore tests
#   more than recovery of the binary perturbation label.
# - Phosphopeptide intensity models are conditional on observing both the
#   phosphopeptide and its parent protein. Detection differences are tested
#   separately using Fisher's exact test.
# - Parent-protein adjustment supports phosphosite-specific regulation but
#   does not establish direct phosphorylation by CDK1.
# - VIM pSer55 and HNRNPU pSer247 are curated site labels reported by the NanoSPLITS study; they are not newly inferred here.
# - No protein-relative residue positions are inferred for other phosphopeptides.
#
# Run scripts 0-6 before this script.
#
# Required packages:
#   tidyverse, svglite
#
# 1. Packages and paths
required_packages <- c("tidyverse", "svglite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop("Install the following required package(s) before running this script: ",
       paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(svglite)
})

get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    context <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)

    if (!is.null(context) && !is.null(context$path) && nzchar(context$path)) {
      return(normalizePath(dirname(context$path), winslash = "/", mustWork = TRUE))
    }
  }

  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)

  if (length(file_arg) == 1) {
    script_path <- sub("^--file=", "", file_arg)
    return(normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE))
  }

  warning("Could not identify the script directory; using the working directory.")
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
processed_dir <- file.path(script_dir, "processed_data")
qc_dir <- file.path(processed_dir, "qc")
shared_dir <- file.path(processed_dir, "shared")
mofa_input_dir <- file.path(processed_dir, "mofa_input")
scot_input_dir <- file.path(processed_dir, "scot_input")
mechanistic_input_dir <- file.path(processed_dir, "mechanistic_input")
exploratory_dir <- file.path(processed_dir, "mofa_downstream", "native_paired58", "exploratory")
gsea_dir <- file.path(processed_dir, "mofa_downstream", "native_paired58", "gsea")

reference_dir <- file.path(processed_dir, "paired_reference", "native_paired58")
geometry_dir <- file.path(reference_dir, "crossmodal_recoverability")
phospho_output_dir <- file.path(reference_dir, "phosphoproteomics")

walk(c(reference_dir, geometry_dir, phospho_output_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

# Core metadata and native MOFA+ outputs
metadata_file <- file.path(processed_dir, "metadata_master_samples.csv")
factor_scores_file <- file.path(exploratory_dir, "mofa_factor_scores_oriented_long.csv")
factor_association_file <- file.path(exploratory_dir, "mofa_factor_condition_associations.csv")
variance_file <- file.path(exploratory_dir, "mofa_variance_explained_per_factor_view.csv")
weights_file <- file.path(exploratory_dir, "mofa_weights_oriented_annotated.csv")

# Frozen GSEA reference
msigdb_membership_file <- file.path(gsea_dir, "msigdb_mouse_gene_set_membership_used.csv")
gsea_all_results_file <- file.path(gsea_dir, "mofa_gsea_all_results.csv")

# Sample QC files
rna_qc_file <- file.path(qc_dir, "rna_sample_qc_metrics.csv")
protein_qc_file <- file.path(qc_dir, "protein_sample_qc_metrics_all68.csv")
phosphopeptide_qc_file <- file.path(qc_dir, "phosphopeptide_sample_qc_metrics_by_cohort.csv")

# Paired58 inputs used for recoverability calculations and the benchmark manifest
rna_pca_file <- file.path(scot_input_dir, "scot_rna_paired58_pca_embedding.csv")
protein_pca_file <- file.path(scot_input_dir, "scot_protein_paired58_pca_embedding.csv")
rna_mofa_file <- file.path(mofa_input_dir, "mofa_rna_paired58_z.csv")
protein_mofa_file <- file.path(mofa_input_dir, "mofa_protein_paired58_z_with_na.csv")
phosphopeptide_mofa_file <- file.path(mofa_input_dir, "mofa_phosphopeptide_paired58_z_with_na.csv")

# All68 inputs for parent-protein-adjusted phosphopeptide analysis
protein_all68_file <- file.path(shared_dir, "protein_log2_proda_median_normalized_all68.csv")
phosphopeptide_mechanistic_file <- file.path(
  mechanistic_input_dir,"phosphopeptide_mechanistic_all68_log2_proda_median_normalized_with_na.csv")
protein_annotation_file <- file.path(shared_dir, "protein_feature_metadata_preprocessed.csv")
phosphopeptide_annotation_file <- file.path(shared_dir, "phosphopeptide_feature_metadata_preprocessed.csv")

required_files <- c(metadata_file, factor_scores_file, factor_association_file, variance_file, weights_file,
                    msigdb_membership_file, gsea_all_results_file, rna_qc_file, protein_qc_file,
                    phosphopeptide_qc_file, rna_pca_file, protein_pca_file, rna_mofa_file, protein_mofa_file,
                    phosphopeptide_mofa_file, protein_all68_file, phosphopeptide_mechanistic_file,
                    protein_annotation_file, phosphopeptide_annotation_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Required file(s) are missing. Run scripts 0-6 first:\n", paste(missing_files, collapse = "\n"))
}

cat("Script directory:  ", script_dir, "\n")
cat("Reference outputs:", reference_dir, "\n\n")

# 2. Analysis settings
set.seed(1)

expected_n_paired <- 58L
expected_n_all68 <- 68L
condition_levels <- c("Control", "Treated")
expected_paired_condition_counts <- c(Control = 27L, Treated = 31L)

# These permutations do not retrain MOFA+. They estimate null distributions
# after preserving condition composition but destroying cell correspondence.
n_geometry_permutations <- 1000L
n_gene_permutations <- 250L
n_module_permutations <- 1000L
knn_k_values <- c(5L, 10L, 15L)

min_gene_complete_overall <- 10L
min_gene_complete_per_condition <- 5L
min_module_genes <- 5L

# At least this many cells with both phosphopeptide and parent-protein values are required
# in each condition for the intensity models.
min_parent_complete_per_condition <- 5L
phosphopeptide_fdr_threshold <- 0.05
phosphopeptide_effect_threshold <- 0.5
top_loading_fraction <- 0.10

# Independent, compact Hallmark modules used to audit cross-modal concordance.
# The exact gene membership is read from the frozen file produced by script 6.
module_pathways <- c("HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                     "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_P53_PATHWAY", 
                     "HALLMARK_XENOBIOTIC_METABOLISM", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")

target_phosphoprotein_genes <- c("Vim", "Hnrnpu")

base_family <- "Arial"
theme_set(theme_bw(base_size = 10, base_family = base_family))

# 3. Helper functions
require_columns <- function(data, required, object_name) {
  missing <- setdiff(required, colnames(data))

  if (length(missing) > 0) {
    stop(object_name, " is missing required column(s): ", paste(missing, collapse = ", "))
  }

  invisible(TRUE)
}

validate_unique_ids <- function(ids, object_name) {
  ids <- as.character(ids)

  if (anyNA(ids) || any(!nzchar(trimws(ids)))) {
    stop(object_name, " contains missing or blank identifiers.")
  }

  duplicated_ids <- unique(ids[duplicated(ids)])

  if (length(duplicated_ids) > 0) {
    stop(object_name, " contains duplicated identifiers: ", paste(head(duplicated_ids, 20), collapse = ", "))
  }

  invisible(TRUE)
}

read_feature_by_sample_matrix <- function(path, id_column, object_name) {
  data <- readr::read_csv(path, show_col_types = FALSE, na = c("", "NA"))
  require_columns(data, id_column, object_name)
  validate_unique_ids(data[[id_column]], paste(object_name, "feature IDs"))

  sample_columns <- setdiff(colnames(data), id_column)

  if (length(sample_columns) == 0) {
    stop(object_name, " contains no sample columns.")
  }

  mat <- as.matrix(data[, sample_columns, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  rownames(mat) <- as.character(data[[id_column]])
  colnames(mat) <- sample_columns

  if (any(!is.finite(mat[!is.na(mat)]))) {
    stop(object_name, " contains non-finite observed values.")
  }

  mat
}

read_sample_by_dimension_matrix <- function(path, id_column, object_name) {
  data <- readr::read_csv(path, show_col_types = FALSE)
  require_columns(data, id_column, object_name)
  validate_unique_ids(data[[id_column]], paste(object_name, "sample IDs"))

  dimension_columns <- setdiff(colnames(data), id_column)

  if (length(dimension_columns) < 2) {
    stop(object_name, " must contain at least two numeric dimensions.")
  }

  mat <- as.matrix(data[, dimension_columns, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  rownames(mat) <- as.character(data[[id_column]])

  if (anyNA(mat) || any(!is.finite(mat))) {
    stop(object_name, " contains missing or non-finite values.")
  }

  mat
}

safe_cor <- function(x, y, method = "spearman", min_complete = 3L) {
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < min_complete || sd(x) == 0 || sd(y) == 0) {
    return(NA_real_)
  }

  suppressWarnings(cor(x, y, method = method))
}

empirical_two_sided_p <- function(observed, null_values) {
  null_values <- null_values[is.finite(null_values)]

  if (!is.finite(observed) || length(null_values) == 0) {
    return(NA_real_)
  }

  (1 + sum(abs(null_values) >= abs(observed))) / (length(null_values) + 1)
}

summarise_against_null <- function(observed, null_values) {
  null_values <- null_values[is.finite(null_values)]

  tibble(Observed = observed, NullN = length(null_values),
         NullMean = ifelse(length(null_values) > 0, mean(null_values), NA_real_),
         NullSD = ifelse(length(null_values) > 1, sd(null_values), NA_real_),
         NullQ025 = ifelse(length(null_values) > 0, quantile(null_values, 0.025), NA_real_),
         NullQ500 = ifelse(length(null_values) > 0, quantile(null_values, 0.500), NA_real_),
         NullQ975 = ifelse(length(null_values) > 0, quantile(null_values, 0.975), NA_real_),
         EmpiricalP = empirical_two_sided_p(observed, null_values))
}

l2_normalize_rows <- function(mat, object_name) {
  row_norms <- sqrt(rowSums(mat^2))

  if (any(!is.finite(row_norms)) || any(row_norms <= 0)) {
    stop(object_name, " contains a sample with a zero or non-finite L2 norm.")
  }

  sweep(mat, 1, row_norms, FUN = "/")
}

make_condition_stratified_permutation <- function(condition) {
  condition <- factor(as.character(condition), levels = condition_levels)
  permutation <- seq_along(condition)

  for (condition_name in levels(droplevels(condition))) {
    index <- which(condition == condition_name)
    permutation[index] <- sample(index, length(index), replace = FALSE)
  }

  permutation
}

distance_concordance <- function(x, y) {
  safe_cor(as.vector(dist(x)), as.vector(dist(y)), method = "spearman", min_complete = 3L)
}

knn_neighbor_list <- function(mat, k) {
  if (k >= nrow(mat)) {
    stop("k must be smaller than the number of samples.")
  }

  distance_matrix <- as.matrix(dist(mat))
  diag(distance_matrix) <- Inf

  setNames(lapply(seq_len(nrow(mat)), function(i) {
      rownames(mat)[order(distance_matrix[i, ], decreasing = FALSE)[seq_len(k)]]
    }), rownames(mat))
}

knn_overlap_fraction <- function(x, y, k) {
  x_neighbors <- knn_neighbor_list(x, k)
  y_neighbors <- knn_neighbor_list(y, k)

  mean(vapply(names(x_neighbors), function(sample_id) {
    length(intersect(x_neighbors[[sample_id]], y_neighbors[[sample_id]])) / k
  }, numeric(1)))
}

collapse_feature_matrix_to_genes <- function(mat, annotation, feature_id_column, gene_column, object_name) {
  require_columns(annotation, c(feature_id_column, gene_column), object_name)

  mapping <- annotation %>%
    transmute(FeatureID = as.character(.data[[feature_id_column]]),
              Gene = trimws(as.character(.data[[gene_column]]))) %>%
    filter(FeatureID %in% rownames(mat), !is.na(Gene), nzchar(Gene)) %>% distinct(FeatureID, Gene)

  if (nrow(mapping) == 0) {
    stop("No features in ", object_name, " could be mapped to genes.")
  }

  feature_groups <- split(mapping$FeatureID, mapping$Gene)

  collapsed <- do.call(rbind, lapply(feature_groups, function(feature_ids) {
    values <- mat[feature_ids, , drop = FALSE]
    output <- colMeans(values, na.rm = TRUE)
    output[is.nan(output)] <- NA_real_
    output
  }))

  rownames(collapsed) <- names(feature_groups)
  colnames(collapsed) <- colnames(mat)
  collapsed
}

residualize_vector <- function(values, metadata, minimum_n = 8L) {
  data <- metadata %>%
    transmute(Value = values, Condition = factor(as.character(Condition), levels = condition_levels),
              Block = factor(Block))

  keep <- complete.cases(data)

  if (sum(keep) < minimum_n) {
    return(rep(NA_real_, length(values)))
  }

  model_matrix <- model.matrix(~ Condition + Block, data = data[keep, , drop = FALSE])

  if (qr(model_matrix)$rank < ncol(model_matrix)) {
    return(rep(NA_real_, length(values)))
  }

  fit <- lm(Value ~ Condition + Block, data = data[keep, , drop = FALSE])
  output <- rep(NA_real_, length(values))
  output[keep] <- residuals(fit)
  output
}

residualize_gene_matrix <- function(samples_by_genes, metadata, minimum_n) {
  output <- matrix(NA_real_, nrow = nrow(samples_by_genes), ncol = ncol(samples_by_genes),
                   dimnames = dimnames(samples_by_genes))

  for (j in seq_len(ncol(samples_by_genes))) {
    output[, j] <- residualize_vector(samples_by_genes[, j], metadata, minimum_n = minimum_n)
  }

  output
}

median_absolute_column_correlations <- function(x, y, minimum_n) {
  correlations <- vapply(seq_len(ncol(x)), function(j) {
    safe_cor(x[, j], y[, j], method = "spearman", min_complete = minimum_n)
  }, numeric(1))

  correlations <- correlations[is.finite(correlations)]

  if (length(correlations) == 0) {
    return(NA_real_)
  }

  median(abs(correlations))
}

extract_lm_coefficient <- function(fit, coefficient_name) {
  coefficient_table <- summary(fit)$coefficients

  if (!coefficient_name %in% rownames(coefficient_table)) {
    return(tibble(Estimate = NA_real_, StandardError = NA_real_, PValue = NA_real_))
  }

  tibble(Estimate = unname(coefficient_table[coefficient_name, "Estimate"]),
         StandardError = unname(coefficient_table[coefficient_name, "Std. Error"]),
         PValue = unname(coefficient_table[coefficient_name, "Pr(>|t|)"]))
}

save_svg <- function(directory, filename, plot, width, height) {
  ggsave(filename = file.path(directory, filename), plot = plot, device = svglite::svglite, width = width,
         height = height)
}

# 4. Load and validate metadata and native-reference outputs
cat("Loading paired-reference metadata and MOFA+ outputs...\n")

master_metadata <- readr::read_csv(metadata_file, show_col_types = FALSE)
factor_scores <- readr::read_csv(factor_scores_file, show_col_types = FALSE)
factor_associations <- readr::read_csv(factor_association_file, show_col_types = FALSE)
variance_per_factor <- readr::read_csv(variance_file, show_col_types = FALSE)
weights <- readr::read_csv(weights_file, show_col_types = FALSE)

require_columns(master_metadata, c("SampleID", "Condition", "ConditionLabel", "Block", "WellIndex", "ProteinQC68",
                                   "Paired_RNA_ProteinQC", "Paired_RNA_ProteinQC_Peptide"), "Master metadata")
require_columns(factor_scores, c("SampleID", "Factor", "value", "OrientedValue"), "Native factor scores")
require_columns(factor_associations, c("Factor", "TreatedMinusControl", "OrientationSign", "PValue", "FDR",
                                       "PrioritizedForInterpretation"), "Native factor associations")
require_columns(variance_per_factor, c("View", "Factor", "VarianceExplainedPercent"),
                "Native variance-explained table")
require_columns(weights, c("View", "Factor", "SourceFeatureID", "RawWeight", "OrientedWeight",
                           "PrioritizedForInterpretation"), "Native MOFA+ weights")

validate_unique_ids(master_metadata$SampleID, "Master metadata SampleID")

paired_metadata <- master_metadata %>% filter(Paired_RNA_ProteinQC_Peptide) %>%
  mutate(Condition = factor(as.character(Condition), levels = condition_levels), Block = factor(Block))

if (nrow(paired_metadata) != expected_n_paired) {
  stop("Expected ", expected_n_paired, " strict paired samples but found ", nrow(paired_metadata), ".")
}

observed_condition_counts <- table(paired_metadata$Condition)

if (!identical(as.integer(observed_condition_counts[condition_levels]), 
               as.integer(expected_paired_condition_counts[condition_levels]))) {
  stop("Unexpected condition composition in the paired58 reference.")
}

paired_ids <- paired_metadata$SampleID

if (!setequal(unique(factor_scores$SampleID), paired_ids)) {
  stop("Native factor-score SampleIDs do not match the paired58 reference.")
}

cat("Paired reference: ", nrow(paired_metadata), " cells (27 Control; 31 Treated)\n\n", sep = "")

# 5. Native factor associations with QC and block
cat("Testing native MOFA+ factors against QC and block...\n")

rna_qc <- readr::read_csv(rna_qc_file, show_col_types = FALSE) %>%
  select(SampleID, RNA_norm_sum, RNA_detected_genes_gt0, RNA_detected_genes_ge5_norm, RNA_zero_fraction)

protein_qc <- readr::read_csv(protein_qc_file, show_col_types = FALSE) %>%
  transmute(SampleID, ProteinDetected = DetectedProteins, ProteinMissingFraction = MissingFraction,
            ProteinNormalizationOffset = NormalizationOffset)

phosphopeptide_qc <- readr::read_csv(phosphopeptide_qc_file, show_col_types = FALSE) %>%
  filter(Cohort == "paired58") %>% transmute(SampleID, PhosphopeptidesDetected, PhosphopeptideMissingFraction,
                                             PhosphopeptideNormalizationOffset = NormalizationOffset)

sample_qc <- paired_metadata %>% select(SampleID, Condition, ConditionLabel, Block, WellIndex) %>%
  left_join(rna_qc, by = "SampleID") %>% left_join(protein_qc, by = "SampleID") %>%
  left_join(phosphopeptide_qc, by = "SampleID")

if (anyNA(sample_qc$Condition)) {
  stop("At least one paired58 QC row lacks condition metadata.")
}

write_csv(sample_qc, file.path(reference_dir, "paired58_combined_sample_qc.csv"))

qc_metrics <- c("RNA_norm_sum", "RNA_detected_genes_gt0", "RNA_detected_genes_ge5_norm", "RNA_zero_fraction",
                "ProteinDetected", "ProteinMissingFraction", "ProteinNormalizationOffset", "PhosphopeptidesDetected",
                "PhosphopeptideMissingFraction", "PhosphopeptideNormalizationOffset")

factor_qc_data <- factor_scores %>% select(SampleID, Factor, RawFactorScore = value) %>%
  left_join(sample_qc, by = "SampleID")

factor_qc_correlations <- expand_grid(Factor = unique(factor_qc_data$Factor), Metric = qc_metrics) %>%
  mutate(N = map2_int(Factor, Metric, function(factor_name, metric_name) {
    data <- factor_qc_data %>% filter(Factor == factor_name)
    sum(complete.cases(data$RawFactorScore, data[[metric_name]]))
    }), SpearmanRho = map2_dbl(Factor, Metric, function(factor_name, metric_name) {
      data <- factor_qc_data %>% filter(Factor == factor_name)
      safe_cor(data$RawFactorScore, data[[metric_name]], method = "spearman")
      }))

write_csv(factor_qc_correlations, file.path(reference_dir, "native_factor_qc_correlations.csv"))

factor_block_associations <- factor_qc_data %>% group_by(Factor) %>%
  group_modify(function(data, key) {
    data <- data %>%  mutate(Condition = factor(as.character(Condition), levels = condition_levels),
                             Block = factor(Block)) %>% drop_na(RawFactorScore, Condition, Block)
    reduced_fit <- lm(RawFactorScore ~ Condition, data = data)
    full_fit <- lm(RawFactorScore ~ Condition + Block, data = data)
    comparison <- anova(reduced_fit, full_fit)
    rss_reduced <- deviance(reduced_fit)
    rss_full <- deviance(full_fit)
    partial_r2 <- ifelse(is.finite(rss_reduced) && rss_reduced > 0, (rss_reduced - rss_full) / rss_reduced,
                         NA_real_)
    tibble(N = nrow(data), BlockDF = unname(comparison[2, "Df"]), ConditionAdjustedBlockPartialR2 = partial_r2,
           PValue = unname(comparison[2, "Pr(>F)"]))
    }) %>% ungroup() %>% mutate(FDR = p.adjust(PValue, method = "BH"))

write_csv(factor_block_associations, file.path(reference_dir, "native_factor_block_associations.csv"))

fit_adjusted_factor_model <- function(data, orientation_sign) {
  model_data <- data %>% transmute(RawFactorScore, Condition = factor(as.character(Condition),
                                                                      levels = condition_levels),Block = factor(Block),
                                   RNA_detected_genes_gt0, ProteinDetected, PhosphopeptidesDetected) %>% drop_na()
  if (nrow(model_data) < 20) {
    return(tibble(N = nrow(model_data), AdjustedEstimate = NA_real_, AdjustedSE = NA_real_,
                  AdjustedPValue = NA_real_, OrientedAdjustedEstimate = NA_real_,
                  ConditionPartialR2 = NA_real_, ModelAdjustedR2 = NA_real_))
  }

  full_formula <- RawFactorScore ~ Condition + Block + scale(RNA_detected_genes_gt0) + scale(ProteinDetected) +
    scale(PhosphopeptidesDetected)

  reduced_formula <- RawFactorScore ~ Block + scale(RNA_detected_genes_gt0) + scale(ProteinDetected) +
    scale(PhosphopeptidesDetected)

  full_matrix <- model.matrix(full_formula, data = model_data)
  reduced_matrix <- model.matrix(reduced_formula, data = model_data)

  if (qr(full_matrix)$rank < ncol(full_matrix) || qr(reduced_matrix)$rank < ncol(reduced_matrix)) {
    warning("Adjusted factor model is rank deficient; returning NA values.")
    return(tibble(N = nrow(model_data), AdjustedEstimate = NA_real_, AdjustedSE = NA_real_,
                  AdjustedPValue = NA_real_, OrientedAdjustedEstimate = NA_real_,
                  ConditionPartialR2 = NA_real_, ModelAdjustedR2 = NA_real_))
  }

  full_fit <- lm(full_formula, data = model_data)
  reduced_fit <- lm(reduced_formula, data = model_data)
  condition_coefficient <- extract_lm_coefficient(full_fit, "ConditionTreated")
  rss_full <- deviance(full_fit)
  rss_reduced <- deviance(reduced_fit)
  partial_r2 <- ifelse(is.finite(rss_reduced) && rss_reduced > 0, (rss_reduced - rss_full) / rss_reduced, NA_real_)

  tibble(N = nrow(model_data), AdjustedEstimate = condition_coefficient$Estimate,
         AdjustedSE = condition_coefficient$StandardError, AdjustedPValue = condition_coefficient$PValue,
         OrientedAdjustedEstimate = condition_coefficient$Estimate * orientation_sign,
         ConditionPartialR2 = partial_r2, ModelAdjustedR2 = summary(full_fit)$adj.r.squared)
}

adjusted_factor_associations <- factor_qc_data %>% group_by(Factor) %>% group_modify(function(data, key) {
  orientation_sign <- factor_associations$OrientationSign[match(key$Factor, factor_associations$Factor)]
  fit_adjusted_factor_model(data, orientation_sign)
  }) %>% ungroup() %>% left_join(factor_associations %>% select(Factor, UnadjustedEstimate = TreatedMinusControl,
                                                                UnadjustedPValue = PValue, UnadjustedFDR = FDR, OrientationSign),
                                 by = "Factor") %>% mutate(AdjustedFDR = p.adjust(AdjustedPValue, method = "BH"),
                                                           OrientedUnadjustedEstimate = UnadjustedEstimate * OrientationSign) %>%
  arrange(AdjustedPValue)

write_csv(adjusted_factor_associations, file.path(reference_dir, "native_factor_condition_associations_adjusted.csv"))

p_factor_qc <- factor_qc_correlations %>% mutate(Factor = factor(Factor, levels = unique(factor_scores$Factor)),
                                                 Metric = factor(Metric, levels = rev(qc_metrics))) %>%
  ggplot(aes(x = Factor, y = Metric, fill = SpearmanRho)) + geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", SpearmanRho)), size = 2.6) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, fill = "Spearman rho") + theme(panel.grid = element_blank())

save_svg(reference_dir, "native_factor_qc_correlations.svg", p_factor_qc, 7.5, 6.2)

# 6. Recoverable RNA-protein geometry in true paired cells
cat("Quantifying paired RNA-protein geometry...\n")

rna_pca <- read_sample_by_dimension_matrix(rna_pca_file, "SampleID", "Paired58 RNA PCA")
protein_pca <- read_sample_by_dimension_matrix(protein_pca_file, "SampleID", "Paired58 protein PCA")

if (!setequal(rownames(rna_pca), paired_ids) || !setequal(rownames(protein_pca), paired_ids)) {
  stop("Paired58 PCA SampleIDs do not match the paired reference.")
}

rna_pca <- l2_normalize_rows(rna_pca[paired_ids, , drop = FALSE], "RNA PCA")
protein_pca <- l2_normalize_rows(protein_pca[paired_ids, , drop = FALSE], "Protein PCA")

geometry_scopes <- c("All", condition_levels)
geometry_summary_list <- list()
geometry_null_list <- list()

for (scope_name in geometry_scopes) {
  scope_ids <- if (scope_name == "All") {
    paired_ids
  } else {
    paired_metadata$SampleID[paired_metadata$Condition == scope_name]
  }

  scope_metadata <- paired_metadata[match(scope_ids, paired_metadata$SampleID), , drop = FALSE]
  rna_scope <- rna_pca[scope_ids, , drop = FALSE]
  protein_scope <- protein_pca[scope_ids, , drop = FALSE]

  valid_k <- knn_k_values[knn_k_values < length(scope_ids)]
  observed_values <- c(DistanceSpearman = distance_concordance(rna_scope, protein_scope), setNames(vapply(valid_k, function(k) {
    knn_overlap_fraction(rna_scope, protein_scope, k)
    }, numeric(1)), paste0("KNNOverlap_k", valid_k)))

  null_matrix <- matrix(NA_real_, nrow = n_geometry_permutations, ncol = length(observed_values),
                        dimnames = list(NULL, names(observed_values)))

  for (permutation_index in seq_len(n_geometry_permutations)) {
    shuffled_index <- make_condition_stratified_permutation(scope_metadata$Condition)
    protein_permuted <- protein_scope[shuffled_index, , drop = FALSE]
    rownames(protein_permuted) <- scope_ids

    null_matrix[permutation_index, "DistanceSpearman"] <- distance_concordance(rna_scope, protein_permuted)

    for (k in valid_k) {
      null_matrix[permutation_index, paste0("KNNOverlap_k", k)] <- knn_overlap_fraction(rna_scope, protein_permuted, k)
    }
  }

  scope_summary <- map_dfr(names(observed_values), function(metric_name) {
    summarise_against_null(observed_values[[metric_name]], null_matrix[, metric_name]) %>%
      mutate(Scope = scope_name, Metric = metric_name, .before = 1)
  })

  scope_null <- as_tibble(null_matrix) %>% mutate(Permutation = row_number(), Scope = scope_name, .before = 1) %>%
    pivot_longer(cols = -c(Scope, Permutation), names_to = "Metric", values_to = "NullValue")

  geometry_summary_list[[scope_name]] <- scope_summary
  geometry_null_list[[scope_name]] <- scope_null
}

geometry_summary <- bind_rows(geometry_summary_list)
geometry_null <- bind_rows(geometry_null_list)

write_csv(geometry_summary, file.path(geometry_dir, "paired58_crossmodal_geometry_summary.csv"))
write_csv(geometry_null, file.path(geometry_dir, "paired58_crossmodal_geometry_permutation_null.csv"))

p_geometry <- geometry_summary %>% mutate(Metric = str_replace(Metric, "KNNOverlap_k", "kNN overlap, k="),
                                          Scope = factor(Scope, levels = geometry_scopes)) %>%
  ggplot(aes(x = Scope, y = Observed, colour = Metric)) + geom_linerange(aes(ymin = NullQ025, ymax = NullQ975),
                                                                         position = position_dodge(width = 0.55), linewidth = 0.8) +
  geom_point(position = position_dodge(width = 0.55), size = 2.5) + labs(x = NULL,
                                                                         y = "Observed value; line = stratified-null 95% interval",
                                                                         colour = NULL) + theme(legend.position = "top")

save_svg(geometry_dir, "paired58_crossmodal_geometry_vs_null.svg", p_geometry, 8.2, 4.8)

# 7. Paired RNA-protein gene and module concordance
cat("Calculating paired RNA-protein gene and module concordance...\n")

rna_mofa <- read_feature_by_sample_matrix(rna_mofa_file, "Gene", "Paired58 RNA MOFA+ matrix")
protein_mofa <- read_feature_by_sample_matrix(protein_mofa_file, "ProteinRowID", "Paired58 protein MOFA+ matrix")
phosphopeptide_mofa <- read_feature_by_sample_matrix(phosphopeptide_mofa_file, "PeptideRowID", "Paired58 phosphopeptide MOFA+ matrix")

for (matrix_name in c("rna_mofa", "protein_mofa", "phosphopeptide_mofa")) {
  matrix_object <- get(matrix_name)

  if (!setequal(colnames(matrix_object), paired_ids)) {
    stop(matrix_name, " SampleIDs do not match the paired58 reference.")
  }
}

rna_mofa <- rna_mofa[, paired_ids, drop = FALSE]
protein_mofa <- protein_mofa[, paired_ids, drop = FALSE]
phosphopeptide_mofa <- phosphopeptide_mofa[, paired_ids, drop = FALSE]

protein_annotation <- readr::read_csv(protein_annotation_file, show_col_types = FALSE)
phosphopeptide_annotation <- readr::read_csv(phosphopeptide_annotation_file, show_col_types = FALSE)

protein_gene_matrix <- collapse_feature_matrix_to_genes(protein_mofa, protein_annotation, "ProteinRowID", "Gene",
                                                        "Protein feature annotation")

common_genes <- intersect(rownames(rna_mofa), rownames(protein_gene_matrix))

if (length(common_genes) < 20) {
  stop("Fewer than 20 genes overlap the paired58 RNA and protein matrices.")
}

rna_samples_by_genes <- t(rna_mofa[common_genes, paired_ids, drop = FALSE])
protein_samples_by_genes <- t(protein_gene_matrix[common_genes, paired_ids, drop = FALSE])

gene_concordance <- map_dfr(common_genes, function(gene_name) {
  rna_values <- rna_samples_by_genes[, gene_name]
  protein_values <- protein_samples_by_genes[, gene_name]
  control_index <- paired_metadata$Condition == "Control"
  treated_index <- paired_metadata$Condition == "Treated"

  rna_residuals <- residualize_vector(rna_values, paired_metadata, minimum_n = min_gene_complete_overall)
  protein_residuals <- residualize_vector(protein_values, paired_metadata, minimum_n = min_gene_complete_overall)

  tibble(Gene = gene_name, CompleteN = sum(complete.cases(rna_values, protein_values)),
         ControlCompleteN = sum(complete.cases(rna_values[control_index], protein_values[control_index])),
         TreatedCompleteN = sum(complete.cases(rna_values[treated_index], protein_values[treated_index])),
         RhoOverall = safe_cor(rna_values, protein_values, min_complete = min_gene_complete_overall),
         RhoControl = safe_cor(rna_values[control_index], protein_values[control_index],
                               min_complete = min_gene_complete_per_condition),
         RhoTreated = safe_cor(rna_values[treated_index], protein_values[treated_index],
                               min_complete = min_gene_complete_per_condition),
         RhoResidualConditionBlock = safe_cor(rna_residuals, protein_residuals,
                                              min_complete = min_gene_complete_overall))
})

write_csv(gene_concordance, file.path(geometry_dir, "paired58_rna_protein_gene_concordance.csv"))

rna_residual_gene_matrix <- residualize_gene_matrix(rna_samples_by_genes, paired_metadata, minimum_n = min_gene_complete_overall)
protein_residual_gene_matrix <- residualize_gene_matrix(protein_samples_by_genes, paired_metadata,
                                                        minimum_n = min_gene_complete_overall)

observed_gene_metric <- median_absolute_column_correlations(rna_residual_gene_matrix, protein_residual_gene_matrix,
                                                            minimum_n = min_gene_complete_overall)

gene_null_values <- vapply(seq_len(n_gene_permutations), function(i) {
  shuffled_index <- make_condition_stratified_permutation(paired_metadata$Condition)
  protein_permuted <- protein_residual_gene_matrix[shuffled_index, , drop = FALSE]

  median_absolute_column_correlations(rna_residual_gene_matrix, protein_permuted, minimum_n = min_gene_complete_overall)
}, numeric(1))

gene_null_summary <- summarise_against_null(observed_gene_metric, gene_null_values) %>%
  mutate(Metric = "MedianAbsoluteConditionBlockResidualGeneRho", CommonGenes = length(common_genes), .before = 1)

write_csv(gene_null_summary, file.path(geometry_dir, "paired58_gene_concordance_null_summary.csv"))
write_csv(tibble(Permutation = seq_along(gene_null_values), NullValue = gene_null_values),
          file.path(geometry_dir, "paired58_gene_concordance_permutation_null.csv"))

# Hallmark module scores
msigdb_membership <- readr::read_csv(msigdb_membership_file, show_col_types = FALSE)
require_columns(msigdb_membership, c("Library", "Pathway", "Gene"), "Frozen MSigDB membership")

available_modules <- intersect(module_pathways, unique(msigdb_membership$Pathway))
missing_modules <- setdiff(module_pathways, available_modules)

if (length(missing_modules) > 0) {
  warning("The following requested Hallmark modules are absent from the frozen MSigDB file: ",
          paste(missing_modules, collapse = ", "))
}

calculate_module_score <- function(samples_by_genes, pathway_genes) {
  genes_used <- intersect(colnames(samples_by_genes), pathway_genes)

  if (length(genes_used) < min_module_genes) {
    return(list(score = rep(NA_real_, nrow(samples_by_genes)), genes = genes_used))
  }

  score <- rowMeans(samples_by_genes[, genes_used, drop = FALSE], na.rm = TRUE)
  score[is.nan(score)] <- NA_real_
  list(score = score, genes = genes_used)
}

module_score_rows <- list()
module_correlation_rows <- list()
module_null_rows <- list()

for (pathway_name in available_modules) {
  pathway_genes <- msigdb_membership %>% filter(Pathway == pathway_name) %>% pull(Gene) %>% unique()

  rna_module <- calculate_module_score(rna_samples_by_genes, pathway_genes)
  protein_module <- calculate_module_score(protein_samples_by_genes, pathway_genes)

  if (length(rna_module$genes) < min_module_genes || length(protein_module$genes) < min_module_genes) {
    warning("Skipping module ", pathway_name, " because fewer than ", min_module_genes,
            " genes are represented in at least one view.")
    next
  }

  module_data <- paired_metadata %>% select(SampleID, Condition, Block) %>%
    mutate(Pathway = pathway_name, RNAModuleScore = rna_module$score, ProteinModuleScore = protein_module$score,
           RNAGenesUsed = length(rna_module$genes), ProteinGenesUsed = length(protein_module$genes))

  rna_residual <- residualize_vector(module_data$RNAModuleScore, module_data, minimum_n = 10L)
  protein_residual <- residualize_vector(module_data$ProteinModuleScore, module_data, minimum_n = 10L)

  observed_residual_rho <- safe_cor(rna_residual, protein_residual, min_complete = 10L)

  module_null <- vapply(seq_len(n_module_permutations), function(i) {
    shuffled_index <- make_condition_stratified_permutation(module_data$Condition)
    safe_cor(rna_residual, protein_residual[shuffled_index], min_complete = 10L)
  }, numeric(1))

  control_index <- module_data$Condition == "Control"
  treated_index <- module_data$Condition == "Treated"

  module_score_rows[[pathway_name]] <- module_data
  module_correlation_rows[[pathway_name]] <- summarise_against_null(observed_residual_rho, module_null) %>%
    mutate(Pathway = pathway_name, RNAGenesUsed = length(rna_module$genes), ProteinGenesUsed = length(protein_module$genes),
           RhoOverall = safe_cor(module_data$RNAModuleScore, module_data$ProteinModuleScore, min_complete = 10L),
      RhoControl = safe_cor(module_data$RNAModuleScore[control_index], module_data$ProteinModuleScore[control_index],
                            min_complete = 5L), RhoTreated = safe_cor(module_data$RNAModuleScore[treated_index],
                                                                      module_data$ProteinModuleScore[treated_index],
                                                                      min_complete = 5L), .before = 1)

  module_null_rows[[pathway_name]] <- tibble(Pathway = pathway_name, Permutation = seq_along(module_null),
                                             NullResidualRho = module_null)
}

module_scores <- bind_rows(module_score_rows)
module_correlations <- bind_rows(module_correlation_rows)
module_null <- bind_rows(module_null_rows)

write_csv(module_scores, file.path(geometry_dir, "paired58_rna_protein_module_scores.csv"))

# Correct empirical module-level P values across all tested modules.
module_correlations <- module_correlations %>%
  mutate(EmpiricalFDR = p.adjust(EmpiricalP, method = "BH"), PrespecifiedG2MBenchmark = Pathway == "HALLMARK_G2M_CHECKPOINT",
         EvidenceCategory = case_when(EmpiricalFDR < 0.05 ~ "FDR_significant",
                                      PrespecifiedG2MBenchmark & EmpiricalP < 0.05 ~ "Prespecified_G2M_nominal",
                                      TRUE ~ "Not_supported"))

write_csv(module_correlations, file.path(geometry_dir, "paired58_rna_protein_module_concordance.csv"))
write_csv(module_null, file.path(geometry_dir, "paired58_module_concordance_permutation_null.csv"))

if (nrow(module_correlations) > 0) {
  p_modules <- module_correlations %>% mutate(PathwayLabel = Pathway %>% str_remove("^HALLMARK_") %>%
                                                str_replace_all("_", " ") %>% str_to_sentence(), 
                                              PathwayLabel = fct_reorder(PathwayLabel, Observed)) %>%
    ggplot(aes(x = Observed, y = PathwayLabel)) + geom_vline(xintercept = 0, colour = "grey65", linetype = 2) +
    geom_errorbarh(aes(xmin = NullQ025, xmax = NullQ975), height = 0.16, colour = "grey45") +
    geom_point(aes(colour = EvidenceCategory), size = 2.5) +
    scale_colour_manual(values = c(FDR_significant = "#0072B2", Prespecified_G2M_nominal = "#D55E00", Not_supported = "grey35"),
                        name = NULL) +
    labs(x = "Condition/block-residual RNA-protein Spearman rho\n(line = stratified-null 95% interval)", y = NULL)

  save_svg(geometry_dir, "paired58_module_concordance_vs_null.svg", p_modules, 8.0, 5.4)
}

# 8. Parent-protein-adjusted phosphopeptide analysis in all68
cat("Running parent-protein-adjusted phosphopeptide analysis...\n")

all68_metadata <- master_metadata %>% filter(ProteinQC68) %>%
  mutate(Condition = factor(as.character(Condition), levels = condition_levels), Block = factor(Block))

if (nrow(all68_metadata) != expected_n_all68) {
  stop("Expected 68 protein-QC cells for phosphopeptide validation.")
}

all68_ids <- all68_metadata$SampleID

protein_all68 <- read_feature_by_sample_matrix(protein_all68_file, "ProteinRowID", "All68 normalized protein matrix")
phosphopeptide_mechanistic <- read_feature_by_sample_matrix(phosphopeptide_mechanistic_file, "PeptideRowID",
                                                            "All68 mechanistic phosphopeptide matrix")

if (!setequal(colnames(protein_all68), all68_ids) || !setequal(colnames(phosphopeptide_mechanistic), all68_ids)) {
  stop("All68 protein/phosphopeptide SampleIDs do not match the metadata.")
}

protein_all68 <- protein_all68[, all68_ids, drop = FALSE]
phosphopeptide_mechanistic <- phosphopeptide_mechanistic[, all68_ids, drop = FALSE]

require_columns(protein_annotation, c("ProteinRowID", "Entry.Name", "Gene"), "Protein feature annotation")
require_columns(phosphopeptide_annotation, c("PeptideRowID", "Modified.Sequence", "Assigned.Modifications",
                                             "Entry.Name", "Gene", "MechanisticTestEligible_all68"),
                "Phosphopeptide feature annotation")

protein_parent_keys <- protein_annotation %>% transmute(ParentProteinRowID = as.character(ProteinRowID),
                                                        ParentEntryName = trimws(as.character(Entry.Name)),
                                                        ParentGene = trimws(as.character(Gene))) %>%
  filter(!is.na(ParentEntryName), nzchar(ParentEntryName))
if (anyDuplicated(protein_parent_keys$ParentEntryName)) {
  stop("Protein Entry.Name is not unique; parent-protein mapping is ambiguous.")
}

phosphopeptide_parent_map <- phosphopeptide_annotation %>%
  transmute(PeptideRowID = as.character(PeptideRowID), ModifiedSequence = as.character(Modified.Sequence), 
            AssignedModifications = as.character(Assigned.Modifications), PhosphopeptideEntryName = trimws(as.character(Entry.Name)),
            Gene = trimws(as.character(Gene)), MechanisticTestEligible_all68) %>%
  left_join(protein_parent_keys, by = c("PhosphopeptideEntryName" = "ParentEntryName")) %>%
  mutate(PresentInMechanisticMatrix = PeptideRowID %in% rownames(phosphopeptide_mechanistic),
         ParentProteinMapped = !is.na(ParentProteinRowID),
         MappingStatus = case_when(!PresentInMechanisticMatrix ~ "Not_mechanistic_test_eligible",
                                   !ParentProteinMapped ~ "Parent_protein_not_mapped", TRUE ~ "Mapped_and_test_eligible"))

write_csv(phosphopeptide_parent_map, file.path(phospho_output_dir, "phosphopeptide_parent_mapping_qc.csv"))

testable_parent_map <- phosphopeptide_parent_map %>% filter(PresentInMechanisticMatrix, ParentProteinMapped)

if (nrow(testable_parent_map) == 0) {
  stop("No mechanistic phosphopeptides could be mapped to parent proteins.")
}

fit_one_phosphopeptide <- function(peptide_row_id, parent_protein_row_id) {
  phosphopeptide_values <- phosphopeptide_mechanistic[peptide_row_id, all68_ids]
  parent_values <- protein_all68[parent_protein_row_id, all68_ids]

  model_data <- all68_metadata %>% select(SampleID, Condition, Block) %>%
    mutate(PhosphopeptideIntensity = as.numeric(phosphopeptide_values[SampleID]),
           ParentProteinIntensity = as.numeric(parent_values[SampleID]), Detected = !is.na(PhosphopeptideIntensity))

  control_detected <- sum(model_data$Detected[model_data$Condition == "Control"], na.rm = TRUE)
  treated_detected <- sum(model_data$Detected[model_data$Condition == "Treated"], na.rm = TRUE)
  control_total <- sum(model_data$Condition == "Control")
  treated_total <- sum(model_data$Condition == "Treated")

  detection_table <- matrix(c(control_detected, control_total - control_detected, treated_detected,
                              treated_total - treated_detected), nrow = 2, byrow = TRUE,
                            dimnames = list(condition_levels, c("Detected", "Missing")))

  fisher_result <- fisher.test(detection_table)

  complete_data <- model_data %>% filter(!is.na(PhosphopeptideIntensity), !is.na(ParentProteinIntensity)) %>%
    mutate(Condition = factor(as.character(Condition), levels = condition_levels), Block = factor(Block),
           ParentProteinZ = as.numeric(scale(ParentProteinIntensity)))

  control_complete <- sum(complete_data$Condition == "Control")
  treated_complete <- sum(complete_data$Condition == "Treated")
  intensity_test_eligible <- control_complete >= min_parent_complete_per_condition &&
    treated_complete >= min_parent_complete_per_condition && is.finite(sd(complete_data$ParentProteinIntensity)) &&
    sd(complete_data$ParentProteinIntensity) > 0
  empty_intensity_result <- tibble(CompleteN = nrow(complete_data), ControlCompleteN = control_complete,
                                   TreatedCompleteN = treated_complete, IntensityTestEligible = intensity_test_eligible,
                                   UnadjustedConditionEffect = NA_real_, UnadjustedConditionSE = NA_real_,
                                   UnadjustedConditionPValue = NA_real_, AdjustedConditionEffect = NA_real_,
                                   AdjustedConditionSE = NA_real_, AdjustedConditionPValue = NA_real_,
                                   ParentProteinEffectPerSD = NA_real_, ParentProteinPValue = NA_real_, AdjustedModelR2 = NA_real_,
                                   AdjustedModelAdjustedR2 = NA_real_)

  if (!intensity_test_eligible) {
    intensity_result <- empty_intensity_result
  } else {
    unadjusted_formula <- PhosphopeptideIntensity ~ Condition + Block
    adjusted_formula <- PhosphopeptideIntensity ~ Condition + ParentProteinZ + Block

    unadjusted_matrix <- model.matrix(unadjusted_formula, data = complete_data)
    adjusted_matrix <- model.matrix(adjusted_formula, data = complete_data)

    if (qr(unadjusted_matrix)$rank < ncol(unadjusted_matrix) || qr(adjusted_matrix)$rank < ncol(adjusted_matrix) ||
        nrow(complete_data) - ncol(adjusted_matrix) < 3) {
      intensity_result <- empty_intensity_result %>% mutate(IntensityTestEligible = FALSE)
      } else {
        unadjusted_fit <- lm(unadjusted_formula, data = complete_data)
        adjusted_fit <- lm(adjusted_formula, data = complete_data)
        unadjusted_condition <- extract_lm_coefficient(unadjusted_fit, "ConditionTreated")
        adjusted_condition <- extract_lm_coefficient(adjusted_fit, "ConditionTreated")
        parent_effect <- extract_lm_coefficient(adjusted_fit, "ParentProteinZ")
        
        intensity_result <- tibble(CompleteN = nrow(complete_data), ControlCompleteN = control_complete,
                                   TreatedCompleteN = treated_complete, IntensityTestEligible = TRUE,
                                   UnadjustedConditionEffect = unadjusted_condition$Estimate,
                                   UnadjustedConditionSE = unadjusted_condition$StandardError,
                                   UnadjustedConditionPValue = unadjusted_condition$PValue,
                                   AdjustedConditionEffect = adjusted_condition$Estimate,
                                   AdjustedConditionSE = adjusted_condition$StandardError,
                                   AdjustedConditionPValue = adjusted_condition$PValue,
                                   ParentProteinEffectPerSD = parent_effect$Estimate, ParentProteinPValue = parent_effect$PValue,
                                   AdjustedModelR2 = summary(adjusted_fit)$r.squared,
                                   AdjustedModelAdjustedR2 = summary(adjusted_fit)$adj.r.squared)
    }
  }

  bind_cols(tibble(PeptideRowID = peptide_row_id, ParentProteinRowID = parent_protein_row_id,
                   ControlDetectedN = control_detected, TreatedDetectedN = treated_detected,
                   ControlDetectedFraction = control_detected / control_total, 
                   TreatedDetectedFraction = treated_detected / treated_total,
                   DetectionFractionDifference = treated_detected / treated_total - control_detected / control_total,
                   DetectionOddsRatio = unname(fisher_result$estimate), DetectionPValue = fisher_result$p.value), intensity_result)
}

phosphopeptide_tests <- map2_dfr(testable_parent_map$PeptideRowID, testable_parent_map$ParentProteinRowID, fit_one_phosphopeptide) %>%
  left_join(testable_parent_map %>% select(PeptideRowID, ModifiedSequence, AssignedModifications, PhosphopeptideEntryName, Gene,
                                           ParentProteinRowID, ParentGene), by = c("PeptideRowID", "ParentProteinRowID")) %>%
  mutate(DetectionFDR = p.adjust(DetectionPValue, method = "BH"),
         UnadjustedConditionFDR = p.adjust(UnadjustedConditionPValue, method = "BH"),
         AdjustedConditionFDR = p.adjust(AdjustedConditionPValue, method = "BH"),
         ParentProteinFDR = p.adjust(ParentProteinPValue, method = "BH"),
         AdjustedDirection = case_when(AdjustedConditionEffect > 0 ~ "Treated-associated",
                                       AdjustedConditionEffect < 0 ~ "Control-associated",
                                       TRUE ~ NA_character_),
         ParentAdjustedPriority = !is.na(AdjustedConditionFDR) &  AdjustedConditionFDR <= phosphopeptide_fdr_threshold & 
           abs(AdjustedConditionEffect) >= phosphopeptide_effect_threshold) %>%
  arrange(AdjustedConditionFDR, desc(abs(AdjustedConditionEffect)))

write_csv(phosphopeptide_tests, file.path(phospho_output_dir, "phosphopeptide_parent_adjusted_condition_results_all68.csv"), na = "")

# Integrate MOFA+ loadings with parent-adjusted phosphopeptide evidence
phosphopeptide_mofa_mechanistic <- weights %>% filter(View == "Phosphopeptide", PrioritizedForInterpretation) %>%
  transmute(Factor, PeptideRowID = SourceFeatureID, RawWeight, OrientedWeight) %>% group_by(Factor) %>%
  mutate(AbsoluteWeight = abs(OrientedWeight), AbsoluteWeightRank = min_rank(desc(AbsoluteWeight)),
         AbsoluteWeightPercentile = percent_rank(AbsoluteWeight)) %>% ungroup() %>%
  left_join(phosphopeptide_tests, by = "PeptideRowID") %>%
  mutate(WeightEffectDirectionConcordant = case_when(is.na(AdjustedConditionEffect) ~ NA,
                                                     OrientedWeight == 0 | AdjustedConditionEffect == 0 ~ NA,
                                                     TRUE ~ sign(OrientedWeight) == sign(AdjustedConditionEffect)),
         TopLoading = AbsoluteWeightPercentile >= (1 - top_loading_fraction),
         JointMechanisticPriority = ParentAdjustedPriority & TopLoading & WeightEffectDirectionConcordant) %>%
  arrange(Factor, desc(JointMechanisticPriority), desc(AbsoluteWeight))

write_csv(phosphopeptide_mofa_mechanistic, file.path(phospho_output_dir, "phosphopeptide_MOFA_parent_adjusted_integration.csv"),
          na = "")

# Explicit lookup for phosphoproteins highlighted in the NanoSPLITS paper.
# An absent row is retained to distinguish "not present" from "not significant".
targeted_phosphoprotein_lookup <- tibble(TargetGene = target_phosphoprotein_genes) %>%
  left_join(phosphopeptide_parent_map %>%  rename(TargetGene = Gene), by = "TargetGene") %>%
  mutate(PresentInSourcePhosphopeptideTable = !is.na(PeptideRowID)) %>%
  left_join(phosphopeptide_tests %>% select(PeptideRowID, ParentProteinRowID, ControlDetectedN, TreatedDetectedN,
                                            DetectionFractionDifference, DetectionFDR, CompleteN, ControlCompleteN,
                                            TreatedCompleteN, IntensityTestEligible, UnadjustedConditionEffect,
                                            UnadjustedConditionFDR, AdjustedConditionEffect, AdjustedConditionFDR,
                                            ParentProteinEffectPerSD, ParentProteinFDR, ParentAdjustedPriority),
            by = c("PeptideRowID", "ParentProteinRowID")) %>% arrange(TargetGene, PeptideRowID)

# Curated mapping of targeted phosphopeptides to sites reported in the NanoSPLITS study.
target_site_annotations <- tribble(~PeptideRowID, ~CanonicalUniProtAccession, ~CanonicalEntryName, ~NanoSPLITSReportedSite,
                                   ~SiteAnnotationStatus, "PEP_023136", "P20152", "VIME_MOUSE", "pSer55",
                                   "Matched_to_NanoSPLITS_reported_site", "PEP_000999", "Q8VEK3", "HNRPU_MOUSE", "pSer247",
                                   "Matched_to_NanoSPLITS_reported_site")

targeted_phosphoprotein_lookup <- targeted_phosphoprotein_lookup %>% left_join(target_site_annotations, by = "PeptideRowID")

# Ensure that both curated target peptides were recovered.
expected_target_peptides <- c("PEP_023136", "PEP_000999")

if (!all(expected_target_peptides %in% targeted_phosphoprotein_lookup$PeptideRowID)) {
  missing_target_peptides <- setdiff(expected_target_peptides, targeted_phosphoprotein_lookup$PeptideRowID)
  
  stop("Curated VIM/HNRNPU phosphopeptide(s) not found: ", paste(missing_target_peptides, collapse = ", "))
}

write_csv(targeted_phosphoprotein_lookup, file.path(phospho_output_dir, "targeted_Vim_Hnrnpu_phosphopeptide_lookup.csv"), na = "")


phosphopeptide_plot_data <- phosphopeptide_tests %>% filter(IntensityTestEligible, is.finite(AdjustedConditionEffect)) %>%
  mutate(PlotFDR = pmax(AdjustedConditionFDR, .Machine$double.xmin),
         Label = if_else(ParentAdjustedPriority | min_rank(AdjustedConditionFDR) <= 8, paste0(Gene, " | ", PeptideRowID), NA_character_))

if (nrow(phosphopeptide_plot_data) > 0) {
  p_phosphopeptide <- ggplot(phosphopeptide_plot_data, aes(x = AdjustedConditionEffect, y = -log10(PlotFDR))) +
    geom_vline(xintercept = c(-phosphopeptide_effect_threshold, phosphopeptide_effect_threshold),
               colour = "grey70", linetype = 2) + geom_hline(yintercept = -log10(phosphopeptide_fdr_threshold),
                                                             colour = "grey70", linetype = 2) +
    geom_point(aes(colour = ParentAdjustedPriority), alpha = 0.8, size = 2) +
    geom_text(aes(label = Label), check_overlap = TRUE, size = 2.5, vjust = -0.6) +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey45"),
                        labels = c(`TRUE` = "Prioritized", `FALSE` = "Not prioritized")) +
    labs(x = "Parent-protein-adjusted Treated - Control effect\n(log2 median-normalized phosphopeptide intensity)",
         y = expression(-log[10](FDR)), colour = NULL) + theme(legend.position = "top")

  save_svg(phospho_output_dir, "phosphopeptide_parent_adjusted_volcano.svg", p_phosphopeptide, 8.2, 5.8)
}

phosphopeptide_loading_plot_data <- phosphopeptide_mofa_mechanistic %>%
  filter(IntensityTestEligible, is.finite(OrientedWeight), is.finite(AdjustedConditionEffect))

if (nrow(phosphopeptide_loading_plot_data) > 0) {
  p_loading_effect <- ggplot(phosphopeptide_loading_plot_data, aes(x = OrientedWeight, y = AdjustedConditionEffect)) +
    geom_hline(yintercept = 0, colour = "grey75") + geom_vline(xintercept = 0, colour = "grey75") +
    geom_point(aes(colour = WeightEffectDirectionConcordant), alpha = 0.8, size = 2) +
    facet_wrap(~ Factor, scales = "free_x") + scale_colour_manual(values = c(`TRUE` = "#009E73", `FALSE` = "#CC79A7"),
                                                                  na.value = "grey55") +
    labs(x = "Condition-oriented MOFA+ phosphopeptide weight", y = "Parent-adjusted condition effect",
         colour = "Direction agrees") + theme(legend.position = "top")

  save_svg(phospho_output_dir, "phosphopeptide_MOFA_weight_vs_parent_adjusted_effect.svg", p_loading_effect, 8.5, 4.8)
}

# 9. Freeze the paired-reference benchmark
cat("Writing paired-reference manifests and recovery targets...\n")

paired_sample_manifest <- paired_metadata %>% mutate(SampleOrder = match(SampleID, paired_ids), .before = 1) %>% arrange(SampleOrder)

write_csv(paired_sample_manifest, file.path(reference_dir, "paired58_sample_manifest.csv"))

feature_manifest <- bind_rows(tibble(View = "RNA", FeatureID = rownames(rna_mofa)),
                              tibble(View = "Protein", FeatureID = rownames(protein_mofa)),
                              tibble(View = "Phosphopeptide", FeatureID = rownames(phosphopeptide_mofa))) %>% group_by(View) %>%
  mutate(FeatureOrder = row_number()) %>% ungroup() %>% select(View, FeatureOrder, FeatureID)

write_csv(feature_manifest, file.path(reference_dir, "paired58_feature_manifest.csv"))

input_file_manifest <- tibble(File = required_files) %>%
  mutate(RelativeFile = if_else(startsWith(File, paste0(script_dir, .Platform$file.sep)), substring(File, nchar(script_dir) + 2L),
                                File), Bytes = file.info(File)$size, MD5 = unname(tools::md5sum(File))) %>%
  select(RelativeFile, Bytes, MD5)

write_csv(input_file_manifest, file.path(reference_dir, "paired58_reference_input_file_manifest.csv"))

native_factor_reference_summary <- factor_associations %>% 
  left_join(adjusted_factor_associations, by = "Factor", suffix = c("", "_AdjustedTable")) %>% 
  left_join(variance_per_factor %>% select(Factor, View, VarianceExplainedPercent) %>%
              pivot_wider(names_from = View, values_from = VarianceExplainedPercent,
                          names_glue = "VarianceExplained_{View}_percent"), by = "Factor") %>% arrange(PValue)

write_csv(native_factor_reference_summary, file.path(reference_dir, "native_factor_reference_summary.csv"), na = "")

recovery_targets <- tribble(~TargetID, ~TargetType, ~NativeDefinition, ~PrimaryAlignedComparison,
                            "T01", "Latent_factor", "Dominant condition-associated native factor",
                            "Raw factor-score correlation after factor matching", "T02", "View_contribution",
                            "Native per-factor variance explained in RNA, Protein and Phosphopeptide",
                            "Matched-factor per-view variance-explained concordance", "T03", "Feature_loading",
                            "Complete native raw and condition-oriented loading ranks",
                            "Spearman rank correlation and top-feature overlap on identical features", "T04",
                            "Pathway", "Complete native Hallmark and Reactome NES profiles",
                            "NES correlation and direction agreement using frozen gene sets", "T05",
                            "Crossmodal_geometry", "True-pair RNA-protein distance and kNN concordance",
                            "SCOT+ improvement over unaligned and condition-stratified null", "T06",
                            "Within_condition",
                            "Control- and Treated-specific geometry plus condition/block-residual gene and module concordance",
                            "Recovery after preserving condition but hiding cell correspondence", "T07", "Phosphoproteomics",
                            "Native phosphopeptide loading plus parent-protein-adjusted all68 condition effect",
                            "Retention of loading direction and mechanistic candidate ranking", "T08", "QC_robustness",
                            "Native factor condition effect before and after QC/block adjustment",
                            "Comparable QC diagnostics in the aligned model")

write_csv(recovery_targets, file.path(reference_dir, "paired58_reference_recovery_targets.csv"))

parameter_report <- tribble(~Parameter, ~Value, "random_seed", "1", "n_geometry_permutations", as.character(n_geometry_permutations),
                            "n_gene_permutations", as.character(n_gene_permutations), "n_module_permutations",
                            as.character(n_module_permutations), "permutation_strategy", "condition_stratified",
                            "PCA_row_normalization", "L2", "knn_k_values", paste(knn_k_values, collapse = ";"),
                            "minimum_gene_complete_overall", as.character(min_gene_complete_overall),
                            "minimum_gene_complete_per_condition", as.character(min_gene_complete_per_condition),
                            "minimum_module_genes", as.character(min_module_genes), "minimum_parent_complete_per_condition",
                            as.character(min_parent_complete_per_condition), "phosphopeptide_FDR_threshold",
                            as.character(phosphopeptide_fdr_threshold), "phosphopeptide_effect_threshold",
                            as.character(phosphopeptide_effect_threshold), "phosphopeptide_intensity_model",
                            "phosphopeptide ~ condition + parent_protein_z + block", "phosphopeptide_detection_test",
                            "Fisher_exact_Control_vs_Treated")

write_csv(parameter_report, file.path(reference_dir, "paired58_reference_diagnostic_parameters.csv"))

analysis_notes <- c("This script does not retrain or modify the native MOFA+ model.",
                    "Condition-stratified permutations preserve Control/Treated labels while destroying cell correspondence.",
                    "Factor matching to a future aligned model must use raw factor scores and absolute correlation; factor numbers and signs are arbitrary.",
                    "Treatment-based sign orientation is only stable for condition-associated factors.",
                    "The aligned model must use the same paired58 SampleIDs and feature manifests for the controlled benchmark.",
                    "The frozen MSigDB membership file from script 6 must be reused for aligned-model GSEA.",
                    "Phosphopeptide intensity models are conditional on observed phosphopeptide and parent-protein values.",
                    "Detection differences are tested separately and should be interpreted together with conditional intensity effects.",
                    "Parent-protein-adjusted association does not establish direct CDK1 phosphorylation.",
                    "VIM pSer55 and HNRNPU pSer247 are curated NanoSPLITS-reported site labels; no other protein-relative positions are inferred.")

write_lines(analysis_notes, file.path(reference_dir, "paired58_reference_analysis_notes.txt"))
write_lines(capture.output(sessionInfo()), file.path(reference_dir, "paired58_reference_sessionInfo.txt"))

# 10. Final checks and completion message
stopifnot(nrow(paired_sample_manifest) == expected_n_paired, setequal(feature_manifest$View, c("RNA", "Protein", "Phosphopeptide")),
          nrow(geometry_summary) > 0, nrow(gene_concordance) > 0, nrow(phosphopeptide_tests) > 0,
          all(is.finite(geometry_summary$Observed)), all(is.finite(phosphopeptide_tests$DetectionPValue)))

# Freeze paired-reference benchmarks before SCOT+ alignment
paired_reference_benchmark_manifest <- tribble(~BenchmarkID, ~Priority, ~Cohort, ~ReferenceTarget, ~PlannedAlignedModelEvaluation,
                                               ~InterpretationRule, "F1_treatment_program", "Primary", "Paired58_all",
                                               "Native Factor 1",
                                               paste("Match the aligned factor to native Factor 1 using absolute",
                                                     "factor-score and weight correlations; report oriented score",
                                                     "correlation, per-view R2, top weights and pathway agreement."),
                                               paste("Recovery requires concordant treatment direction and", 
                                                     "multi-view biological interpretation."), "Treated_geometry_k5", "Primary",
                                               "Treated_only", "RNA-protein k=5 neighbourhood correspondence",
                                               paste("Compare SCOT+ alignment with true pairing and a",
                                                     "condition-stratified permutation null."),
                                               "Improvement should be assessed relative to the stratified null.",
                                               "Treated_geometry_k10", "Primary", "Treated_only",
                                               "RNA-protein k=10 neighbourhood correspondence",
                                               paste("Compare SCOT+ alignment with true pairing and a",
                                                     "condition-stratified permutation null."),
                                               "Improvement should be assessed relative to the stratified null.",
                                               "G2M_module_concordance", "Secondary_prespecified", "Paired58_all_and_treated",
                                               "HALLMARK_G2M_CHECKPOINT RNA-protein concordance",
                                               paste("Compare module-score concordance between native and", "aligned models."),
                                               paste("Treat as a prespecified biological benchmark;", 
                                                     "report both empirical P value and BH FDR."), "Gene_residual_concordance",
                                               "Secondary", "Paired58_all",
                                               paste("Condition- and block-adjusted RNA-protein", "gene-level concordance"),
                                               paste("Compare the median absolute residual correlation with", 
                                                     "the same condition-stratified null."),
                                               "Interpret as weak distributed correspondence, not strong gene coupling.",
                                               "Factor1_phosphopeptide_program", "Mechanistic_validation", "Paired58", 
                                               "Factor 1 phosphopeptide weights and 11 joint-priority phosphopeptides",
                                               paste("Test whether the aligned MOFA factor preserves loading", 
                                                     "directions and prioritization."),
                                               "Phosphoproteomics validates the aligned factor but is not directly SCOT+-aligned.",
                                               "VIM_HNRNPU_targets", "Targeted_validation", "All68_phosphoproteomics",
                                               "NanoSPLITS-reported VIM pSer55 and HNRNPU pSer247 phosphopeptides",
                                               paste("Confirm that the aligned-model-associated phosphoproteomic", 
                                                     "factor retains the expected directions."), 
                                               "Evaluate these as targeted external biological validation.", "Control_geometry",
                                               "Negative_control", "Control_only", "Native control-cell RNA-protein geometry",
                                               "Report SCOT+ results but do not require above-null recovery.",
                                               "The native paired reference showed no detectable control-cell geometry.",
                                               "Factor3_technical", "Negative_control", "Paired58_all", "Native Factor 3",
                                               "Assess correlation with RNA detection metrics in the aligned model.",
                                               "Do not interpret recovery of Factor 3 as biological success.")
write_csv(paired_reference_benchmark_manifest, file.path(reference_dir, "paired_reference_benchmark_manifest.csv"))
cat("\nPaired-reference diagnostics complete.\n\n")
cat("Primary reference directory:\n  ", reference_dir, "\n\n", sep = "")
cat("Key outputs:\n")
cat("  native_factor_condition_associations_adjusted.csv\n")
cat("  crossmodal_recoverability/paired58_crossmodal_geometry_summary.csv\n")
cat("  crossmodal_recoverability/paired58_rna_protein_gene_concordance.csv\n")
cat("  crossmodal_recoverability/paired58_rna_protein_module_concordance.csv\n")
cat("  phosphoproteomics/phosphopeptide_parent_adjusted_condition_results_all68.csv\n")
cat("  phosphoproteomics/phosphopeptide_MOFA_parent_adjusted_integration.csv\n")
cat("  paired58_reference_recovery_targets.csv\n")