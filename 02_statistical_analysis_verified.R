# =============================================================================
# 02_statistical_analysis_verified.R
# =============================================================================
# Verified downstream statistical analysis.
# This model reproduced all 164 published treatment-group CpGs when applied
# to the retained batch-adjusted M-value matrix.
# =============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c("limma")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install missing packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(limma))

# Required objects:
# beta_before : 895,978 x 48 normalized/QC-filtered beta matrix
# M_before    : 895,978 x 48 normalized, unadjusted M-value matrix
# M_after     : 895,978 x 48 batch-adjusted M-value matrix used in primary analysis
# pheno       : metadata with sentrix_key, sample_id, subject_id, Group, Time,
#               treatment_period, device, sex, age, Epithelial, Immune

required_pheno <- c(
  "sentrix_key", "sample_id", "subject_id", "Group", "Time",
  "treatment_period", "device", "sex", "age", "Epithelial", "Immune"
)
missing_pheno <- setdiff(required_pheno, names(pheno))
if (length(missing_pheno) > 0L) {
  stop("Missing phenotype columns: ", paste(missing_pheno, collapse = ", "))
}

pheno <- pheno[match(colnames(M_after), pheno$sentrix_key), , drop = FALSE]

stopifnot(
  identical(pheno$sentrix_key, colnames(M_after)),
  nrow(M_after) == 895978L,
  ncol(M_after) == 48L
)

pheno$Group <- factor(pheno$Group, levels = c("Treatment", "Control"))
pheno$Time <- factor(pheno$Time, levels = c("Pre", "Post"))
pheno$sex <- factor(pheno$sex, levels = c("F", "M"))
pheno$device <- factor(pheno$device, levels = c("No", "Yes"))
pheno$age <- as.numeric(pheno$age)
pheno$treatment_period <- as.numeric(pheno$treatment_period)

run_primary_model <- function(M_matrix, pheno_group, design, sample_index) {
  M_group <- M_matrix[, sample_index, drop = FALSE]

  corfit <- limma::duplicateCorrelation(
    M_group,
    design,
    block = pheno_group$subject_id
  )

  fit <- limma::lmFit(
    M_group,
    design,
    block = pheno_group$subject_id,
    correlation = corfit$consensus
  )

  fit <- limma::eBayes(fit)

  table_all <- limma::topTable(
    fit,
    coef = "TimePost",
    number = Inf,
    sort.by = "P"
  )

  table_all$CpG <- rownames(table_all)
  table_all$FDR <- p.adjust(table_all$P.Value, method = "BH")

  list(
    fit = fit,
    corfit = corfit,
    table = table_all,
    significant = table_all[table_all$FDR < 0.05, , drop = FALSE]
  )
}

run_dmp_one_group <- function(
    group_label,
    phenotype,
    M_matrix,
    beta_matrix
) {
  sample_index <- which(phenotype$Group == group_label)
  if (length(sample_index) < 3L) {
    stop("Too few samples for group: ", group_label)
  }

  pdg <- droplevels(phenotype[sample_index, , drop = FALSE])
  beta_group <- beta_matrix[, sample_index, drop = FALSE]

  if (group_label == "Treatment") {
    design <- model.matrix(
      ~ Time + age + sex + treatment_period + device +
        Epithelial + Immune,
      data = pdg
    )
  } else if (group_label == "Control") {
    design <- model.matrix(
      ~ Time + age + Epithelial + Immune,
      data = pdg
    )
  } else {
    stop("group_label must be Treatment or Control")
  }

  result <- run_primary_model(
    M_matrix = M_matrix,
    pheno_group = pdg,
    design = design,
    sample_index = sample_index
  )

  pre_ix <- which(pdg$Time == "Pre")
  post_ix <- which(pdg$Time == "Post")

  mean_pre <- rowMeans(beta_group[, pre_ix, drop = FALSE], na.rm = TRUE)
  mean_post <- rowMeans(beta_group[, post_ix, drop = FALSE], na.rm = TRUE)

  result$table$MeanBeta_Pre <- mean_pre[result$table$CpG]
  result$table$MeanBeta_Post <- mean_post[result$table$CpG]
  result$table$DeltaBeta <- (
    result$table$MeanBeta_Post - result$table$MeanBeta_Pre
  )

  result$significant <- result$table[
    result$table$FDR < 0.05,
    ,
    drop = FALSE
  ]

  result$group <- group_label
  result$design <- design
  result
}

treatment_result <- run_dmp_one_group(
  group_label = "Treatment",
  phenotype = pheno,
  M_matrix = M_after,
  beta_matrix = beta_before
)

control_result <- run_dmp_one_group(
  group_label = "Control",
  phenotype = pheno,
  M_matrix = M_after,
  beta_matrix = beta_before
)

# Absolute DeltaBeta thresholds
delta_threshold_summary <- data.frame(
  threshold = c(0.01, 0.02, 0.05),
  n_CpGs = vapply(
    c(0.01, 0.02, 0.05),
    function(x) {
      sum(
        abs(treatment_result$significant$DeltaBeta) > x,
        na.rm = TRUE
      )
    },
    integer(1)
  )
)

# Alternative batch-adjustment sensitivity analysis
run_position_covariate_sensitivity <- function(
    M_unadjusted,
    phenotype
) {
  sample_index <- which(phenotype$Group == "Treatment")
  pdg <- droplevels(phenotype[sample_index, , drop = FALSE])
  M_group <- M_unadjusted[, sample_index, drop = FALSE]

  pdg$slide_position <- factor(sub("^.*_", "", pdg$sentrix_key))

  design <- model.matrix(
    ~ Time + age + sex + treatment_period + device +
      Epithelial + Immune + slide_position,
    data = pdg
  )

  corfit <- limma::duplicateCorrelation(
    M_group,
    design,
    block = pdg$subject_id
  )

  fit <- limma::lmFit(
    M_group,
    design,
    block = pdg$subject_id,
    correlation = corfit$consensus
  )
  fit <- limma::eBayes(fit)

  table_all <- limma::topTable(
    fit,
    coef = "TimePost",
    number = Inf,
    sort.by = "P"
  )
  table_all$CpG <- rownames(table_all)
  table_all$FDR <- p.adjust(table_all$P.Value, method = "BH")

  list(
    fit = fit,
    corfit = corfit,
    design = design,
    table = table_all,
    significant = table_all[table_all$FDR < 0.05, , drop = FALSE]
  )
}

position_sensitivity <- run_position_covariate_sensitivity(
  M_unadjusted = M_before,
  phenotype = pheno
)

# Validation against published results
validation <- c(
  treatment_164 = nrow(treatment_result$significant) == 164L,
  control_0 = nrow(control_result$significant) == 0L,
  hyper_73 = sum(
    treatment_result$significant$DeltaBeta > 0,
    na.rm = TRUE
  ) == 73L,
  hypo_91 = sum(
    treatment_result$significant$DeltaBeta < 0,
    na.rm = TRUE
  ) == 91L,
  delta_001_55 = delta_threshold_summary$n_CpGs[1] == 55L,
  delta_002_21 = delta_threshold_summary$n_CpGs[2] == 21L,
  delta_005_1 = delta_threshold_summary$n_CpGs[3] == 1L,
  sensitivity_3 = nrow(position_sensitivity$significant) == 3L
)

print(validation)

output_dir <- "analysis_results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  treatment_result$table,
  file.path(output_dir, "Treatment_all_CpGs.csv"),
  row.names = FALSE
)
write.csv(
  treatment_result$significant,
  file.path(output_dir, "Treatment_164_CpGs.csv"),
  row.names = FALSE
)
write.csv(
  control_result$table,
  file.path(output_dir, "Control_all_CpGs.csv"),
  row.names = FALSE
)
write.csv(
  position_sensitivity$significant,
  file.path(output_dir, "Sensitivity_position_3_CpGs.csv"),
  row.names = FALSE
)
write.csv(
  delta_threshold_summary,
  file.path(output_dir, "DeltaBeta_threshold_summary.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo_analysis.txt")
)
