# =============================================================================
# 01_preprocessing_reconstructed.R
# =============================================================================
# Reconstructed preprocessing workflow for:
# DNA Methylation Changes Following Rapid Maxillary Expansion in Children
#
# IMPORTANT:
# The original command history for the exact ENmix preprocessing calls was not
# preserved. This script combines the manuscript Methods with recovered helper
# functions. Verify the ENmix block against retained intermediate matrices.
# =============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c(
  "minfi",
  "ENmix",
  "ChAMP",
  "SummarizedExperiment",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop("Install missing packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(minfi)
  library(ENmix)
  library(ChAMP)
  library(SummarizedExperiment)
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
})

# -----------------------------------------------------------------------------
# 1. Input files
# -----------------------------------------------------------------------------
idat_dir <- "PATH/TO/IDAT_FILES"
metadata_file <- "sample_metadata.csv"
output_dir <- "processed_data"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pheno <- read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)

required_metadata <- c(
  "sentrix_key",
  "sample_id",
  "subject_id",
  "Group",
  "Time",
  "treatment_period",
  "device",
  "sex",
  "age"
)

missing_metadata <- setdiff(required_metadata, names(pheno))
if (length(missing_metadata) > 0L) {
  stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
}

# -----------------------------------------------------------------------------
# 2. IDAT import
# -----------------------------------------------------------------------------
targets <- data.frame(
  Sample_Name = pheno$sample_id,
  Sentrix_ID = sub("_.*$", "", pheno$sentrix_key),
  Sentrix_Position = sub("^.*_", "", pheno$sentrix_key),
  Basename = file.path(idat_dir, pheno$sentrix_key),
  stringsAsFactors = FALSE
)

rg <- minfi::read.metharray.exp(
  targets = targets,
  extended = TRUE,
  force = TRUE
)
colnames(rg) <- pheno$sentrix_key

# Retained workspace reference:
# class(rg) = RGChannelSetExtended
# dim(rg)   = 1,105,209 x 48

mset_raw <- minfi::preprocessRaw(rg)

# Retained workspace reference:
# class(mset_raw) = MethylSet
# dim(mset_raw)   = 936,990 x 48

# -----------------------------------------------------------------------------
# 3. ENmix preprocessing
# -----------------------------------------------------------------------------
# Manuscript-described steps:
# - out-of-band background correction
# - RELIC dye-bias correction
# - RCP probe-design correction
# - quantile normalization
#
# The exact historical function calls were not retained. Replace this block
# only after verifying the commands against retained objects:
#   mset_enmix : 936,990 x 48
#   mset_rcp   : 936,990 x 48
#   beta_start : 936,990 x 48
#
# Example placeholder:
# mset_enmix <- ...
# mset_rcp   <- ...
# beta_start <- ...
#
stop(
  "The ENmix block must be verified against retained intermediate objects ",
  "before executing the remainder of this script."
)

# -----------------------------------------------------------------------------
# 4. Detection P values and bead counts
# -----------------------------------------------------------------------------
detP <- minfi::detectionP(rg)
beads <- SummarizedExperiment::assay(rg, "NBeads")

# -----------------------------------------------------------------------------
# 5. Recovered filtering helpers
# -----------------------------------------------------------------------------
clean_cg <- function(x) {
  sub("^((cg|rs)[0-9]+).*", "\\1", trimws(as.character(x)))
}

step_filter <- function(
    beta_in,
    doDetP = FALSE,
    doBeads = FALSE,
    doNoCG = FALSE,
    doSNP = FALSE,
    doMulti = FALSE
) {
  cf <- ChAMP::champ.filter(
    beta = beta_in,
    pd = as.data.frame(SummarizedExperiment::colData(rg)),
    detP = if (doDetP) detP else NULL,
    beadcount = if (doBeads) beads else NULL,
    autoimpute = FALSE,
    filterXY = FALSE,
    filterDetP = doDetP,
    filterBeads = doBeads,
    filterNoCG = doNoCG,
    filterSNPs = doSNP,
    filterMultiHit = doMulti,
    arraytype = "EPICv2",
    detPcut = 0.01,
    beadCutoff = 0.05,
    detSamplecut = 1,
    ProbeCutoff = 0,
    SampleCutoff = 0
  )
  cf$beta
}

run_champ_single <- function(
    beta_in,
    doNoCG = FALSE,
    doSNP = FALSE,
    doMulti = FALSE
) {
  if (nrow(beta_in) == 0L) return(beta_in)

  cf <- ChAMP::champ.filter(
    beta = beta_in,
    pd = as.data.frame(SummarizedExperiment::colData(rg)),
    detP = NULL,
    beadcount = NULL,
    autoimpute = FALSE,
    filterXY = FALSE,
    filterDetP = FALSE,
    filterBeads = FALSE,
    filterNoCG = doNoCG,
    filterSNPs = doSNP,
    filterMultiHit = doMulti,
    arraytype = "EPICv2",
    ProbeCutoff = 0,
    SampleCutoff = 0
  )
  cf$beta
}

run_champ_on_clean <- function(
    beta_in,
    doNoCG = FALSE,
    doSNP = FALSE,
    doMulti = FALSE
) {
  if (nrow(beta_in) == 0L) return(beta_in)

  orig_ids <- rownames(beta_in)
  clean_ids <- clean_cg(orig_ids)

  if (anyDuplicated(clean_ids)) {
    stop(
      "Probe IDs are duplicated after suffix removal. Apply the same ",
      "representative-probe selection used in the retained analysis."
    )
  }

  beta_tmp <- beta_in
  rownames(beta_tmp) <- clean_ids

  cf <- ChAMP::champ.filter(
    beta = beta_tmp,
    pd = as.data.frame(SummarizedExperiment::colData(rg)),
    detP = NULL,
    beadcount = NULL,
    autoimpute = FALSE,
    filterXY = FALSE,
    filterDetP = FALSE,
    filterBeads = FALSE,
    filterNoCG = doNoCG,
    filterSNPs = doSNP,
    filterMultiHit = doMulti,
    arraytype = "EPICv2",
    ProbeCutoff = 0,
    SampleCutoff = 0
  )

  keep_clean <- rownames(cf$beta)
  keep_orig <- orig_ids[match(keep_clean, clean_ids)]
  beta_in[keep_orig, , drop = FALSE]
}

# -----------------------------------------------------------------------------
# 6. Reconstructed staged filtering
# -----------------------------------------------------------------------------
# beta_detP <- step_filter(beta_start, doDetP = TRUE)
# beta_bead <- step_filter(beta_detP, doBeads = TRUE)
#
# Resolve duplicated cleaned EPIC v2 probe IDs using the same representative
# selection rule used in the retained analysis.
#
# beta_noCG  <- run_champ_on_clean(beta_unique, doNoCG = TRUE)
# beta_snp   <- run_champ_on_clean(beta_noCG, doSNP = TRUE)
# beta_multi <- run_champ_on_clean(beta_snp, doMulti = TRUE)
#
# Selective sex-chromosome filtering should then be applied using the
# EPIC v2 hg38 annotation.
#
# Expected final matrix:
# 895,978 CpGs x 48 samples

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo_preprocessing.txt")
)
