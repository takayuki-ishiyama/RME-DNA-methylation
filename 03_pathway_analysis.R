# =============================================================================
# 03_pathway_analysis.R
# =============================================================================
# KEGG over-representation analysis using missMethyl.
# =============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c("missMethyl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install missing packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(missMethyl))

sig_cpg <- treatment_result$significant$CpG
all_tested_cpg <- treatment_result$table$CpG

kegg_result <- missMethyl::gometh(
  sig.cpg = sig_cpg,
  all.cpg = all_tested_cpg,
  collection = "KEGG",
  array.type = "EPIC",
  prior.prob = TRUE
)

kegg_result$FDR <- p.adjust(kegg_result$P.DE, method = "BH")
kegg_result <- kegg_result[
  order(kegg_result$FDR, kegg_result$P.DE),
  ,
  drop = FALSE
]

output_dir <- "analysis_results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  kegg_result,
  file.path(output_dir, "KEGG_ORA_all_pathways.csv"),
  row.names = FALSE
)

write.csv(
  head(kegg_result, 20),
  file.path(output_dir, "KEGG_ORA_top20_pathways.csv"),
  row.names = FALSE
)

print(head(kegg_result, 3))

# Expected top-ranked pathways:
# 1. Autophagy - animal
# 2. Oxidative phosphorylation
# 3. Cell cycle
#
# No KEGG pathway reached FDR < 0.05.
#
# KEGG imagery requires formal permission from Kanehisa Laboratories.
