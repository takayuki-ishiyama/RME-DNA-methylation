# DNA Methylation Changes Following Rapid Maxillary Expansion in Children

This repository contains the R code associated with the manuscript:

**DNA Methylation Changes Following Rapid Maxillary Expansion in Children**

## Data availability

Raw IDAT files and processed DNA-methylation beta values are deposited in NCBI GEO:

**GSE317740**

The GEO record must be publicly accessible before final acceptance.

## Repository structure

- `01_preprocessing_reconstructed.R`  
  Reconstructed preprocessing and quality-control workflow based on the manuscript Methods, retained R objects, and recovered ChAMP helper functions.

- `02_statistical_analysis_verified.R`  
  Verified downstream limma workflow. When applied to the retained batch-adjusted M-value matrix, this model reproduced all 164 published treatment-group CpGs.

- `03_pathway_analysis.R`  
  KEGG over-representation analysis using `missMethyl`.

- `sample_metadata_template.csv`  
  Required metadata structure.

- `CODE_AVAILABILITY.md`  
  Scope and limitations of the recovered code.

- `editor_response_comment1.txt`  
  Draft response to the Scientific Reports editor.

## Software versions

- R 4.5.1
- minfi 1.54.1
- ENmix 1.44.3
- ChAMP 2.38.1
- limma 3.64.3
- EpiDISH 2.24.0
- missMethyl 1.42.0

## Statistical models

### Treatment group

The model was:

```r
~ Time + age + sex + treatment_period + device + Epithelial + Immune
```

Reference levels:

- Time: `Pre`
- Sex: `F`
- Device: `No`

Repeated samples were handled using `limma::duplicateCorrelation` with `subject_id` as the blocking variable.

### Control group

The model was:

```r
~ Time + age + Epithelial + Immune
```

Sex was excluded because all control participants were female.

## Expected published results

- 895,978 analyzed CpGs
- 164 FDR-significant CpGs in the treatment group
- 0 FDR-significant CpGs in the control group
- 73 hypermethylated CpGs
- 91 hypomethylated CpGs
- Absolute Δβ thresholds:
  - >0.01: 55 CpGs
  - >0.02: 21 CpGs
  - >0.05: 1 CpG
- Alternative position-covariate sensitivity analysis: 3 CpGs
- Device interaction analysis: 13 CpGs
- No KEGG pathway significant at FDR < 0.05

## Important limitation

The original command history for the exact ENmix preprocessing calls and the exact `removeBatchEffect` call was not preserved. Therefore, the preprocessing script is explicitly labeled as reconstructed. The downstream limma model was verified against the retained analysis objects and reproduced the published 164-CpG result.

## KEGG permission

KEGG pathway imagery is copyrighted by Kanehisa Laboratories. Formal permission is required for publication in *Scientific Reports* under an Open Access license.

## License

MIT License.
