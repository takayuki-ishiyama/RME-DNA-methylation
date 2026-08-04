# Code availability and scope

## Recovered code

The following analysis components were recovered directly from the retained R workspace:

- `run_primary_model()`
- `run_dmp_one_group()`
- `run_allterms()`
- `run_oras()`
- `step_filter()`
- `run_champ_single()`
- `run_champ_on_clean()`

These functions were used as the basis for the scripts in this repository.

## Verified downstream analysis

The retained batch-adjusted M-value matrix reproduced all 164 published CpGs using:

- `limma::duplicateCorrelation`
- `limma::lmFit`
- `limma::eBayes`
- Benjamini-Hochberg FDR correction
- the manuscript-specified covariates

## Reconstructed preprocessing

The manuscript states that preprocessing included:

1. out-of-band background correction;
2. RELIC dye-bias correction;
3. RCP correction of Type I/II probe-design bias;
4. quantile normalization;
5. ChAMP-based filtering;
6. selective sex-chromosome filtering;
7. batch adjustment using `limma::removeBatchEffect`.

The exact original ENmix command sequence and exact `removeBatchEffect` call were not preserved in the command history. These parts are therefore documented as reconstructed rather than exact historical code.

## Public release recommendation

Do not describe `01_preprocessing_reconstructed.R` as the verbatim original script. It is a transparent reconstruction based on the Methods and retained objects. `02_statistical_analysis_verified.R` is the strongest verified component.
