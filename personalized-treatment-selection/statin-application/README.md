# Population-Based Observational Study of Statin Use — R Application

## Scope

This folder contains the author's original R script for the statin-use application associated with the paper *Personalized Treatment Selection Using Observational Data*. The script is preserved **as supplied**; it has not been rewritten, corrected, or independently validated as executable end-to-end.

## Files

- `R/R-code.R` — original research script, including estimation functions and statin application code.

## Methods represented in the script

- Adaptive LASSO variable selection using ridge-derived penalty factors and cross-validation (`xnames`).
- Gradient-boosting outcome modeling and relative-influence-based variable selection (`var.gbm`).
- Multinomial treatment propensity modeling and propensity-score weighting, including a GBM-based approach through `twang` (`match.gbm`, `match.lasso`).
- Weighted Nadaraya–Watson and local-linear estimators (`alpha.hat.w`).
- Treatment-specific score construction and comparison of treatment assignments.
- A five-fold split stratified by the three observed treatment groups, followed by training/test construction.

## Data

The original script refers to `statin.csv` and cites a public Dryad data source:

- https://doi.org/10.5061/dryad.6qs53

**The dataset is not included in this repository.** Users must independently obtain the appropriate data and confirm the variable definitions and any applicable reuse conditions. The original script expects a CSV with variables including `Ethnicity`, `Education`, `Albuminuria`, `z`, and `y`; it also depends on the source dataset's column ordering.

## R packages

The script calls `MASS`, `nnet`, `glmnet`, `twang`, `gbm`, `tibble`, `dplyr`, `kedd`, and `fastDummies`. Package availability and API compatibility may vary by R/package version.

## Running the research script

1. Obtain the appropriate dataset from the cited source and review its data dictionary.
2. Install the packages listed above.
3. Update the original local path in the line beginning `statin <- read.csv(...)` to point to your local `statin.csv` file.
4. Inspect and resolve the execution dependencies noted below before running the analysis.

## Important execution and reproducibility notes

- The original script contains a hard-coded Windows data path; it is **not portable without editing that path**.
- Several functions refer to `m` and `yn` from their calling environment rather than receiving them as explicit arguments. The supplied wrapper `resultgbm` defines these locally, but `match.gbm` is called as a separate function, so the original script may fail unless the scope is addressed.
- The script constructs five train/test folds but calls `resultgbm(test1, train1)` only; it does not automatically aggregate results over all five folds.
- After constructing a dummy-encoded dataset, the script immediately assigns `dat <- statin`, overriding the preceding constructed dataset. This should be checked against the original analysis intent.
- Treatment categories, data types, column ordering, and model prediction inputs should be checked against the original data before execution.
- The script includes experimental/intermediate model objects and code that has not been independently checked for compatibility with current package versions.
- No numeric findings or exact reproduction of the published application are claimed from this release.

## Research provenance

This repository makes an existing research script available publicly. The GitHub upload date is the date of public code release, **not** the date the research was originally performed. Cite the associated journal publication when using the methodology or discussing the findings.
