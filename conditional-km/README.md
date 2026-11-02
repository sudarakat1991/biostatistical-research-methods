# Biostatistical Research Methods — Conditional Kaplan–Meier Estimator

This repository contains R research scripts and saved simulation results associated with doctoral work on the **conditional Kaplan–Meier estimator with functional covariates**. These are original research scripts and research outputs, **not** an installable R package or a fully validated, one-command reproduction pipeline.

**Related publication:** Tholkage, S., Zheng, Q., & Kulasekera, K. B. (2022). *Conditional Kaplan–Meier Estimator with Functional Covariates for Time-to-Event Data*. *Stats*, 5(4), 1113–1129.

## Repository contents

```text
conditional-km/
├── simulations/
│   ├── Scenario1/    simu.R, true1.txt, brier1.txt
│   ├── Scenario2/    simu.R, true1.txt, brier1.txt
│   ├── Scenario3/    simu.R, true1.txt, brier1.txt
│   └── Scenario4/    simu.R, true1.txt, brier1.txt
├── simulation-results/
│   ├── NewMeansdwt.R
│   └── figures/Results.docx
└── application/
    └── Brier_comp.R
```

## Simulation workflow

Each `simulations/Scenario*/simu.R` contains scenario-specific simulation and evaluation code, including simulated functional covariates and censored time-to-event outcomes, conditional survival estimation, bandwidth evaluation, and comparisons with conventional Kaplan–Meier estimation. Each script writes two output files in its **current R working directory**:

- `true1.txt` — saved results associated with bandwidth selection using error against known survival probabilities in the simulation setting.
- `brier1.txt` — saved results associated with Brier-score-based bandwidth selection.

**All four scenarios have their own script and saved output pair.** The updated `Scenario2/simu.R` specifies sample sizes **n = 100 and n = 400**, consistent with the other supplied scenario scripts. The scenarios differ in their simulation models and/or tuning grids; do not interchange their scripts or output files.

The saved output files are included to permit inspection without rerunning the computationally intensive simulations. The simulations use parallel workers (8 or 10, depending on the script) and may require substantial computing time. No complete rerun or numerical equivalence check has been performed for this release.

### Summary analysis and figures

`conditional-km/simulation-results/NewMeansdwt.R` reads each scenario's `true1.txt` and `brier1.txt` and performs comparative summaries and plotting. Its original `setwd()` commands point to local Windows directories from the research environment, including an earlier name for Scenario 2 (`Scenario2_re`). To run it from this repository, **update each `setwd()` path** to the corresponding directory under `conditional-km/simulations/` (or replace them with repository-relative paths). The original script has been preserved without algorithmic edits.

`conditional-km/simulation-results/figures/Results.docx` contains saved figures for the four scenarios, including comparisons at **n = 100** and **n = 400** with approximately **15% and 25% censoring**, bandwidth-selection plots, and additional figures (including n = 1000 and an application comparison). These figures are research records; the supplied scripts have **not** been verified to recreate every figure in the document.

## Real-data application

`conditional-km/application/Brier_comp.R` contains code for an application involving functional covariates and comparisons of Brier scores for conditional survival estimation and other survival-analysis approaches. It references the data files `hippocampussurface.dat` and `clinical.dat`, which are available through supporting information for the following publication:

- [Biometrics article (DOI: 10.1111/biom.12748)](https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.12748)
- [Publisher's supporting-information ZIP (S2)](https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fbiom.12748&file=biom12748-sup-0002-SuppInfo-S2.zip)

These external datasets are **not redistributed** here. Download them from the publisher, review any applicable terms, and update the original local data path in `Brier_comp.R` before attempting to run the application.

## Software dependencies and execution notes

The R scripts use packages including `survival`, `robustHD`, `MASS`, and `parallel` for simulations; `ggplot2`, `tidyr`, `cowplot`, `plyr`, `ggpubr`, and `gridExtra` for results analysis; and `locpol`, `survAUC`, `survcomp`, and `DescTools` for parts of the application workflow. Consult the individual scripts for exact package calls and dependencies.

1. Open the desired `simulations/Scenario*/simu.R` and inspect the parameter settings and parallel-worker count before running it.
2. Set the R working directory to that scenario's folder if you want its generated `true1.txt` and `brier1.txt` saved there. **Back up the supplied output files first**, as rerunning the script may overwrite them.
3. Update the original Windows-specific paths in `NewMeansdwt.R` to load the four saved output pairs from their respective scenario folders.
4. For the application, obtain the supporting data from the publisher and update `Brier_comp.R` to point to the extracted files.

**Reproducibility status:** The repository documents the code and saved results, but a clean-environment end-to-end execution, exact reproduction of all plots, and software validation have not been completed. Do not treat the scripts as production-ready or clinically validated software.

## Provenance and public release

The research scripts and results are associated with the doctoral research and 2022 publication cited above. The repository organization and this README were prepared for later public dissemination. GitHub commit dates show when files were uploaded, not necessarily when the underlying research was originally performed.

## Licensing and use

No software license has been selected pending confirmation of appropriate permissions from any relevant coauthors and institutions. This code is provided for scholarly inspection and is not intended for clinical decision-making without independent validation. Do not include confidential patient data, restricted third-party data, or proprietary employer materials.
