# Biostatistical Research Methods

This repository shares selected R research code developed in connection with doctoral research in biostatistics. The projects address survival analysis with complex covariates and individualized treatment selection using observational data.

The repository is organized by research project. Each project has its own README describing the available scripts, data requirements, and scope of the released materials.

## Projects

### 1. Conditional Kaplan–Meier estimation

**Folder:** [`conditional-km/`](conditional-km/)

Research code associated with a conditional Kaplan–Meier estimator for survival analysis with functional covariates. The project includes simulation scripts, simulation-results processing, and an application involving longitudinal or functional predictors. See the project README for the original study details, file descriptions, and data dependencies.

### 2. Personalized treatment selection using observational data

**Folder:** [`personalized-treatment-selection/`](personalized-treatment-selection/)

Research code associated with personalized treatment selection using observational data, including machine-learning approaches and propensity-score weighting. **The materials currently shared focus on the population-based observational study of statin use**; they are not presented as a complete release of all simulations or analyses from the paper.

The application code includes gradient boosting, adaptive LASSO, weighting, and nonparametric estimation components. Consult the project README for its data source, software dependencies, and execution considerations.

## Repository structure

```text
biostatistical-research-methods/
├── README.md
├── conditional-km/
│   └── ...
└── personalized-treatment-selection/
    └── ...
```

## Scope and reproducibility

- The code is shared for research transparency and methodological reference.
- Scripts may require project-specific data files, package versions, working-directory adjustments, or other setup described in their respective READMEs.
- Availability of source code does not imply that every analysis can be run end-to-end without additional data or configuration.
- Original research code and later repository documentation should be distinguished; the GitHub publication date is not the original date of the research.

## Publications and attribution

These materials relate to the author's doctoral research publications on conditional survival estimation and personalized treatment selection. For complete bibliographic information, methodological assumptions, and interpretation of results, consult the corresponding published papers and project-specific documentation.

## Use and contact

Please cite the corresponding publication when referring to the research methods. No software license is asserted in this README; reuse permissions depend on any license separately provided in the repository.
