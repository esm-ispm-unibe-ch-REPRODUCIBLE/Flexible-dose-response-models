
# Flexible dose-response models

This repository contains the code and analysis files for estimating
flexible dose-response models using patient-visit data.

## Description

The analysis focuses on treatment dose assignment over follow-up visits.
The main clinical outcome is the Hamilton Depression Rating Scale (HAMD)
score, and treatment response is defined as improvement from baseline.

The project includes models for:

  - dose-assignment probabilities
  - stabilized inverse probability of treatment weights
  - stabilized inverse probability of censoring weights
  - total stabilized weights
  - marginal structural model analyses

## Repository structure

``` text
.
├── Outputs/                  # Rendered reports and output files
├── R/                        # R functions used in the analysis
├── R-scripts/                # Main analysis scripts
├── .gitattributes
├── Flexible dose-response.Rproj
└── README.md
```

## Requirements

The analysis uses R and the following packages:

``` r
install.packages(c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "knitr",
  "MASS",
  "scales",
  "forcats",
  "scales",
  "MASS",
  "rms",
  "purr"
))
```

## How to run the analysis

Open the project file in RStudio and run the main analysis script:

``` r
source("R-scripts/00_master.R")
```

Alternatively, run the individual scripts in the order indicated in the
`R-scripts/` folder.

## Main outputs

The analysis produces:

  - descriptive summaries of patient-visits, dose levels, efficacy, and
    side effects
  - dose-assignment model coefficient tables
  - stabilized treatment-weight summaries
  - stabilized censoring-weight summaries
  - total-weight summaries before and after truncation
  - diagnostic plots for dose assignment and weight stability

## Notes

The data are not available in this repository.

## Authors

Konstantina Chalkou
