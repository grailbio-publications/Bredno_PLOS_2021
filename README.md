## Overview

This project provides data and code to reproduce the analysis results
presented in Bredno J et al., "Clinical determinants of cfDNA tumor fraction from
studies of a multi-cancer early detection test", 2021.

### Background

Liquid biopsy applications of cell-free DNA (cfDNA) are often limited by the amount of circulating tumor DNA (ctDNA) and the fraction of cfDNA in a blood sample that is derived from tumor cells (tumor fraction). This tumor fraction varies widely between individuals and cancer types. Clinical factors that influence tumor fraction have not been completely elucidated.

### Methods and Findings

Tumor fraction was determined for breast, lung, and colorectal cancer participant samples in the first substudy of the Circulating Cell-free Genome Atlas study CCGA (NCT02889978), and was related to clinical tumor and patient characteristics. Linear models were created to determine the influence of tumor size combined with mitotic or metabolic activity (as tumor mitotic volume or excessive lesion glycolysis, respectively), histologic type, histologic grade, and lymph node status on tumor fraction. For breast cancer and lung cancer, tumor mitotic volume and excessive lesion glycolysis (primary lesion volume scaled by %Ki-67+ or PET SUV - 1.0, respectively) were the only statistically significant covariates. For colorectal cancer, the surface area of tumors invading beyond the subserosa was the only significant covariate. The models were validated with cases from the second CCGA substudy, and show that these clinical determinants of tumor fraction can predict and explain the performance of a multi-cancer early detection test.

### Conclusions

Prognostic clinical variables, including mitotic or metabolic activity and depth of invasion, were identified as determinants of ctDNA by linear models that relate clinical covariates to tumor fraction. The identified correlations indicate that faster growing tumors have higher tumor fractions. Early cancer detection from liquid biopsy is determined by tumor fraction. The results therefore support that early detection is particularly sensitive for faster growing, aggressive tumors with high mortality, many of which have no available screening today.

## System requirements and installation

Follow instructions at https://cran.r-project.org/bin/linux/ubuntu/README.html
to install R (version 3.6.2 was used for this analysis). The packages needed to
reproduce the results have been captured with the package `packrat`.

To add the required packages, open R in the project directory:

```
packrat::restore()
```

This pulls the versioned dependencies from CRAN and Bioconductor.

## Data Access

All data required to reproduce the results is provided in the directory
`data`. The data used in the presented analyses are provided in tables separate
for breast, lung, and colorectal cancers and for participants used in training and validation. 

- `supplemental_data_bc_training.tsv`: Training data to generate the biophysical
  model for breast cancer.
- `supplemental_data_bc_validation.tsv`: Validation data for the breast cancer
  model.
- `supplemental_data_crc_training.tsv`: Training data to generate the biophysical
  model for colorectal cancer.
- `supplemental_data_crc_validation.tsv`: Validation data for the colorectal cancer
  model.
- `supplemental_data_luc_training.tsv`: Training data to generate the biophysical
  model for lung cancer.
- `supplemental_data_luc_validation.tsv`: Validation data for the lung cancer
  model.
  
The columns and variables are described in `data_dictionary.txt`.

## Code

The results of biophysical model generation and validation are generated
separately for breast, lung, and colorectal cancer using resp. markdowns in the
folder `Rmd`. Shared code and functionality is contained in
`R/biophysical_modeling.R`, which is sourced in each of the markdowns.

For example, to generate tables, lists and figures for the lung cancer model,
render the document `report_luc_models.Rmd` in the directory `Rmd`:

```
setwd("Rmd")
rmarkdown::render("report_luc_model.Rmd")
```

This generates tables and lists in the file `report_luc_model.html` and also
creates figures in eps and svg formats in the directory `figures` and saves
models for prediction in the directory `cmd`.

A shell script allows to predict cfDNA tumor fraction and total number
of cell-free tumor-derived genome equivalents (GE) in the circulation of a patient.
The prediction for all models is available in one shell script:

```
cd cmd
./predict_ectf.R
```

This shell script queries user input to select a cancer type, then enter number
of primary lesions, the size for each lesion, and cancer-type-specific clinical
information to estimate cfDNA tumor fraction and total amount of ctDNA.
