# W_Kenya_ANC_fitting

Using ANC data to reconstruct community prevalence trends \# Assessing
Malaria Prevalence: Community vs Antenatal Care in Rarieda Sub-county,
Western Kenya

This project analyzes the relationship between malaria prevalence within
antenatal care (ANC) and the community in Rarieda sub-county, Western
Kenya. The analysis includes generating synthetic data, fitting models
of differing complexity, making out-of-sample predictions, and
quantifying the incremental value of ANC data.

You can view the full report [here](https://patrickgtwalker.github.io/W_Kenya_ANC_fitting/malaria_analysis.html).
## Table of Contents

-   [Introduction](#introduction)
-   [Generating Synthetic Data](#generating-synthetic-data)
-   [Fitting Models](#fitting-models)
-   [Model Evaluation](#model-evaluation)
-   [Out-of-Sample Predictions](#out-of-sample-predictions)
-   [Quantitative Assessment](#quantitative-assessment)
-   [Incremental Value of ANC Data](#incremental-value-of-anc-data)
-   [Usage](#usage)
-   [Repository Structure](#repository-structure)
-   [Contributing](#contributing)
-   [License](#license)

## Introduction 

This analysis summarizes the key features of our approach to assessing
the relationship between malaria prevalence within ANC relative to the
community in Rarieda sub-county, Western Kenya. We will first show how
we fit and compare models of differing complexity. Then, we will
demonstrate how to generate out-of-sample predictions for settings based
on ANC prevalence alone and quantify the incremental value of the data
compared to information typically available via population-based
surveys.

## Generating Synthetic Data 

We generate synthetic data representing the simplest and most complex
models to provide toy examples which then allow us to illustrate our approach.

## Fitting Models 

We fit the generated data to models of varying complexity to compare
their performance.

## Model Evaluation 

We evaluate the fitted models by comparing the parameters to the
simulated values and assess their predictive accuracy.

## Out-of-Sample Predictions

We make out-of-sample predictions to evaluate how well the models
reconstruct ANC prevalence between surveys spaced over a few months
every few years.

## Quantitative Assessment

We use two metrics to assess the predictive accuracy: Residual Mean
Square Error (RMSE) and Continuous Ranked Probability Score (CRPS).

## Incremental Value of ANC Data

We compare the predictive power of ANC data with that of survey data
alone to highlight the added value of ANC data.

## Usage

### Prerequisites

-   R (version 3.4.0 or higher)
-   R packages: `pacman`

### Running the Analysis

1.  Clone the repository:

2.  Install the required R packages: The setup.R script includes pacman
    to automatically install and load the required packages.

3.  If you have never used the 'greta' R package before you will need to run the function greta::install_greta_deps() before running the RMarkdown script to install the required dependencies.

4.  Run the RMarkdown script:
    rmarkdown::render("docs/malaria_analysis.Rmd")

Note: The RMarkdown script may take a long time to run. We have included
a pre-existing run in the repository for convenience. You can find the
generated HTML file in the docs folder.
The file is also available here: 

4.  View the generated HTML file in the docs folder.

## Repository Structure
W_Kenya_ANC_fitting/
├── docs/
│   ├── malaria_analysis.Rmd
│   ├── malaria_analysis.html
├── R/
│   ├── setup.R
│   ├── data_generation.R
│   ├── models.R
│   ├── summary_functions.R
├── .gitignore
├── README.md
├── LICENSE

-   **`docs/`**: Contains the RMarkdown script and generated output.

-   **`R/`**: Contains R scripts for setting up, generating data,
    fitting models, and summarizing model fits.

-   **`.gitignore`**: Specifies files and directories to be ignored by
    Git.

-   **`README.md`**: This file, providing an overview and usage
    instructions for the project.

## **Contributing**

Contributions are welcome! Please open an issue or submit a pull request
for any improvements or additions.

## **License**

This project is licensed under the MIT License - see the
License file for details.
