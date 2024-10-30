# The Prediction of Serum C-Reactive Protein Concentration Using Nonlinear Mixed-Effects Model

This repository contains the code and a dummy dataset for the paper **"The Prediction of Serum C-Reactive Protein Concentration Using Nonlinear Mixed-Effects Model"** by Suk Joo Bae, Gyu Ri Kim, Sun Geu Chae, and Yeesuk Kim. The study proposes a nonlinear mixed-effects (NME) model to describe and predict the nonlinear patterns of serum C-reactive protein (CRP) concentration in patients following hip surgery.

## Repository Contents

- `data/`: This folder includes a **dummy dataset** (`diagnosisAndCrptest2.csv`) with essential columns used in the study. The dummy data is generated to simulate CRP concentration and other patient-related variables.
- `script.R`: The primary R script for data preprocessing, model building, evaluation, and visualization as outlined in the paper.

## Overview of the Study

In clinical settings, CRP is a critical biomarker for inflammation and postoperative infection. This study employs a bi-exponential model with random effects to predict CRP levels in hip surgery patients, capturing individual patient variations effectively. Using Monte Carlo simulations, the model estimates the distribution of CRP levels over time, supporting clinical decisions related to infection risk and patient recovery.

## Requirements

The following R packages are required for running the analysis:

- `dplyr`
- `readr`
- `nlme`
- `MASS`
- `stringr`
- `genalg`
- `GA`
- `R.utils`
- `data.table`
- `doSNOW`
- `car`
- `multcompView`

To install these packages, run:

```R
packages <- c('dplyr', 'readr', 'nlme', 'MASS', 'stringr', 'genalg', 'GA', 'R.utils', 'data.table', 'doSNOW', 'car', 'multcompView')
install.packages(packages)
```

## Usage

To run the analysis:

1. Download the repository and ensure the dummy dataset is in the `data/` directory.
2. Open `script.R` in R or RStudio.
3. Run the code to preprocess the data, fit the model, and generate visualizations.

**Note**: The included dataset is a dummy dataset designed for demonstration purposes. For research applications, please use actual clinical data as described in the paper.

## Citation

If you use this code or data in your research, please cite this repository:

> Bae, S.J., Kim, G.R., Chae, S.G., Kim, Y. (2024). The Prediction of Serum C-Reactive Protein Concentration Using Nonlinear Mixed-Effects Model. GitHub repository, https://github.com/ekiben-c/The-Prediction-of-Serum-C-Reactive-Protein-Concentration-Using-Nonlinear-Mixed-Effects-Model
> 
## License

This repository is licensed under the MIT License.
