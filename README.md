# The Prediction of Serum C-Reactive Protein Concentration Using Nonlinear Mixed-Effects Model

This repository contains the code and a dummy dataset for the paper **"The Prediction of Serum C-Reactive Protein Concentration Using Nonlinear Mixed-Effects Model"** by Suk Joo Bae, Gyu Ri Kim, Sun Geu Chae, and Yeesuk Kim. The study proposes a nonlinear mixed-effects (NME) model to describe and predict the nonlinear patterns of serum C-reactive protein (CRP) concentration in patients following hip surgery.

## Repository Contents

- `data/`: This folder includes a **dummy dataset** (`diagnosisAndCrptest2.csv`) with essential columns used in the study. The dummy data is generated to simulate CRP concentration and other patient-related variables.
- `script.R`: The primary R script for data preprocessing, model building, evaluation, and visualization as outlined in the paper.
- `README.md`: This file.

## Overview of the Study

In clinical settings, CRP is a critical biomarker for inflammation and postoperative infection. This study employs a bi-exponential model with random effects to predict CRP levels in hip surgery patients, capturing individual patient variations effectively. Using Monte Carlo simulations, the model estimates the distribution of CRP levels over time, supporting clinical decisions related to infection risk and patient recovery.

## Code and Analysis Steps

1. **Data Preprocessing**: The dataset is filtered and transformed to include relevant CRP concentration data from patients meeting inclusion criteria. The R script (`script.R`) loads, cleans, and prepares the data.

2. **Modeling**: A nonlinear mixed-effects (NME) model with bi-exponential functions is applied to capture CRP concentration changes over time. 
   
   - **Model Fitting**: The model is fitted using the `nlme` R package, leveraging Lindstrom and Bates' (LB) algorithm for parameter estimation.
   - **Simulation**: Monte Carlo simulation generates simulated CRP paths to evaluate prediction accuracy and analyze concentration patterns.
   - **Metrics Calculation**: Model performance is evaluated using metrics such as Mean Squared Error (MSE), Root Mean Squared Error (RMSE), and Mean Absolute Error (MAE).

3. **Statistical Analysis**: An ANOVA test compares CRP levels across age groups, followed by Tukey-Kramer post-hoc comparisons to analyze statistically significant differences.

4. **Visualization**: The R script includes code to plot CRP concentration trends, model predictions, and simulation results.

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

If you use this code or data in your research, please cite the paper:

> Bae, S.J., Kim, G.R., Chae, S.G., Kim, Y. (2024). The Prediction of Serum C-Reactive Protein Concentration Using Nonlinear Mixed-Effects Model. *IEEE Access*. DOI: 10.1109/ACCESS.2024.0429000

## License

This repository is licensed under the MIT License.
