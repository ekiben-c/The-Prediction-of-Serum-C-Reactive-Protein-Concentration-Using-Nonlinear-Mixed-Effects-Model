# # The Prediction of Serum C-Reactive Protein Concentration Using Nonlinear Mixed-Effects Model

# This repository contains the code and a dummy dataset for the paper **"The Prediction of Serum C-Reactive Protein Concentration Using Nonlinear Mixed-Effects Model"** by Suk Joo Bae, Gyu Ri Kim, Sun Geu Chae, and Yeesuk Kim. The study proposes a nonlinear mixed-effects (NME) model to describe and predict the nonlinear patterns of serum C-reactive protein (CRP) concentration in patients following hip surgery.

# ## Repository Contents

# - `data/`: This folder includes a **dummy dataset** (`diagnosisAndCrptest2.csv`) with essential columns used in the study. The dummy data is generated to simulate CRP concentration and other patient-related variables.
# - `script.R`: The primary R script for data preprocessing, model building, evaluation, and visualization as outlined in the paper.
# - `README.md`: This file.

# ## Overview of the Study

# In clinical settings, CRP is a critical biomarker for inflammation and postoperative infection. This study employs a bi-exponential model with random effects to predict CRP levels in hip surgery patients, capturing individual patient variations effectively. Using Monte Carlo simulations, the model estimates the distribution of CRP levels over time, supporting clinical decisions related to infection risk and patient recovery.

# ## Requirements

# The following R packages are required for running the analysis:

# - `dplyr`
# - `readr`
# - `nlme`
# - `MASS`
# - `stringr`
# - `genalg`
# - `GA`
# - `R.utils`
# - `data.table`
# - `doSNOW`
# - `car`
# - `multcompView`

# To install these packages, run:

# ```R
# packages <- c('dplyr', 'readr', 'nlme', 'MASS', 'stringr', 'genalg', 'GA', 'R.utils', 'data.table', 'doSNOW', 'car', 'multcompView')
# install.packages(packages)
# ```

# ## Usage

# To run the analysis:

# 1. Download the repository and ensure the dummy dataset is in the `data/` directory.
# 2. Open `script.R` in R or RStudio.
# 3. Run the code to preprocess the data, fit the model, and generate visualizations.

# **Note**: The included dataset is a dummy dataset designed for demonstration purposes. For research applications, please use actual clinical data as described in the paper.

# ## Citation

# If you use this code or data in your research, please cite the paper:

# > Bae, S.J., Kim, G.R., Chae, S.G., Kim, Y. (2024). The Prediction of Serum C-Reactive Protein Concentration Using Nonlinear Mixed-Effects Model. *IEEE Access*. DOI: 10.1109/ACCESS.2024.0429000

# ## License

# This repository is licensed under the MIT License.


# ================================================================================
# 1. Setup
# ================================================================================

# Set the working directory to the script's location and clear the environment
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls()); cat('\014'); gc()

# Define the packages needed for the script
packages <- c('dplyr', 'readr', 'nlme', 'MASS', 'stringr', 'rstudioapi',
              'genalg', 'GA', 'R.utils', 'data.table', 'doSNOW', 'car', 'multcompView')

# Install and load packages if necessary
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  require(pkg, character.only = TRUE)
})

# Set up a progress bar for later tasks
pb <- txtProgressBar(max = 100, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# ================================================================================
# 2. Data Preparation
# ================================================================================

# Load the dataset and select relevant columns
df <- read.csv("diagnosisAndCrptest2.csv") %>%
  dplyr::select(uniqueid, result, surgeryToTestDays, ageAtSurgery)

# Filter data for completeness and relevant time range, and remove duplicates
df <- df[complete.cases(df),] %>%
  filter(between(surgeryToTestDays, 0, 30)) %>%
  unique()

# Retain subjects with at least 5 measurements
df <- df %>%
  group_by(uniqueid) %>%
  filter(n() >= 5) %>%
  ungroup()

# Rename columns for consistency
colnames(df) <- c('Subject', 'y', 'Time', 'age')

# Categorize age into age groups and order data by time
df <- df %>%
  mutate(ageGroup = cut(age, breaks = seq(0, 110, by = 10), right = FALSE,
                        labels = paste0(seq(0, 100, by = 10)))) %>%
  arrange(Time)

# Convert the data to a groupedData object
set.seed(1234)
gdf <- groupedData(y ~ Time | Subject, df, order.groups = FALSE)

# ================================================================================
# 3. Model Fitting
# ================================================================================

# Define initial parameter guesses for nonlinear model
init <- mvrnorm(1e4, c(-8.23, 0.3442, 12.845, -2.423), diag(0.01, 4))

# Iterative model fitting with error handling
nlm <- 'Reset'
while (nlm == 'Reset') {
  nlm <- tryCatch({
    best <- init[sample(1:1e4, 1), ]
    nlm <- nlme(y ~ SSbiexp(Time, A1, lrc1, A2, lrc2), data = gdf,
                fixed = A1 + lrc1 + A2 + lrc2 ~ 1,
                random = pdDiag(A1 + lrc1 + A2 + lrc2 ~ 1),
                start = best, na.action = na.omit,
                control = list(msVerbose = TRUE, opt = "nlminb", msMaxIter = 700,
                               pnlsMaxIter = 10, tolerance = 0.5))
  }, error = function(e) {
    print('Reset')
    'Reset'
  })
}
bp <- c(best[1], -exp(best[2]), best[3], -exp(best[4]))
bp

# ================================================================================
# 4. Parsimonious Model
# ================================================================================

# Simplify model by updating the random effect structure
nlm3 <- update(nlm, random = pdDiag(A1 + A2 + lrc2 ~ 1),
               control = list(tolerance = 0.1))
nlm2 <- update(nlm, random = pdDiag(A2 + lrc2 ~ 1),
               control = list(tolerance = 1e-100))

# Select the final model and calculate variance
fm <- nlm2
anova(nlm, fm)
summary(fm)$sigma^2
data.frame(VarCorr(fm)[])

# ================================================================================
# 5. Define Metrics for Model Evaluation
# ================================================================================

# Define a function to calculate evaluation metrics for predictions
Metrics <- function(gdf, pred) {
  MSE <- mean((gdf$y - pred)^2)
  RMSE <- sqrt(MSE)
  MAE <- mean(abs(gdf$y - pred))
  RAE <- mean(abs(gdf$y - pred) / abs(gdf$y - mean(gdf$y)))
  MRE <- mean(abs(gdf$y - pred) / gdf$y)
  data.frame(MSE = MSE, RMSE = RMSE, MAE = MAE, RAE = RAE, MRE = MRE)
}

# Calculate and display metrics for the final model
final_mu <- fm$coefficients$fixed
round(Metrics(gdf, fm$fitted[, 1]), 3)
round(Metrics(gdf, fm$fitted[, 2]), 3)

# ================================================================================
# 6. Visualization of Raw Data
# ================================================================================

# Plot raw data for selected subjects
id <- c(133, 134, 154, 343, 444, 200, 503, 453, 84, 347)
type <- c(2:8, 10:12)
col <- c(2:8, 'darkgreen', 'darkred', 'darkblue')

par(mfrow = c(2, 5))
sapply(1:length(id), function(i) {
  pdf <- df[df$Subject == id[i],]
  plot(pdf$Time, pdf$y, col = col[i], pch = type[i], cex = 1,
       xlim = c(0, 30), ylim = c(0, 20), lty = type[i], type = 'b',
       xlab = 'Time of measurements (Days)', ylab = 'Concentration (mg/ml)')
  grid()
})
par(mfrow = c(1, 1))

# ================================================================================
# 7. Model Simulation and Confidence Interval
# ================================================================================

# Define function to predict CRP concentration and simulate model
pred <- function(t, par, b) {
  (par[1] + b[1]) * exp(-exp(par[2] + b[2]) * t) + (par[3] + b[3]) * exp(-exp(par[4] + b[4]) * t)
}

# Simulate model and plot prediction intervals
simul <- function(obj) {
  set.seed(1); n <- 100
  mu <- obj$coefficients$fixed
  varfix <- obj$varFix
  sim <- mvrnorm(n, mu, varfix)
  histo <- sapply(1:n, function(j) {
    optim(par = 0, fn = function(x) abs(pred(x, sim[j, ], rep(0, 4)) - 0.3),
          lower = 0, upper = 100, method = 'Brent')$par
  })
  list(sim = sim, histo = histo)
}
result <- simul(fm)

# Histogram of simulated results
hist(result$histo, breaks = 50, col = 1, xlim = c(min(result$histo), 50),
     prob = TRUE, xlab = 'Time of measurements (Days)', main = '')

# ================================================================================
# 8. Tukey's Post Hoc Test and Age Group Visualization
# ================================================================================

# Run ANOVA and Tukey's test on age groups
ph <- transform(df, ageGroup = factor(ageGroup))
aovm <- aov(y ~ ageGroup, data = ph)
TukeyHSD <- TukeyHSD(aovm, unbalanced = TRUE)
Tx <- data.frame(round(TukeyHSD$ageGroup, 3))
comp <- unlist(lapply(str_split(rownames(Tx), '-'), function(x) paste0(x[2], '-', x[1])))
Tx <- cbind(comp, Tx)

# Age group boxplot with Tukey results
boxplot(ph$y ~ ph$ageGroup, outline = FALSE, variation = 'IQR',
        ylim = c(0, 25), col = 2:8, xlab = 'Age group', ylab = 'Concentration (mg/ml)')
legend('topright', inset = c(0.01, 0.01),
       paste0(sort(unique(df$ageGroup)), 's'), col = 2:8, pch = 15)

# ================================================================================
# End of Script
# ================================================================================
