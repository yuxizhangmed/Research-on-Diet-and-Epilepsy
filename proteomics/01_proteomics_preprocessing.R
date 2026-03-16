# ==============================================================================
# Script: 01_proteomics_preprocessing.R
# Description: Dynamic median imputation, Z-score standardization for proteins 
#              and covariates, and survival time calculation.
# ==============================================================================

library(data.table)

# 1. Load Data
dt <- fread("data/step1_standarded_coviarates.csv")
setDT(dt)

# 2. Define Variables Dynamically (No hardcoded indices)
continuous_vars <- c("Age.at.recruitment", "Townsend.deprivation.index.at.recruitment", 
                     "Sleep.duration", "BMI", "PRS", "MET")
categorical_vars <- c("Sex", "Qualifications", "Ethnicity", 
                      "Smoking.status", "Alcohol.intake.frequency")

# Define all metadata and covariate columns
meta_cols <- c("new_baseline_date", "Epilepsy_date", "Epilepsy", "MIND_total_std",
               continuous_vars, categorical_vars)

# Automatically identify protein columns (everything else)
prot_cols <- setdiff(names(dt), meta_cols)

# 3. Protein Data Imputation and Standardization
for (col in prot_cols) {
  m_val <- median(dt[[col]], na.rm = TRUE)
  dt[is.na(get(col)), (col) := m_val]
  
  std_val <- sd(dt[[col]], na.rm = TRUE)
  if(std_val > 0) {
    dt[, (col) := (get(col) - mean(get(col))) / std_val]
  } else {
    dt[, (col) := 0]
  }
}

# 4. Covariate Standardization and Factorization
for (var in continuous_vars) {
  dt[, (var) := as.vector(scale(get(var)))]
}
dt[, (categorical_vars) := lapply(.SD, as.factor), .SDcols = categorical_vars]

# 5. Survival Time Calculation (Years)
dt[, new_baseline_date := as.Date(new_baseline_date)]
dt[, Epilepsy_date := as.Date(Epilepsy_date)]
study_end_date <- as.Date("2023-03-31")

dt[, times := as.numeric(study_end_date - new_baseline_date) / 365.25]
dt[Epilepsy == 1, times := as.numeric(Epilepsy_date - new_baseline_date) / 365.25]

fwrite(dt, "data/step2_proteomics_preprocessed.csv")
