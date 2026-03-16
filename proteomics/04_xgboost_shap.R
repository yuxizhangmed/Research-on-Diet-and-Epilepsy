# ==============================================================================
# Script: 04_xgboost_shap.R
# Description: SHAP summary plots for top XGBoost predictors.
# ==============================================================================

library(data.table)
library(Matrix)
library(SHAPforxgboost)
library(ggplot2)
library(dplyr)

# 1. Feature Matrix Construction
all_features <- c(protein_cols, covariate_cols)
X_matrix <- model.matrix(~ . - 1, data = dt_clean[, ..all_features])

# (Assuming xgb_model is pre-trained)

# 2. SHAP Calculation
shap_long_data <- shap.prep(xgb_model = xgb_model, X_train = X_matrix, top_n = 20)

# 3. Clean Variable Names
shap_long_data <- shap_long_data %>%
  mutate(variable = as.character(variable)) %>%
  mutate(variable = if_else(variable == "Townsend.deprivation.index.at.recruitment", "TDI", variable)) %>%
  mutate(variable = factor(variable, levels = rev(unique(variable))))

# 4. Plotting
final_shap_plot <- shap.plot.summary(shap_long_data) +
  scale_color_gradient(low = "#FFCC33", high = "#660099", 
                       guide = guide_colorbar(title = "Feature Value", barwidth = 1.2, barheight = 10)) +
  labs(title = "SHAP Analysis: Top 20 Predictors of Epilepsy",
       x = "SHAP value (Impact on Epilepsy Risk)") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", size = 10))

ggsave("output/Final_SHAP_Fixed_TDI.png", final_shap_plot, width = 11, height = 8, dpi = 300)
