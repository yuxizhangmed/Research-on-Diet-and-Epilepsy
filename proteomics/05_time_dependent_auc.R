# ==============================================================================
# Script: 05_time_dependent_auc.R
# Description: ROC comparison between clinical models and protein-enhanced models.
# ==============================================================================

library(timeROC)
library(ggplot2)

eval_times <- c(1, 3, 5, 8, 10, 12)

# 1. Compute ROC
res_enhanced <- timeROC(T = dt_clean$times, delta = dt_clean$Epilepsy,
                        marker = dt_clean$score_enhanced, cause = 1,
                        times = eval_times, iid = TRUE)

res_conv <- timeROC(T = dt_clean$times, delta = dt_clean$Epilepsy,
                    marker = dt_clean$score_conv, cause = 1,
                    times = eval_times, iid = TRUE)

summary_table <- data.frame(
  Year = eval_times,
  AUC_Enhanced = res_enhanced$AUC * 100,
  AUC_Conv = res_conv$AUC * 100
)
summary_table$Delta <- summary_table$AUC_Enhanced - summary_table$AUC_Conv

# 2. 10-Year ROC Plot
target_idx <- 5 
df_roc <- data.frame(
  FP_E = res_enhanced$FP[, target_idx], TP_E = res_enhanced$TP[, target_idx],
  FP_C = res_conv$FP[, target_idx], TP_C = res_conv$TP[, target_idx]
)

final_roc_plot <- ggplot(df_roc) +
  geom_line(aes(x = FP_E, y = TP_E, color = "Protein + Clinical"), linewidth = 1.2) +
  geom_line(aes(x = FP_C, y = TP_C, color = "Clinical Only (Model 3)"), linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Protein + Clinical" = "#5D87B1", "Clinical Only (Model 3)" = "#E68A5D")) +
  labs(title = "10-Year Epilepsy Risk Prediction",
       subtitle = paste0("AUC Improvement: ", round(summary_table$Delta[target_idx], 2), "%"),
       x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", color = "Models") +
  theme_bw() +
  theme(legend.position = c(0.75, 0.2), legend.background = element_blank())

ggsave("output/Final_ROC_Comparison_10y.png", final_roc_plot, width = 6, height = 6, dpi = 300)
