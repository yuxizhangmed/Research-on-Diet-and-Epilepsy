#' @title Analyze Diet Scores and Epilepsy Risk
#' @description Perform Cox regression for continuous, quintile, and trend analysis.
#' @author Yuxi Zhang
#' @date 2026-03-16

library(survival)
library(dplyr)
library(data.table)

# --- 1. 核心分析函数 ---
run_cox_analysis <- function(data, score_name, event_col, time_col, cov_formula, model_label) {
  
  results <- list()
  
  # A. 连续变量分析 (Per 1 SD)
  std_col <- paste0(score_name, "_std")
  data[[std_col]] <- as.numeric(scale(data[[score_name]]))
  
  f_cont <- as.formula(paste("Surv(", time_col, ",", event_col, ") ~", std_col, "+", cov_formula))
  fit_cont <- coxph(f_cont, data = data)
  res_cont <- summary(fit_cont)$coefficients[std_col, ]
  
  results$cont <- data.frame(
    Diet = score_name, Model = model_label, Type = "Continuous", Level = "Per 1 SD",
    HR = exp(res_cont["coef"]),
    LCI = exp(res_cont["coef"] - 1.96 * res_cont["se(coef)"]),
    UCI = exp(res_cont["coef"] + 1.96 * res_cont["se(coef)"]),
    P_Value = res_cont["Pr(>|z|)"]
  )
  
  # B. 五分位分类与趋势分析
  rank_col <- paste0(score_name, "_rank")
  data[[rank_col]] <- dplyr::ntile(data[[score_name]], 5)
  data$q_factor <- factor(data[[rank_col]], levels = 1:5, labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))
  
  # 分类回归 (Categorical)
  f_cat <- as.formula(paste("Surv(", time_col, ",", event_col, ") ~ q_factor +", cov_formula))
  fit_cat <- coxph(f_cat, data = data)
  
  # 提取 Q2-Q5 (Q1 为 Reference)
  for (i in 2:5) {
    term <- paste0("q_factorQ", i)
    res_q <- summary(fit_cat)$coefficients[term, ]
    results[[paste0("Q", i)]] <- data.frame(
      Diet = score_name, Model = model_label, Type = "Quintiles", Level = paste0("Q", i),
      HR = exp(res_q["coef"]),
      LCI = exp(res_q["coef"] - 1.96 * res_q["se(coef)"]),
      UCI = exp(res_q["coef"] + 1.96 * res_q["se(coef)"]),
      P_Value = res_q["Pr(>|z|)"]
    )
  }
  
  # C. 趋势检验 (Trend)
  f_trend <- as.formula(paste("Surv(", time_col, ",", event_col, ") ~", rank_col, "+", cov_formula))
  fit_trend <- coxph(f_trend, data = data)
  res_trend <- summary(fit_trend)$coefficients[rank_col, ]
  
  results$trend <- data.frame(
    Diet = score_name, Model = model_label, Type = "P for Trend", Level = "Per Quintile",
    HR = exp(res_trend["coef"]),
    LCI = exp(res_trend["coef"] - 1.96 * res_trend["se(coef)"]),
    UCI = exp(res_trend["coef"] + 1.96 * res_trend["se(coef)"]),
    P_Value = res_trend["Pr(>|z|)"]
  )
  
  return(bind_rows(results))
}

# --- 2. 主循环执行部分 ---

# 设定变量（建议在单独的 config 脚本或脚本头部定义）
# diet_scores <- c("MIND_Score", "MedDiet_Score") 
# models <- list(Model1 = "age + sex", Model2 = "age + sex + smoking + bmi")

all_output <- list()

for (score in diet_scores) {
  message(">>> Processing: ", score)
  for (mod_name in names(models)) {
    tmp_res <- run_cox_analysis(NMR, score, "Epilepsy", "time_years", models[[mod_name]], mod_name)
    all_output[[paste(score, mod_name)]] <- tmp_res
  }
}

# --- 3. 结果美化与保存 ---
final_table <- bind_rows(all_output) %>%
  mutate(
    across(c(HR, LCI, UCI), ~ round(.x, 3)),
    P_Value = format.pval(P_Value, eps = 0.001, digits = 3)
  )

write.csv(final_table, "Cox_Results_Summary.csv", row.names = FALSE)
message("Done! Results saved.")
