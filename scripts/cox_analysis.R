library(survival)
library(dplyr)
library(data.table)

# --- Configuration ---
# diet_scores <- c("MIND", "MedDiet")
# models <- list(m1 = "age + sex", m2 = "age + sex + bmi + smoke")

all_results <- list()

for (score in diet_scores) {
  message("Analyzing score: ", score)
  
  # Normalize and define quintiles
  std_col <- paste0(score, "_std")
  NMR[, (std_col) := scale(get(score))]
  
  rank_col <- paste0(score, "_rank")
  NMR[, (rank_col) := ntile(get(std_col), 5)]
  
  q_col <- paste0(score, "_cat")
  NMR[, (q_col) := factor(get(rank_col), levels = 1:5, labels = paste0("Q", 1:5))]
  
  # Loop through adjustment models
  for (mod_name in names(models)) {
    covs <- models[[mod_name]]
    
    # 1. Continuous (Per SD)
    f_cont <- as.formula(paste("Surv(time_years, Epilepsy) ~", std_col, "+", covs))
    fit_cont <- coxph(f_cont, data = NMR)
    res_cont <- summary(fit_cont)$coefficients[std_col, ]
    
    all_results[[paste(score, mod_name, "cont")]] <- data.frame(
      diet = score, model = mod_name, type = "Continuous", level = "Per 1-SD",
      hr = exp(res_cont["coef"]),
      lci = exp(res_cont["coef"] - 1.96 * res_cont["se(coef)"]),
      uci = exp(res_cont["coef"] + 1.96 * res_cont["se(coef)"]),
      p = res_cont["Pr(>|z|)"]
    )
    
    # 2. Quintiles (Categorical)
    f_cat <- as.formula(paste("Surv(time_years, Epilepsy) ~", q_col, "+", covs))
    fit_cat <- coxph(f_cat, data = NMR)
    
    # Reference group (Q1)
    all_results[[paste(score, mod_name, "Q1")]] <- data.frame(
      diet = score, model = mod_name, type = "Quintiles", level = "Q1 (Ref)",
      hr = 1.0, lci = NA, uci = NA, p = NA
    )
    
    # Q2-Q5
    for (i in 2:5) {
      term <- paste0(q_col, "Q", i)
      if (term %in% rownames(summary(fit_cat)$coefficients)) {
        res_q <- summary(fit_cat)$coefficients[term, ]
        all_results[[paste(score, mod_name, "Q", i)]] <- data.frame(
          diet = score, model = mod_name, type = "Quintiles", level = paste0("Q", i),
          hr = exp(res_q["coef"]),
          lci = exp(res_q["coef"] - 1.96 * res_q["se(coef)"]),
          uci = exp(res_q["coef"] + 1.96 * res_q["se(coef)"]),
          p = res_q["Pr(>|z|)"]
        )
      }
    }
    
    # 3. P for Trend
    f_trend <- as.formula(paste("Surv(time_years, Epilepsy) ~", rank_col, "+", covs))
    fit_trend <- coxph(f_trend, data = NMR)
    res_trend <- summary(fit_trend)$coefficients[rank_col, ]
    
    all_results[[paste(score, mod_name, "trend")]] <- data.frame(
      diet = score, model = mod_name, type = "Trend", level = "Per-quintile",
      hr = exp(res_trend["coef"]),
      lci = exp(res_trend["coef"] - 1.96 * res_trend["se(coef)"]),
      uci = exp(res_trend["coef"] + 1.96 * res_trend["se(coef)"]),
      p = res_trend["Pr(>|z|)"]
    )
  }
}

# --- Export ---
final_res <- bind_rows(all_results) %>%
  mutate(across(c(hr, lci, uci), ~ round(.x, 3)),
         p = format.pval(p, eps = 0.001, digits = 3))

write.csv(final_res, "diet_epilepsy_cox_results.csv", row.names = FALSE)
