library(survival)
library(dplyr)
library(data.table)

# ---------------------------------------------------------
# Setup and Helper Functions
# ---------------------------------------------------------
# diet_scores <- c("MIND", "MedDiet")
# models <- list(m1 = "Age + Sex", m2 = "Age + Sex + BMI + Smoking")

setDT(NMR)

# Function to extract HR, 95% CI, and P-value
extract_cox_result <- function(coef_row) {
  data.frame(
    hr  = exp(coef_row["coef"]),
    lci = exp(coef_row["coef"] - 1.96 * coef_row["se(coef)"]),
    uci = exp(coef_row["coef"] + 1.96 * coef_row["se(coef)"]),
    p   = coef_row["Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------
# Main Analysis Pipeline
# ---------------------------------------------------------
all_results <- list()

for (score in diet_scores) {
  message("Processing: ", score)

  tryCatch({
    # 1. Variable Preparation
    if (!score %in% names(NMR)) stop("Score not found: ", score)

    std_col  <- paste0(score, "_std")
    rank_col <- paste0(score, "_rank")
    q_col    <- paste0(score, "_cat")

    # Standardization and Quintile distribution
    NMR[, (std_col) := as.numeric(scale(get(score)))]
    NMR[, (rank_col) := dplyr::ntile(get(std_col), 5)]
    NMR[, (q_col) := factor(get(rank_col), levels = 1:5, labels = paste0("Q", 1:5))]

    # 2. Model Iteration
    for (mod_name in names(models)) {
      covs <- models[[mod_name]]

      # Subset logic to calculate N and Events
      vars_needed <- unique(c("time_years", "Epilepsy", std_col, rank_col, q_col,
                            all.vars(as.formula(paste("~", covs)))))
      vars_needed <- vars_needed[vars_needed %in% names(NMR)]
      dat_sub <- NMR[, ..vars_needed]

      n_total <- nrow(dat_sub)
      n_event <- sum(dat_sub$Epilepsy == 1, na.rm = TRUE)

      # A. Continuous Analysis (Per 1-SD)
      f_cont <- as.formula(paste("Surv(time_years, Epilepsy) ~", std_col, "+", covs))
      fit_cont  <- coxph(f_cont, data = NMR)
      coef_cont <- summary(fit_cont)$coefficients

      if (std_col %in% rownames(coef_cont)) {
        res <- extract_cox_result(coef_cont[std_col, ])
        all_results[[paste(score, mod_name, "cont", sep = "_")]] <- data.frame(
          diet = score, model = mod_name, type = "Continuous", level = "Per 1-SD",
          n = n_total, n_event = n_event, hr = res$hr, lci = res$lci, uci = res$uci, p = res$p
        )
      }

      # B. Quintile Analysis (Q1 as Reference)
      f_cat <- as.formula(paste("Surv(time_years, Epilepsy) ~", q_col, "+", covs))
      fit_cat  <- coxph(f_cat, data = NMR)
      coef_cat <- summary(fit_cat)$coefficients

      # Q1 Reference Row
      all_results[[paste(score, mod_name, "Q1", sep = "_")]] <- data.frame(
        diet = score, model = mod_name, type = "Quintiles", level = "Q1 (Ref)",
        n = n_total, n_event = n_event, hr = 1, lci = NA, uci = NA, p = NA
      )

      # Q2-Q5 Estimates
      for (i in 2:5) {
        term <- paste0(q_col, "Q", i)
        if (term %in% rownames(coef_cat)) {
          res <- extract_cox_result(coef_cat[term, ])
          all_results[[paste(score, mod_name, paste0("Q", i), sep = "_")]] <- data.frame(
            diet = score, model = mod_name, type = "Quintiles", level = paste0("Q", i),
            n = n_total, n_event = n_event, hr = res$hr, lci = res$lci, uci = res$uci, p = res$p
          )
        }
      }

      # C. Trend Test (Linearity across quintiles)
      f_trend <- as.formula(paste("Surv(time_years, Epilepsy) ~", rank_col, "+", covs))
      fit_trend  <- coxph(f_trend, data = NMR)
      coef_trend <- summary(fit_trend)$coefficients

      if (rank_col %in% rownames(coef_trend)) {
        res <- extract_cox_result(coef_trend[rank_col, ])
        all_results[[paste(score, mod_name, "trend", sep = "_")]] <- data.frame(
          diet = score, model = mod_name, type = "Trend", level = "Per-quintile",
          n = n_total, n_event = n_event, hr = res$hr, lci = res$lci, uci = res$uci, p = res$p
        )
      }
    }
  }, error = function(e) message("Error in ", score, ": ", e$message))
}

# ---------------------------------------------------------
# Formatting and Export
# ---------------------------------------------------------
final_res <- bind_rows(all_results) %>%
  mutate(
    type = factor(type, levels = c("Continuous", "Quintiles", "Trend")),
    level = factor(level, levels = c("Per 1-SD", "Q1 (Ref)", "Q2", "Q3", "Q4", "Q5", "Per-quintile")),
    hr = round(hr, 3), lci = round(lci, 3), uci = round(uci, 3),
    p = ifelse(is.na(p), NA, format.pval(p, eps = 0.001, digits = 3))
  ) %>%
  arrange(diet, model, type, level)

write.csv(final_res, "cox_results_summary.csv", row.names = FALSE)
