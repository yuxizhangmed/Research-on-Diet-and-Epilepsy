# ==============================================================================
# Script: 03_triple_cox_filtering.R
# Description: Validates independent proteins across three progressive Cox models.
# ==============================================================================

library(data.table)
library(survival)
library(broom)
library(pbapply)
library(dplyr)

# Assume cov_m1, cov_m2, cov_m3 are predefined based on your study design
# protein_list <- fread("output/Independent_Proteins_List.csv")$protein

run_triple_cox <- function(prot) {
  f1 <- as.formula(paste("Surv(times, Epilepsy) ~", prot, "+", paste(cov_m1, collapse = " + ")))
  f2 <- as.formula(paste("Surv(times, Epilepsy) ~", prot, "+", paste(cov_m2, collapse = " + ")))
  f3 <- as.formula(paste("Surv(times, Epilepsy) ~", prot, "+", paste(cov_m3, collapse = " + ")))
  
  res1 <- tidy(coxph(f1, data = dt), exponentiate = TRUE) %>% filter(term == prot) %>% mutate(model = "Model1")
  res2 <- tidy(coxph(f2, data = dt), exponentiate = TRUE) %>% filter(term == prot) %>% mutate(model = "Model2")
  res3 <- tidy(coxph(f3, data = dt), exponentiate = TRUE) %>% filter(term == prot) %>% mutate(model = "Model3")
  
  output <- rbind(res1, res2, res3)
  output$protein <- prot
  return(output)
}

pboptions(type = "txt", char = "=")
cox_triple_results <- pblapply(protein_list, run_triple_cox)

cox_results_final <- rbindlist(cox_triple_results)
final_intersection <- dcast(cox_results_final, protein ~ model, value.var = "p.value")

final_star_proteins <- final_intersection[Model1 < 0.05 & Model2 < 0.05 & Model3 < 0.05, protein]
write.csv(data.frame(protein = final_star_proteins), "output/Robust_Proteins_List.csv", row.names = FALSE)
