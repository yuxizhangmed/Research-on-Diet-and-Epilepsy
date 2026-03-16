# ==============================================================================
# Script: 02_protein_screening.R
# Description: High-throughput linear regression for differential expression 
#              and collinearity filtering (Pearson r <= 0.6).
# ==============================================================================

library(data.table)
library(pbapply)
library(broom)

dt <- fread("data/step2_proteomics_preprocessed.csv")
setDT(dt)

# Re-identify protein columns dynamically
meta_cols <- c("new_baseline_date", "Epilepsy_date", "Epilepsy", "MIND_total_std", "times")
covariates <- c("Age.at.recruitment", "Townsend.deprivation.index.at.recruitment", 
                "Sleep.duration", "BMI", "PRS", "MET", "Sex", "Qualifications", 
                "Ethnicity", "Smoking.status", "Alcohol.intake.frequency")
prot_cols <- setdiff(names(dt), c(meta_cols, covariates))

# 1. Differential Expression Analysis
cov_formula <- paste(c("MIND_total_std", covariates), collapse = " + ")

run_prot_de <- function(prot) {
  formula_res <- as.formula(paste(prot, "~", cov_formula))
  fit <- lm(formula_res, data = dt)
  res <- tidy(fit)
  target <- res[res$term == "MIND_total_std", ]
  target$protein <- prot
  return(target)
}

pboptions(type = "txt", char = ">")
de_results_list <- pblapply(prot_cols, run_prot_de)

de_results <- rbindlist(de_results_list)
de_results[, p.adj := p.adjust(p.value, method = "BH")]

fwrite(de_results, "output/Protein_DE_Final_Results.csv", row.names = FALSE)

# 2. Collinearity Filtering
setorder(de_results, p.value)
sig_proteins <- de_results[p.adj < 0.05, protein]

kept_proteins <- c()
pb <- txtProgressBar(min = 0, max = length(sig_proteins), style = 3)

for (i in seq_along(sig_proteins)) {
  current_p <- sig_proteins[i]
  if (length(kept_proteins) == 0) {
    kept_proteins <- c(kept_proteins, current_p)
  } else {
    cor_vals <- cor(dt[[current_p]], dt[, ..kept_proteins], use = "pairwise.complete.obs")
    if (max(abs(cor_vals), na.rm = TRUE) <= 0.6) {
      kept_proteins <- c(kept_proteins, current_p)
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

write.csv(data.frame(protein = kept_proteins), "output/Independent_Proteins_List.csv", row.names = FALSE)
