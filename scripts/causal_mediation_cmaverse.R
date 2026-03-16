# ==============================================================================
# Script: 4-Way Decomposition Causal Mediation Analysis (CMAverse)
# Description: Automated scanning of multiple mediators with Cox proportional 
#              hazards models, 500 bootstraps, and checkpointing.
# ==============================================================================

library(CMAverse)
library(data.table)

# ==============================================================================
# 1. Configuration & Data Preparation
# ==============================================================================
# Assuming 'dt' is already loaded in your environment
setDT(dt)

# Define Variables
exposure_name <- "MIND_total_std"
outcome_name  <- "follow_up_years"
event_name    <- "Epilepsy"

cont_covs_raw <- c("Age.at.recruitment", "Townsend.deprivation.index.at.recruitment", 
                   "Sleep.duration", "MET", "BMI", "PRS")
cat_covs_raw  <- c("Sex", "Qualifications", "Ethnicity", "Smoking.status")

# Define Mediators to Scan (Replace indices with your actual mediator column names if needed)
old_names <- names(dt)
mediators_to_scan_raw <- old_names[18:57] # Captures total score + metabolites + inflammation

# Clean names for CMAverse compatibility (handles spaces/special characters)
new_names <- make.names(old_names)
setnames(dt, new_names)

all_covs_clean <- make.names(c(cont_covs_raw, cat_covs_raw))
cat_covs_clean <- make.names(cat_covs_raw)
mediators_to_scan_clean <- make.names(mediators_to_scan_raw)

# Ensure categorical covariates are factors
dt[, (cat_covs_clean) := lapply(.SD, as.factor), .SDcols = cat_covs_clean]

# ==============================================================================
# 2. Checkpoint & Resume Logic
# ==============================================================================
progress_file <- "MIND_Mediation_4way_Progress.csv"

if (file.exists(progress_file)) {
  done_data <- fread(progress_file)
  done_vars_raw <- done_data$Variable_Name
  
  # Filter out already processed mediators
  pending_idx <- !(mediators_to_scan_raw %in% done_vars_raw)
  vars_to_run_clean <- mediators_to_scan_clean[pending_idx]
  vars_to_run_raw   <- mediators_to_scan_raw[pending_idx]
  
  scan_results <- split(done_data, by = "Variable_Name")
  message(sprintf("Found checkpoint. Processed: %d | Remaining: %d", 
                  length(done_vars_raw), length(vars_to_run_clean)))
} else {
  vars_to_run_clean <- mediators_to_scan_clean
  vars_to_run_raw   <- mediators_to_scan_raw
  scan_results <- list()
  message(sprintf("Starting new analysis. Total mediators to process: %d", 
                  length(vars_to_run_clean)))
}

# ==============================================================================
# 3. Execution Loop
# ==============================================================================
message("--- Initialization complete. Starting robust single-core scan ---")

for (i in seq_along(vars_to_run_clean)) {
  current_var_clean <- vars_to_run_clean[i]
  current_var_raw   <- vars_to_run_raw[i]
  
  message(sprintf("[%d/%d] Processing: %s | Time: %s", 
                  i, length(vars_to_run_clean), current_var_raw, Sys.time()))
  
  tryCatch({
    # Run CMAverse 4-way decomposition
    fit <- cmest(
      data = dt,
      exposure = exposure_name,
      mediator = current_var_clean,
      outcome = outcome_name,
      event = event_name,
      basec = all_covs_clean,
      model = "rb",
      mreg = list("linear"),
      yreg = "coxph",
      EMint = TRUE,       # Enable Exposure-Mediator interaction
      mval = list(0),
      inference = "bootstrap",
      nboot = 500         # Bootstrap iterations
    )
    
    s <- summary(fit)$summary
    
    # Extract key metrics
    res_row <- data.frame(
      Variable_Name        = current_var_raw,
      Total_Effect_HR      = s["Rte", "Estimate"],
      Total_P              = s["Rte", "P.val"],
      Direct_Effect_HR     = s["Rcde", "Estimate"],
      Direct_P             = s["Rcde", "P.val"],
      Int_Med_Effect       = s["ERintmed", "Estimate"],
      Int_Med_P            = s["ERintmed", "P.val"],
      Pure_Indirect_Effect = s["ERpnie", "Estimate"],
      Prop_Mediated        = s["pm", "Estimate"],
      Prop_Interaction     = s["int", "Estimate"],
      stringsAsFactors     = FALSE
    )
    
    # Append to list and save checkpoint
    scan_results[[current_var_raw]] <- res_row
    current_progress <- do.call(rbind, scan_results)
    write.csv(current_progress, progress_file, row.names = FALSE)
    
  }, error = function(e) {
    message(sprintf("!!! Error processing %s: %s", current_var_raw, e$message))
  })
}

# ==============================================================================
# 4. Final Export
# ==============================================================================
if (length(scan_results) > 0) {
  final_scan_table <- do.call(rbind, scan_results)
  # Order by Total Effect P-value
  final_scan_table <- final_scan_table[order(final_scan_table$Total_P), ]
  
  final_out_file <- "MIND_Mediation_4way_Final_500Boot.csv"
  write.csv(final_scan_table, final_out_file, row.names = FALSE)
  message("--- 🎉 Analysis complete! Final results saved to: ", final_out_file, " ---")
} else {
  message("No results were generated. Check for errors in the run loop.")
}
