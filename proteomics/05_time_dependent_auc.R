#!/usr/bin/env Rscript

# =========================================================
# Time-dependent ROC analysis for incident epilepsy prediction
# Conventional covariates vs. protein-enriched model
#
# Author: Yuxi Zhang
# Description:
#   This script fits Cox models using conventional covariates alone
#   and conventional covariates plus 20 selected proteins, then compares
#   time-dependent AUCs from 1 to 12 years with bootstrap 95% CIs.
#
# Input:
#   A CSV file containing follow-up time, epilepsy status, covariates,
#   and 20 protein variables.
#
# Output:
#   1) PDF and PNG files of time-dependent ROC curves
#   2) CSV table containing AUCs and bootstrap confidence intervals
#
# Example:
#   Rscript scripts/timeROC_20proteins.R \
#     --input data/step3_20proteins_data.csv \
#     --output results/timeROC_20proteins \
#     --protein-start 16 \
#     --protein-end 35 \
#     --bootstrap 1000
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(timeROC)
  library(ggplot2)
  library(dplyr)
  library(pbapply)
})

# Cairo is optional. If unavailable, ggsave() will still export PNG,
# and PDF will use the default device.
has_cairo <- requireNamespace("Cairo", quietly = TRUE)

# =========================================================
# 0. Helper functions
# =========================================================

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  pos <- match(flag, args)
  if (!is.na(pos) && pos < length(args)) {
    return(args[pos + 1])
  }
  default
}

message_line <- function(...) {
  message(paste0(...))
}

to_numeric_if_exists <- function(dt, vars) {
  for (v in vars) {
    if (v %in% names(dt)) {
      dt[, (v) := as.numeric(get(v))]
    }
  }
  invisible(dt)
}

to_factor_if_exists <- function(dt, vars) {
  for (v in vars) {
    if (v %in% names(dt)) {
      dt[, (v) := as.factor(get(v))]
    }
  }
  invisible(dt)
}

check_required_columns <- function(dt, required_vars) {
  missing_vars <- setdiff(required_vars, names(dt))
  if (length(missing_vars) > 0) {
    stop(
      "The following required columns are missing from the input file:\n",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
}

get_ci_df <- function(prefix, times_vec, boot_mat) {
  bind_rows(lapply(times_vec, function(t) {
    x <- boot_mat[paste0(prefix, t), ]
    data.frame(
      Time_Val = t,
      low = unname(quantile(x, 0.025, na.rm = TRUE)),
      high = unname(quantile(x, 0.975, na.rm = TRUE))
    )
  }))
}

get_roc_df <- function(roc_obj, times_vec, model_name) {
  bind_rows(lapply(times_vec, function(t) {
    idx <- match(t, roc_obj$times)
    data.frame(
      FPR = roc_obj$FP[, idx],
      TPR = roc_obj$TP[, idx],
      Time = paste0(t, "-Year"),
      Model = model_name
    )
  }))
}

boot_auc_once <- function(data, times_vec, formula_baseline, formula_protein) {
  idx <- sample.int(nrow(data), size = nrow(data), replace = TRUE)
  d <- data[idx]

  empty_out <- c(
    setNames(rep(NA_real_, length(times_vec)), paste0("cov_", times_vec)),
    setNames(rep(NA_real_, length(times_vec)), paste0("prot_", times_vec)),
    setNames(rep(NA_real_, length(times_vec)), paste0("delta_", times_vec))
  )

  fit_b <- tryCatch(
    coxph(formula_baseline, data = d, ties = "breslow"),
    error = function(e) NULL
  )

  fit_p <- tryCatch(
    coxph(formula_protein, data = d, ties = "breslow"),
    error = function(e) NULL
  )

  if (is.null(fit_b) || is.null(fit_p)) {
    return(empty_out)
  }

  marker_b <- tryCatch(predict(fit_b, type = "lp"), error = function(e) NULL)
  marker_p <- tryCatch(predict(fit_p, type = "lp"), error = function(e) NULL)

  if (is.null(marker_b) || is.null(marker_p)) {
    return(empty_out)
  }

  roc_b <- tryCatch(
    timeROC(
      T = d$times,
      delta = d$Epilepsy,
      marker = marker_b,
      cause = 1,
      weighting = "marginal",
      times = times_vec,
      iid = FALSE
    ),
    error = function(e) NULL
  )

  roc_p <- tryCatch(
    timeROC(
      T = d$times,
      delta = d$Epilepsy,
      marker = marker_p,
      cause = 1,
      weighting = "marginal",
      times = times_vec,
      iid = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(roc_b) || is.null(roc_p)) {
    return(empty_out)
  }

  auc_b <- roc_b$AUC
  auc_p <- roc_p$AUC
  auc_d <- auc_p - auc_b

  c(
    setNames(auc_b, paste0("cov_", times_vec)),
    setNames(auc_p, paste0("prot_", times_vec)),
    setNames(auc_d, paste0("delta_", times_vec))
  )
}

# =========================================================
# 1. Command-line arguments
# =========================================================

input_file <- get_arg("--input", "data/step3_20proteins_data.csv")
output_dir <- get_arg("--output", "results/timeROC_20proteins")
protein_start <- as.integer(get_arg("--protein-start", "16"))
protein_end <- as.integer(get_arg("--protein-end", "35"))
B <- as.integer(get_arg("--bootstrap", "1000"))
seed <- as.integer(get_arg("--seed", "20260419"))

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file, "\n",
    "Please provide the file path with --input, for example:\n",
    "Rscript scripts/timeROC_20proteins.R --input data/step3_20proteins_data.csv",
    call. = FALSE
  )
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# =========================================================
# 2. Load data
# =========================================================

protein <- fread(input_file)
setDT(protein)

message_line("Input file: ", input_file)
message_line("Output directory: ", output_dir)
message_line("Rows: ", nrow(protein), "; columns: ", ncol(protein))

# =========================================================
# 3. Variable processing
# =========================================================

factor_vars <- c(
  "Sex",
  "Qualifications",
  "Ethnicity",
  "Smoking.status"
)

numeric_vars <- c(
  "Age.at.recruitment",
  "Townsend.deprivation.index.at.recruitment",
  "Sleep.duration",
  "MET",
  "BMI",
  "PRS",
  "times"
)

covariates <- c(
  "Sex",
  "Age.at.recruitment",
  "Townsend.deprivation.index.at.recruitment",
  "Qualifications",
  "Ethnicity",
  "Sleep.duration",
  "Smoking.status",
  "MET",
  "BMI",
  "PRS"
)

if (protein_start < 1 || protein_end > ncol(protein) || protein_start > protein_end) {
  stop(
    "Invalid protein column range. Please check --protein-start and --protein-end.",
    call. = FALSE
  )
}

protein_vars <- names(protein)[protein_start:protein_end]

required_vars <- c("times", "Epilepsy", covariates, protein_vars)
check_required_columns(protein, required_vars)

to_factor_if_exists(protein, factor_vars)
to_numeric_if_exists(protein, c(numeric_vars, protein_vars))

protein[, Epilepsy := as.integer(Epilepsy)]

# Ensure that Epilepsy is coded as 0/1.
if (!all(na.omit(unique(protein$Epilepsy)) %in% c(0L, 1L))) {
  stop("The outcome variable 'Epilepsy' must be coded as 0/1.", call. = FALSE)
}

# =========================================================
# 4. Complete-case analysis dataset
# =========================================================

vars_needed <- c("times", "Epilepsy", covariates, protein_vars)
dt <- na.omit(protein[, ..vars_needed])

message_line("Complete-case sample size: ", nrow(dt))
message_line("Incident epilepsy cases: ", sum(dt$Epilepsy == 1))
message_line("Non-cases: ", sum(dt$Epilepsy == 0))

if (sum(dt$Epilepsy == 1) == 0) {
  stop("No incident epilepsy cases were found after complete-case filtering.", call. = FALSE)
}

# =========================================================
# 5. Fit Cox models
# =========================================================

formula_baseline <- as.formula(
  paste0("Surv(times, Epilepsy) ~ ", paste(covariates, collapse = " + "))
)

formula_protein <- as.formula(
  paste0("Surv(times, Epilepsy) ~ ", paste(c(covariates, protein_vars), collapse = " + "))
)

message_line("Fitting conventional model...")
fit_baseline <- coxph(formula_baseline, data = dt, ties = "breslow")

message_line("Fitting protein-enriched model...")
fit_protein <- coxph(formula_protein, data = dt, ties = "breslow")

marker_baseline <- predict(fit_baseline, type = "lp")
marker_protein <- predict(fit_protein, type = "lp")

# =========================================================
# 6. Time-dependent ROC analysis
# =========================================================

eval_times <- 1:12

message_line("Calculating time-dependent ROC from 1 to 12 years...")

roc_baseline <- timeROC(
  T = dt$times,
  delta = dt$Epilepsy,
  marker = marker_baseline,
  cause = 1,
  weighting = "marginal",
  times = eval_times,
  iid = FALSE
)

roc_protein <- timeROC(
  T = dt$times,
  delta = dt$Epilepsy,
  marker = marker_protein,
  cause = 1,
  weighting = "marginal",
  times = eval_times,
  iid = FALSE
)

message_line("timeROC calculation completed.")

# =========================================================
# 7. Bootstrap 95% CIs
# =========================================================

set.seed(seed)
message_line("Running bootstrap for 95% CIs. B = ", B)
pboptions(type = "timer")

boot_list <- pblapply(seq_len(B), function(i) {
  boot_auc_once(
    data = dt,
    times_vec = eval_times,
    formula_baseline = formula_baseline,
    formula_protein = formula_protein
  )
})

boot_mat <- do.call(cbind, boot_list)
message_line("Bootstrap completed.")

ci_cov <- get_ci_df("cov_", eval_times, boot_mat)
ci_prot <- get_ci_df("prot_", eval_times, boot_mat)
ci_delta <- get_ci_df("delta_", eval_times, boot_mat)

# =========================================================
# 8. Prepare plotting data and AUC table
# =========================================================

plot_data <- bind_rows(
  get_roc_df(roc_baseline, eval_times, "Conventional Factors Only"),
  get_roc_df(roc_protein, eval_times, "Protein-enriched Model")
)

plot_data$Time <- factor(plot_data$Time, levels = paste0(eval_times, "-Year"))
plot_data$Model <- factor(
  plot_data$Model,
  levels = c("Conventional Factors Only", "Protein-enriched Model")
)

auc_labels <- data.frame(
  Time_Val = eval_times,
  AUC_Conv = roc_baseline$AUC,
  AUC_Integ = roc_protein$AUC
) %>%
  left_join(ci_cov %>% rename(CI_low_cov = low, CI_high_cov = high), by = "Time_Val") %>%
  left_join(ci_prot %>% rename(CI_low_int = low, CI_high_int = high), by = "Time_Val") %>%
  left_join(ci_delta %>% rename(CI_low_delta = low, CI_high_delta = high), by = "Time_Val") %>%
  mutate(
    Delta_AUC = AUC_Integ - AUC_Conv,
    Time = factor(paste0(Time_Val, "-Year"), levels = paste0(eval_times, "-Year")),
    Full_Label = paste0(
      "Conv AUC: ", sprintf("%.3f", AUC_Conv),
      " (", sprintf("%.3f", CI_low_cov), ", ", sprintf("%.3f", CI_high_cov), ")",
      "\nProt AUC: ", sprintf("%.3f", AUC_Integ),
      " (", sprintf("%.3f", CI_low_int), ", ", sprintf("%.3f", CI_high_int), ")",
      "\n\u0394AUC: ",
      ifelse(Delta_AUC >= 0, "+", ""),
      sprintf("%.3f", Delta_AUC),
      " (", sprintf("%.3f", CI_low_delta), ", ", sprintf("%.3f", CI_high_delta), ")"
    )
  )

auc_table <- auc_labels %>%
  transmute(
    Year = Time_Val,
    AUC_Conv,
    CI_low_cov,
    CI_high_cov,
    AUC_Prot = AUC_Integ,
    CI_low_prot = CI_low_int,
    CI_high_prot = CI_high_int,
    Delta_AUC,
    CI_low_delta,
    CI_high_delta
  )

# =========================================================
# 9. Plot
# =========================================================

final_plot <- ggplot(plot_data, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_text(
    data = auc_labels,
    aes(x = 0.98, y = 0.02, label = Full_Label),
    inherit.aes = FALSE,
    size = 2.8,
    hjust = 1,
    vjust = 0,
    lineheight = 1.03,
    color = "black"
  ) +
  facet_wrap(~Time, ncol = 4) +
  scale_color_manual(
    name = "Predictive Models",
    values = c(
      "Conventional Factors Only" = "#E64B35FF",
      "Protein-enriched Model" = "#4DBBD5FF"
    )
  ) +
  labs(
    title = "Time-dependent ROC Analysis: 1 to 12 Year Prediction",
    x = "1 - Specificity (FPR)",
    y = "Sensitivity (TPR)"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

print(final_plot)

# =========================================================
# 10. Export outputs
# =========================================================

pdf_file <- file.path(output_dir, "Figure_timeROC_1to12years_20proteins_4perrow_final.pdf")
png_file <- file.path(output_dir, "Figure_timeROC_1to12years_20proteins_4perrow_final.png")
csv_file <- file.path(output_dir, "TimeROC_AUC_1to12years_20proteins_with_CI.csv")

if (has_cairo) {
  ggsave(
    pdf_file,
    plot = final_plot,
    width = 16,
    height = 12,
    device = Cairo::cairo_pdf
  )
} else {
  ggsave(
    pdf_file,
    plot = final_plot,
    width = 16,
    height = 12
  )
}

ggsave(
  png_file,
  plot = final_plot,
  width = 16,
  height = 12,
  dpi = 600
)

fwrite(as.data.table(auc_table), csv_file)

message_line("Done. Outputs saved to: ", output_dir)
