# ==============================================================================
# Script: Pre-diagnostic Metabolite Heatmap (Nested Case-Control)
# Description: Generates heatmaps for 2-year and 3-year pre-diagnostic windows
#              filtered by top features selected via Elastic Net.
# ==============================================================================

# Set working directory to your project root
# setwd("path/to/your/working_directory")

library(data.table)
library(ggplot2)
library(patchwork)

# ===========================================================
# 0. Load Data
# ===========================================================
# Define relative paths
data_path <- "data/step1_30metabolites_data.csv"
elastic_net_path <- "data/ElasticNet_alpha0.8_100seeds_lambda_1se_topN_metabolites.csv"

dt <- fread(data_path)
metabolite_df <- fread(elastic_net_path)

setDT(dt)
setDT(metabolite_df)

# ===========================================================
# 1. Extract Core Metabolites & Define Y-axis Order
# ===========================================================
# Assuming the metabolite column in Elastic Net results is named "Metabolite"
top_mets <- as.character(metabolite_df$Metabolite)

# Clean formatting (trim spaces, replace underscores) to match main dataset
top_mets <- trimws(gsub("_", " ", top_mets, fixed = TRUE))
top_mets <- top_mets[top_mets != ""]

# Reverse the order so the most important feature appears at the top of the y-axis
met_levels <- rev(top_mets) 

# ===========================================================
# 2. Extract 2-year and 3-year Result Tables
# ===========================================================
# ---- 2-year windows (columns 1:9) ----
dt2 <- as.data.table(dt[, 1:9])
setnames(dt2, c(
  "Metabolite", "time_frame", "n_case", "n_control",
  "case_mean", "control_mean", "diff_mean", "t_stat", "p_value"
))

# ---- 3-year windows (columns 11:19) ----
dt3 <- as.data.table(dt[, 11:19])
setnames(dt3, c(
  "Metabolite", "time_frame", "n_case", "n_control",
  "case_mean", "control_mean", "diff_mean", "t_stat", "p_value"
))

# ===========================================================
# 3. Clean and Filter Metabolites/Time Frames
# ===========================================================
dt2[, Metabolite := trimws(gsub("_", " ", as.character(Metabolite), fixed = TRUE))]
dt3[, Metabolite := trimws(gsub("_", " ", as.character(Metabolite), fixed = TRUE))]
dt2[, time_frame := trimws(as.character(time_frame))]
dt3[, time_frame := trimws(as.character(time_frame))]

# Filter by Top Elastic Net Metabolites ONLY
dt2 <- dt2[Metabolite %in% top_mets]
dt3 <- dt3[Metabolite %in% top_mets]

dt2 <- dt2[!is.na(time_frame) & time_frame != ""]
dt3 <- dt3[!is.na(time_frame) & time_frame != ""]

# ===========================================================
# 4. Define Complete Time Window Levels
# ===========================================================
tf2_levels <- c("10-12", "8-10", "6-8", "4-6", "2-4", "0-2")
tf3_levels <- c("9-12", "6-9", "3-6", "0-3")

dt2 <- dt2[time_frame %in% tf2_levels]
dt3 <- dt3[time_frame %in% tf3_levels]

# ===========================================================
# 5. Build Complete Grids (Handles Missing Values)
# ===========================================================
grid2 <- CJ(
  Metabolite = factor(met_levels, levels = met_levels),
  time_frame = tf2_levels,
  unique = TRUE
)

grid3 <- CJ(
  Metabolite = factor(met_levels, levels = met_levels),
  time_frame = tf3_levels,
  unique = TRUE
)

# Convert current Metabolite to factor
dt2[, Metabolite := factor(Metabolite, levels = met_levels)]
dt3[, Metabolite := factor(Metabolite, levels = met_levels)]

# Merge grids
dt2 <- merge(grid2, dt2, by = c("Metabolite", "time_frame"), all.x = TRUE, sort = FALSE)
dt3 <- merge(grid3, dt3, by = c("Metabolite", "time_frame"), all.x = TRUE, sort = FALSE)

# ===========================================================
# 6. Map Time Frames to X-axis Intervals
# ===========================================================
map2 <- data.table(
  time_frame = tf2_levels,
  xmin = c(-12, -10, -8, -6, -4, -2),
  xmax = c(-10,  -8, -6, -4, -2,  0)
)

map3 <- data.table(
  time_frame = tf3_levels,
  xmin = c(-12, -9, -6, -3),
  xmax = c( -9, -6, -3,  0)
)

dt2 <- merge(dt2, map2, by = "time_frame", all.x = TRUE, sort = FALSE)
dt3 <- merge(dt3, map3, by = "time_frame", all.x = TRUE, sort = FALSE)

# ===========================================================
# 7. Finalize Plotting Variables
# ===========================================================
dt2[, Metabolite := factor(as.character(Metabolite), levels = met_levels)]
dt3[, Metabolite := factor(as.character(Metabolite), levels = met_levels)]

# Assign significance stars (P < 0.05)
dt2[, sig := ifelse(!is.na(p_value) & p_value < 0.05, "*", "")]
dt3[, sig := ifelse(!is.na(p_value) & p_value < 0.05, "*", "")]

# Y-axis numeric positioning
dt2[, y_id := as.numeric(Metabolite)]
dt3[, y_id := as.numeric(Metabolite)]

# Shared color scale limits
fill_lim <- max(abs(c(dt2$diff_mean, dt3$diff_mean)), na.rm = TRUE)

# ===========================================================
# 8. Core Plotting Function
# ===========================================================
make_interval_heatmap <- function(d, xbreaks, xlabels, title_text, show_y = TRUE) {
  
  p <- ggplot(d) +
    geom_rect(aes(
      xmin = xmin, xmax = xmax,
      ymin = y_id - 0.5, ymax = y_id + 0.5,
      fill = diff_mean
    ),
    color = "white",
    linewidth = 0.5
    ) +
    
    geom_text(aes(
      x = (xmin + xmax) / 2,
      y = y_id,
      label = sig
    ),
    size = 4.0,
    fontface = "bold",
    color = "black"
    ) +
    
    scale_y_continuous(
      breaks = seq_along(levels(d$Metabolite)),
      labels = levels(d$Metabolite),
      expand = c(0, 0)
    ) +
    
    scale_x_continuous(
      breaks = xbreaks,
      labels = xlabels,
      expand = c(0, 0)
    ) +
    
    scale_fill_gradient2(
      low = "#3B6FB6",
      mid = "white",
      high = "#C73E3A",
      midpoint = 0,
      limits = c(-fill_lim, fill_lim),
      name = "Case - Control\n(Diff Mean)"
    ) +
    
    coord_cartesian(
      xlim = c(min(xbreaks), 0),
      clip = "off"
    ) +
    
    labs(
      x = "Years before diagnosis",
      y = NULL,
      title = title_text
    ) +
    
    theme_bw(base_size = 12) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.text.x  = element_text(size = 11, color = "black"),
      axis.text.y  = element_text(size = 9, color = "black"),
      panel.grid   = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      legend.title = element_text(face = "bold", size = 10),
      legend.text  = element_text(size = 9),
      plot.margin  = margin(10, 15, 10, 10)
    )
  
  if (!show_y) {
    p <- p + theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  
  return(p)
}

# ===========================================================
# 9. Generate Individual & Combined Plots
# ===========================================================
# ---- 2-year plot ----
p2_comb <- make_interval_heatmap(
  d = dt2,
  xbreaks = c(-12, -10, -8, -6, -4, -2, 0),
  xlabels = c("-12", "-10", "-8", "-6", "-4", "-2", "Onset"),
  title_text = "2-Year Windows",
  show_y = TRUE
) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

# ---- 3-year plot ----
p3_comb <- make_interval_heatmap(
  d = dt3,
  xbreaks = c(-12, -9, -6, -3, 0),
  xlabels = c("-12", "-9", "-6", "-3", "Onset"),
  title_text = "3-Year Windows",
  show_y = FALSE
) + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

# ---- Combined plot ----
p_combined <- p2_comb + p3_comb +
  plot_layout(widths = c(1.25, 1)) +
  plot_annotation(
    title = "Pre-Diagnostic Metabolite Differences Across Time Windows",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))
  )

# ===========================================================
# 10. Export Figures
# ===========================================================
# Optional: Create an output directory if it doesn't exist
# dir.create("output_figures", showWarnings = FALSE)

ggsave("heatmap_combined_top_metabolites.pdf", p_combined, width = 16, height = 9)
ggsave("heatmap_combined_top_metabolites.png", p_combined, width = 16, height = 9, dpi = 600)

print("Heatmap generation complete. Outputs saved successfully.")
