# ============================================================
# Plotting script for BUGS vs BUGSRHS toy-example introduction figure
#
# Creates:
#   1) Intro panel 1: marginal scores
#   2) Intro panel 2: posterior local-scale boxplots (Guided vs Unguided)
#   3) Intro panel 3: posterior exceedance probabilities (Guided vs Unguided)
#   4) Combined 3x1 introduction panel
#
# Expected files in Results/:
#   Output_marginal_scores_n_<n>_p_<p>_rep_<rep>.csv
#   Output_BUGS_shrinkage_draws_n_<n>_p_<p>_rep_<rep>.csv
#   Output_BUGS_variable_summary_n_<n>_p_<p>_rep_<rep>.csv
#   Output_BUGSRHS_shrinkage_draws_n_<n>_p_<p>_rep_<rep>.csv
#   Output_BUGSRHS_variable_summary_n_<n>_p_<p>_rep_<rep>.csv
# ============================================================

rm(list = ls())
setwd("U:/BUGS/Marginal score benefits")

# Display range for local-scale plot (visual zoom only)
local_scale_min <- 1e-2
local_scale_max <- 1e3

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(patchwork)
  library(scales)
  library(cowplot)
  library(grid)
})

# -----------------------------
# User settings
# -----------------------------
n   <- 100
p   <- 10
rep <- 1

results_dir <- "Results"
out_dir <- file.path(results_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Threshold used in posterior exceedance probability
sel_thresh <- 0.01

# -----------------------------
# File paths
# -----------------------------
file_scores <- file.path(
  results_dir,
  sprintf("Output_marginal_scores_n_%d_p_%d_rep_%d.csv", n, p, rep)
)

file_bugs_draws <- file.path(
  results_dir,
  sprintf("Output_BUGS_shrinkage_draws_n_%d_p_%d_rep_%d.csv", n, p, rep)
)
file_bugs_var <- file.path(
  results_dir,
  sprintf("Output_BUGS_variable_summary_n_%d_p_%d_rep_%d.csv", n, p, rep)
)

file_rhs_draws <- file.path(
  results_dir,
  sprintf("Output_BUGSRHS_shrinkage_draws_n_%d_p_%d_rep_%d.csv", n, p, rep)
)
file_rhs_var <- file.path(
  results_dir,
  sprintf("Output_BUGSRHS_variable_summary_n_%d_p_%d_rep_%d.csv", n, p, rep)
)

needed <- c(file_scores, file_bugs_draws, file_bugs_var, file_rhs_draws, file_rhs_var)
missing_files <- needed[!file.exists(needed)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

# -----------------------------
# Read data
# -----------------------------
scores    <- read_csv(file_scores, show_col_types = FALSE)
bugs_draw <- read_csv(file_bugs_draws, show_col_types = FALSE)
bugs_var  <- read_csv(file_bugs_var, show_col_types = FALSE)
rhs_draw  <- read_csv(file_rhs_draws, show_col_types = FALSE)
rhs_var   <- read_csv(file_rhs_var, show_col_types = FALSE)

# -----------------------------
# Harmonize and enrich data
# -----------------------------
truth_df <- scores %>%
  transmute(
    index = as.integer(index),
    is_signal = as.logical(is_signal),
    truth_label = if_else(as.logical(is_signal), "True variable", "Noise variable"),
    marginal_score = marginal_score,
    z_score_clipped = z_score_clipped
  )

# Guided model: use saved z_score_clipped to construct lambda_j^2 exp(eta z_j*)
bugs_draw <- bugs_draw %>%
  mutate(index = as.integer(index)) %>%
  left_join(truth_df %>% select(index, truth_label), by = "index") %>%
  mutate(
    model = "Guided",
    local_scale = lambda2_draw * exp(eta_draw * z_score_clipped)
  )

# Unguided model: baseline uses lambda_j^2
rhs_draw <- rhs_draw %>%
  mutate(index = as.integer(index)) %>%
  left_join(truth_df %>% select(index, truth_label), by = "index") %>%
  mutate(
    model = "Unguided",
    local_scale = lambda2_draw
  )

bugs_var <- bugs_var %>%
  mutate(index = as.integer(index)) %>%
  left_join(truth_df, by = "index") %>%
  mutate(model = "Guided")

rhs_var <- rhs_var %>%
  mutate(index = as.integer(index)) %>%
  left_join(truth_df, by = "index") %>%
  mutate(model = "Unguided")

# Variable labels
var_levels <- sprintf("X%d", 1:p)
var_map <- tibble(
  index = 1:p,
  var = factor(sprintf("X%d", 1:p), levels = var_levels)
)

scores <- scores %>%
  mutate(index = as.integer(index)) %>%
  left_join(var_map, by = "index") %>%
  mutate(truth_label = if_else(as.logical(is_signal), "True variable", "Noise variable"))

bugs_draw <- bugs_draw %>% left_join(var_map, by = "index")
rhs_draw  <- rhs_draw  %>% left_join(var_map, by = "index")
bugs_var  <- bugs_var  %>% left_join(var_map, by = "index")
rhs_var   <- rhs_var   %>% left_join(var_map, by = "index")

# -----------------------------
# Color palette
# -----------------------------
col_guided_true    <- "#B22222"  # deep red
col_guided_false   <- "#F4A6A6"  # light red
col_unguided_true  <- "#1F4E79"  # deep blue
col_unguided_false <- "#A9C8E8"  # light blue
col_marg_true      <- "#5C3D99"  # deep purple
col_marg_false     <- "#CDB7F6"  # light purple

fill_truth_marg <- c(
  "True variable"  = col_marg_true,
  "Noise variable" = col_marg_false
)

fill_model_truth <- c(
  "Guided | True variable"    = col_guided_true,
  "Guided | Noise variable"   = col_guided_false,
  "Unguided | True variable"  = col_unguided_true,
  "Unguided | Noise variable" = col_unguided_false
)

# -----------------------------
# Theme
# -----------------------------
theme_paper <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10.5, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold", color = "black"),
      axis.text.y = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      strip.background = element_rect(fill = "grey95", color = "grey60"),
      strip.text = element_text(face = "bold")
    )
}

# -----------------------------
# 1) Intro panel 1: marginal scores
# -----------------------------
# Marginal-score legend is placed inside the plot (top-right).
p_marg <- ggplot(scores, aes(x = var, y = marginal_score, fill = truth_label)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  scale_fill_manual(
    values = fill_truth_marg,
    name = NULL,
    breaks = c("True variable", "Noise variable")
  ) +
  labs(
    title = "Marginal scores",
    x = "Variable",
    y = expression(s[j] == frac("|" * x[j]^T * y * "|", n))
  ) +
  theme_paper() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = scales::alpha("white", 0.9), color = "grey40"),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size = 12),
    legend.margin = margin(3, 4, 3, 4)
  )

# -----------------------------
# 2) Intro panel 2: boxplots of posterior local scale
# -----------------------------
shrink_df <- bind_rows(
  bugs_draw %>%
    select(index, var, truth_label, value = local_scale) %>%
    mutate(model = "Guided"),
  rhs_draw %>%
    select(index, var, truth_label, value = local_scale) %>%
    mutate(model = "Unguided")
) %>%
  mutate(model_truth = paste(model, truth_label, sep = " | "))

# Remove non-positive values so log scale is valid
shrink_df <- shrink_df %>%
  filter(is.finite(value), value > 0)

p_shrink <- ggplot(shrink_df, aes(x = var, y = value, fill = model_truth)) +
  geom_boxplot(
    position = position_dodge(width = 0.78),
    width = 0.68,
    outlier.size = 0.30,
    outlier.alpha = 0.22,
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = fill_model_truth,
    breaks = c("Guided | True variable",
               "Guided | Noise variable",
               "Unguided | True variable",
               "Unguided | Noise variable"),
    labels = c("Guided: true variable", "Guided: noise variable",
               "Unguided: true variable", "Unguided: noise variable"),
    name = NULL
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  coord_cartesian(ylim = c(local_scale_min, local_scale_max)) +
  annotation_logticks(sides = "l") +
  labs(
    title = "Posterior local-scale distributions (boxplots)",
    x = "Variable",
    y = "Local scale draws"
  ) +
  theme_paper()

# -----------------------------
# 3) Intro panel 3: posterior exceedance probabilities
# -----------------------------
incl_df <- bind_rows(
  bugs_var %>%
    select(index, var, truth_label, post_prob_abs_gt_thresh) %>%
    mutate(model = "Guided"),
  rhs_var %>%
    select(index, var, truth_label, post_prob_abs_gt_thresh) %>%
    mutate(model = "Unguided")
) %>%
  mutate(model_truth = paste(model, truth_label, sep = " | "))

p_incl <- ggplot(incl_df, aes(x = var, y = post_prob_abs_gt_thresh, fill = model_truth)) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.7,
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = fill_model_truth,
    breaks = c("Guided | True variable",
               "Guided | Noise variable",
               "Unguided | True variable",
               "Unguided | Noise variable"),
    labels = c("Guided: true variable", "Guided: noise variable",
               "Unguided: true variable", "Unguided: noise variable"),
    name = NULL
  ) +
  scale_y_continuous(
    breaks = seq(0.2, 1, by = 0.2),
    labels = number_format(accuracy = 0.1),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0.2, 1)) +
  labs(
    title = "Posterior exceedance probability",
    subtitle = bquote(P(group("|", beta[j], "|") > .(sel_thresh) ~ "|" ~ y)),
    x = "Variable",
    y = "Posterior probability"
  ) +
  theme_paper()

# -----------------------------
# Combined 3x1 introduction panel
# -----------------------------
# Shared Guided/Unguided legend at the top
leg_model <- cowplot::get_legend(
  p_incl +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "center",
      legend.spacing.x = unit(0.6, "lines"),
      legend.key.width = unit(1.4, "lines"),
      legend.key.size = unit(0.6, "cm"),
      legend.text = element_text(size = 12, face = "bold"),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
      legend.box.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
      legend.margin = margin(t = 4, r = 6, b = 4, l = 6)
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
)

# Main stacked panel:
# - p_marg keeps its own legend inside the plot
# - p_shrink and p_incl use the shared top legend
p_intro_core <- (
  p_marg
) / (
  p_shrink + theme(legend.position = "none")
) / (
  p_incl + theme(legend.position = "none")
) +
  plot_layout(heights = c(1.0, 1.15, 1.05)) &
  theme(plot.margin = margin(t = 6, r = 5, b = 6, l = 5))

# Stack shared legend above the full figure
p_intro_final <- cowplot::plot_grid(
  leg_model,
  p_intro_core,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

p_intro_final

# -----------------------------
# Save final combined figure
# -----------------------------
ggsave(
  file.path(out_dir, sprintf("Figure_intro_3x1_panel_n_%d_p_%d_rep_%d.jpg", n, p, rep)),
  p_intro_final, width = 10.5, height = 11.0, dpi = 600
)

ggsave(
  file.path(out_dir, sprintf("Figure_intro_3x1_panel_n_%d_p_%d_rep_%d.pdf", n, p, rep)),
  p_intro_final, width = 10.5, height = 11.0, device = cairo_pdf
)

# -----------------------------
# Export figure data used in the combined figure
# -----------------------------
write_csv(incl_df, file.path(out_dir, sprintf("FigureData_intro_exceedance_n_%d_p_%d_rep_%d.csv", n, p, rep)))
write_csv(shrink_df, file.path(out_dir, sprintf("FigureData_intro_localscale_n_%d_p_%d_rep_%d.csv", n, p, rep)))

message("Final combined figure saved in: ", normalizePath(out_dir))