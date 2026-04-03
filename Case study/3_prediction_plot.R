# ============================================================
# Real-scale observed vs predicted plot using top CpGs
# Reads MATLAB-prepared files from:
#   ./Final prediction plot data/
#
# Prediction is based on top use_top_many posterior mean coefficients
# on the original response scale.
# ============================================================

# ----------------------------
# Working directory and clean start
# ----------------------------
setwd("U:/BUGS/Case study")
rm(list = ls())

# ----------------------------
# Required libraries
# ----------------------------
library(readr)
library(dplyr)
library(ggplot2)
library(scales)

# ----------------------------
# User settings
# ----------------------------
# use_top_many:
#   Number of top CpGs to use for prediction and plotting.
#
# savedFromMatlabK:
#   MATLAB had saved the top-K CpGs; this code reads those saved files.
#
# use_scale:
#   "standardized" -> use standardized-scale coefficients
#   "raw"          -> use raw-scale coefficients if available
#
# recompute_metrics:
#   1 -> use beta_mean from metadata (post-processed/recomputed setting)
#   0 -> use the scale chosen via use_scale
#
# obs_cutoff:
#   Used only for naming/output tracking here; corresponds to the cutoff
#   applied when recomputing summaries from posterior draws.
#
# effect_thresh:
#   Threshold used in associated metadata summaries and output naming.
use_top_many <- 10
savedFromMatlabK <- 200
use_scale <- "standardized"   # "standardized" or "raw"

recompute_metrics <- 1        # 1 = recompute summaries from posterior draws after cutoff
obs_cutoff <- 5e-3            # draws with abs(beta_draw) < obs_cutoff are excluded if recompute_metrics == 1
effect_thresh <- 0.01         # threshold used for post_prob_abs_gt_thresh and selection flags

# Output controls
save_plot <- TRUE
save_prediction_csv <- TRUE

# ----------------------------
# Paths
# ----------------------------
# Input directory prepared externally (from MATLAB side)
datadir <- "Final prediction plot data"

# Output directory for R-generated plots and CSV files
outdir  <- "Output"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Input files
file_X    <- file.path(datadir, paste0("selected_X_raw_top", savedFromMatlabK, ".csv"))
file_y    <- file.path(datadir, "y_raw.csv")
file_ysum <- file.path(datadir, "rawY_mean_sd.csv")
file_meta <- file.path(datadir, paste0("selected_top", savedFromMatlabK, "_metadata.csv"))

# Output file: scatter plot
outfile_plot <- file.path(
  outdir,
  paste0(
    "Output_R_observed_vs_predicted_top", use_top_many,
    "_fromK", savedFromMatlabK,
    "_", use_scale,
    "_recompute", recompute_metrics,
    "_cut", format(obs_cutoff, scientific = FALSE),
    "_thr", format(effect_thresh, scientific = FALSE),
    ".png"
  )
)

# Output file: prediction CSV
outfile_pred <- file.path(
  outdir,
  paste0(
    "Output_R_observed_vs_predicted_top", use_top_many,
    "_fromK", savedFromMatlabK,
    "_", use_scale,
    "_recompute", recompute_metrics,
    "_cut", format(obs_cutoff, scientific = FALSE),
    "_thr", format(effect_thresh, scientific = FALSE),
    ".csv"
  )
)

# ----------------------------
# Read files
# ----------------------------
# X_raw_all:
#   Raw methylation/CpG design matrix for the selected top-K CpGs
#
# y_raw:
#   Response on original/raw scale (here age in months)
#
# y_info:
#   Contains rawY_mean and rawY_sd used to map between standardized and raw scale
#
# meta_all:
#   Metadata for the top-K CpGs, including rank, original index,
#   scaling info, and coefficient summaries
X_raw_all <- as.matrix(read_csv(file_X, col_names = FALSE, show_col_types = FALSE))
y_raw <- read_csv(file_y, col_names = FALSE, show_col_types = FALSE)[[1]]
y_info <- read_csv(file_ysum, show_col_types = FALSE)
meta_all <- read_csv(file_meta, show_col_types = FALSE)

# ----------------------------
# Basic checks
# ----------------------------
# Confirm key columns exist in metadata
stopifnot("rank" %in% names(meta_all))
stopifnot("index" %in% names(meta_all))
stopifnot("rawX_mean" %in% names(meta_all))
stopifnot("rawX_sd" %in% names(meta_all))
stopifnot("rawY_mean" %in% names(y_info) || all(c("rawY_mean", "rawY_sd") %in% names(y_info)))

# Require both mean and sd for Y scaling back to raw scale
if (!all(c("rawY_mean", "rawY_sd") %in% names(y_info))) {
  stop("rawY_mean_sd.csv must contain columns rawY_mean and rawY_sd.")
}

# Check sample size consistency between X and y
if (nrow(X_raw_all) != length(y_raw)) {
  stop("Number of rows in selected_X_raw does not match length of y_raw.")
}

# Sort metadata by rank so that top rows correspond to top-ranked CpGs
meta_all <- meta_all %>%
  dplyr::arrange(rank)

# Ensure requested number of CpGs does not exceed available number
use_top_many <- min(use_top_many, nrow(meta_all))

# Extract only the top use_top_many CpGs
X_raw <- X_raw_all[, 1:use_top_many, drop = FALSE]
meta_use <- meta_all %>%
  dplyr::slice(1:use_top_many)

# Y scaling information
muY <- y_info$rawY_mean[1]
sdY <- y_info$rawY_sd[1]

# X scaling information for the selected CpGs
muX <- meta_use$rawX_mean
sdX <- meta_use$rawX_sd

# Protect against division by zero / missing SD
sdX[is.na(sdX) | sdX == 0] <- 1

# ----------------------------
# Choose coefficient source
# ----------------------------
# Cases:
#   (i) recompute_metrics == 1:
#       use beta_mean from metadata
#   (ii) recompute_metrics == 0 and use_scale == "standardized":
#       use beta_mean
#   (iii) recompute_metrics == 0 and use_scale == "raw":
#       use beta_mean_raw
if (recompute_metrics == 1) {
  if (!("beta_mean" %in% names(meta_use))) {
    stop("Metadata file does not contain beta_mean needed for recompute_metrics = 1.")
  }
  beta_std <- meta_use$beta_mean
} else {
  if (use_scale == "standardized") {
    if (!("beta_mean" %in% names(meta_use))) {
      stop("Metadata file does not contain beta_mean.")
    }
    beta_std <- meta_use$beta_mean
  } else if (use_scale == "raw") {
    if (!("beta_mean_raw" %in% names(meta_use))) {
      stop("Metadata file does not contain beta_mean_raw.")
    }
    beta_raw <- meta_use$beta_mean_raw
  } else {
    stop("use_scale must be either 'standardized' or 'raw'.")
  }
}

# ----------------------------
# Build predictions on real scale
# ----------------------------
# If coefficients are standardized (or recomputed in that form):
#   1. Standardize X using raw means/SDs
#   2. Compute standardized linear predictor
#   3. Transform prediction back to raw Y scale
#   4. Derive equivalent raw-scale coefficients and implied intercept
#
# If raw-scale coefficients are used:
#   1. Center X by raw means
#   2. Multiply by raw coefficients
#   3. Add back implied raw intercept via muY
if (use_scale == "standardized" || recompute_metrics == 1) {
  X_std <- sweep(X_raw, 2, muX, FUN = "-")
  X_std <- sweep(X_std, 2, sdX, FUN = "/")
  
  yhat_std <- as.vector(X_std %*% beta_std)
  yhat_raw <- muY + sdY * yhat_std
  
  beta_raw_equiv <- (sdY / sdX) * beta_std
  intercept_raw <- muY - sum(muX * beta_raw_equiv)
} else {
  yhat_raw <- muY + as.vector(sweep(X_raw, 2, muX, FUN = "-") %*% beta_raw)
  intercept_raw <- muY - sum(muX * beta_raw)
  beta_raw_equiv <- beta_raw
}

# ----------------------------
# Performance summaries
# ----------------------------
# RMSE, MAE, correlation, and R^2 on raw/original scale
rmse <- sqrt(mean((y_raw - yhat_raw)^2, na.rm = TRUE))
mae <- mean(abs(y_raw - yhat_raw), na.rm = TRUE)
cor_val <- suppressWarnings(cor(y_raw, yhat_raw, use = "complete.obs"))
r2_val <- suppressWarnings(cor_val^2)

# Observation-level prediction table
pred_df <- tibble(
  obs_id = seq_along(y_raw),
  y_raw = y_raw,
  yhat_raw = yhat_raw,
  residual = y_raw - yhat_raw
)

# Coefficient table for the selected CpGs, augmented with raw-equivalent beta
coef_df <- meta_use %>%
  dplyr::mutate(
    beta_raw_equiv = beta_raw_equiv
  ) %>%
  dplyr::select(rank, index, everything())

# Save prediction and coefficient CSV files if requested
if (save_prediction_csv) {
  write_csv(pred_df, outfile_pred)
  write_csv(
    coef_df,
    file.path(
      outdir,
      paste0(
        "Output_R_prediction_coefficients_top", use_top_many,
        "_fromK", savedFromMatlabK,
        "_", use_scale,
        "_recompute", recompute_metrics,
        "_cut", format(obs_cutoff, scientific = FALSE),
        "_thr", format(effect_thresh, scientific = FALSE),
        ".csv"
      )
    )
  )
}

# ----------------------------
# Plot limits
# ----------------------------
# Shared axis range for observed vs predicted scatter plot
xy_min <- min(c(y_raw, yhat_raw), na.rm = TRUE)
xy_max <- max(c(y_raw, yhat_raw), na.rm = TRUE)
xy_pad <- 0.04 * (xy_max - xy_min)
xy_lim <- c(xy_min - xy_pad, xy_max + xy_pad)

# Annotation text printed inside plots
annot_txt <- paste0(
  "Top ", use_top_many, " CpGs",
  "\nRMSE = ", sprintf("%.3f", rmse),
  "\nMAE = ", sprintf("%.3f", mae),
  "\nCor = ", sprintf("%.3f", cor_val),
  "\nR^2 = ", sprintf("%.3f", r2_val)
)

# ----------------------------
# Plot 1: Scatter plot
# ----------------------------
# Observed raw age vs predicted raw age, with:
#   - 45-degree reference line
#   - in-panel summary statistics
#   - equal coordinate scaling
p <- ggplot(pred_df, aes(x = y_raw, y = yhat_raw)) +
  geom_point(
    shape = 21,
    size = 2.6,
    stroke = 0.5,
    fill = "#5B8CC0",
    color = "black",
    alpha = 0.75
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey45"
  ) +
  annotate(
    "text",
    x = xy_lim[1] + 0.06 * diff(xy_lim),
    y = xy_lim[2] - 0.08 * diff(xy_lim),
    label = annot_txt,
    hjust = 0,
    vjust = 1,
    size = 4.1
  ) +
  scale_x_continuous(
    limits = xy_lim,
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = xy_lim,
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "True age (months)",
    y = "Predicted age (months)"
  ) +
  coord_equal() +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
    axis.line = element_blank(),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    plot.margin = margin(10, 14, 8, 8)
  )

print(p)

# Save scatter plot if requested
if (save_plot) {
  ggsave(
    filename = outfile_plot,
    plot = p,
    width = 7.6,
    height = 7.0,
    dpi = 500,
    bg = "white"
  )
}

# ----------------------------
# Console output after scatter plot
# ----------------------------
message("Saved plot to: ", outfile_plot)
if (save_prediction_csv) {
  message("Saved predictions to: ", outfile_pred)
}
message("Raw-scale intercept implied by selected top CpGs: ", sprintf("%.6f", intercept_raw))
message("RMSE = ", sprintf("%.6f", rmse))
message("MAE  = ", sprintf("%.6f", mae))
message("Cor  = ", sprintf("%.6f", cor_val))
message("R^2  = ", sprintf("%.6f", r2_val))

# ----------------------------
# Libraries repeated in original code
# ----------------------------
# Kept exactly as in your original script.
library(ggplot2)
library(dplyr)
library(scales)

# ----------------------------
# Prepare data for violin plot
# ----------------------------
# pred_df already contains:
#   y_raw     = true/raw age
#   yhat_raw  = predicted/raw age
#
# For violin plot, convert true ages to factor so each observed age level
# gets its own violin/jitter group.
plot_df <- pred_df %>%
  dplyr::mutate(
    true_age = factor(y_raw, levels = sort(unique(y_raw)))
  )

# Numeric labels for x-axis, preserving age order
age_levels_num <- sort(unique(pred_df$y_raw))

# Recompute plot limits for the violin panel
xy_min <- min(c(pred_df$y_raw, pred_df$yhat_raw), na.rm = TRUE)
xy_max <- max(c(pred_df$y_raw, pred_df$yhat_raw), na.rm = TRUE)
xy_pad <- 0.04 * (xy_max - xy_min)
xy_lim <- c(xy_min - xy_pad, xy_max + xy_pad)

# Same annotation block reused in violin plot
annot_txt <- paste0(
  "Top ", use_top_many, " CpGs",
  "\nRMSE = ", sprintf("%.3f", rmse),
  "\nMAE = ", sprintf("%.3f", mae),
  "\nCor = ", sprintf("%.3f", cor_val),
  "\nR² = ", sprintf("%.3f", r2_val)
)

# ----------------------------
# Plot 2: Violin plot by true age
# ----------------------------
# Shows distribution of predicted ages within each true age group.
# Includes:
#   - violin density
#   - jittered individual predictions
#   - group mean marker
#
p <- ggplot(plot_df, aes(x = true_age, y = yhat_raw)) +
  geom_hline(
    yintercept = c(3, 9, 48, 72),
    linetype = "dotted",
    linewidth = 0.6,
    color = "grey40"
  ) +
  geom_violin(
    fill = "#D9E6F2",
    color = "grey40",
    linewidth = 0.5,
    alpha = 0.9,
    width = 0.8,
    trim = FALSE
  ) +
  geom_jitter(
    width = 0.08,
    height = 0,
    shape = 21,
    size = 2.3,
    stroke = 0.4,
    fill = "#5B8CC0",
    color = "black",
    alpha = 0.65
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3.8,
    stroke = 0.6,
    fill = "#C65A46",
    color = "black"
  ) +
  annotate(
    "label",
    x = 0.5,
    y = 7+xy_lim[2] - 0.08 * diff(xy_lim),
    label = annot_txt,
    hjust = 0,
    vjust = 1,
    size = 6.0,
    label.size = 0.99,           # border thickness
    label.padding = unit(0.25, "lines"),
    fill = "white",             # non-transparent background
    color = "black"
  ) +
  scale_x_discrete(
    name = "True age (months)",
    labels = age_levels_num
  ) +
  scale_y_continuous(
    limits = xy_lim,
    breaks = c(3, 9, 48, 72),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.02)),
    name = "Predicted age (months)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 15, margin = margin(t = 10)),
    axis.title.y = element_text(size = 15, margin = margin(r = 10)),
    plot.margin = margin(10, 14, 8, 8)
  )

print(p)

# Save violin plot if requested
if (save_plot) {
  ggsave(
    filename = file.path(
      outdir,
      paste0(
        "Output_R_observed_vs_predicted_violin_top", use_top_many,
        "_fromK", savedFromMatlabK,
        "_", use_scale,
        "_recompute", recompute_metrics,
        "_cut", format(obs_cutoff, scientific = FALSE),
        "_thr", format(effect_thresh, scientific = FALSE),
        ".png"
      )
    ),
    plot = p,
    width = 7.9,
    height = 6.4,
    dpi = 500,
    bg = "white"
  )
}