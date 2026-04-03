# ============================================================
# Companion plots for top CpGs
# 1. Selection-strength plot
# 2. Threshold sensitivity path plot
# Consistent with forest-plot settings
# ============================================================

setwd("U:/BUGS/Case study")
rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(scales)

# ----------------------------
# User settings
# ----------------------------
display_top_many <- 10
save_top_many <- 10
use_scale <- "standardized"   # "standardized" or "raw"

recompute_metrics <- 1        # 1 = recompute summaries from posterior draws after cutoff
obs_cutoff <- 5e-3            # draws with abs(beta_draw) < obs_cutoff are excluded if recompute_metrics == 1
effect_thresh <- 0.01         # threshold used for post_prob_abs_gt_thresh and selection flags

# threshold grid for sensitivity path
threshold_grid <- seq(0, 0.20, length.out = 100)

save_plot <- TRUE

# ----------------------------
# Paths
# ----------------------------
outdir <- "Output"
datadir <- "Age methylation data"

file_summary <- file.path(outdir, "Output_FULL_topK_beta_summary.csv")
file_draws   <- file.path(outdir, "Output_FULL_topK_beta_draws.csv")
file_names   <- file.path(datadir, "X_cpg_names.rds")

outfile_select <- file.path(
  outdir,
  paste0(
    "Output_R_selection_strength_top", save_top_many,
    "_recompute", recompute_metrics,
    "_cut", format(obs_cutoff, scientific = FALSE),
    "_thr", format(effect_thresh, scientific = FALSE),
    ".png"
  )
)

outfile_path <- file.path(
  outdir,
  paste0(
    "Output_R_threshold_sensitivity_top", save_top_many,
    "_recompute", recompute_metrics,
    "_cut", format(obs_cutoff, scientific = FALSE),
    "_thr", format(effect_thresh, scientific = FALSE),
    ".png"
  )
)

outfile_summary <- file.path(
  outdir,
  paste0(
    "Output_R_selection_strength_summary_top", save_top_many,
    "_recompute", recompute_metrics,
    "_cut", format(obs_cutoff, scientific = FALSE),
    "_thr", format(effect_thresh, scientific = FALSE),
    ".csv"
  )
)

# ----------------------------
# Read files
# ----------------------------
top_sum   <- read_csv(file_summary, show_col_types = FALSE)
top_draws <- read_csv(file_draws, show_col_types = FALSE)
cpg_names <- readRDS(file_names)

# ----------------------------
# Basic checks
# ----------------------------
stopifnot("rank" %in% names(top_sum))
stopifnot("index" %in% names(top_sum))
stopifnot(length(cpg_names) >= max(top_sum$index, na.rm = TRUE))

display_top_many <- min(display_top_many, nrow(top_sum))
save_top_many <- min(save_top_many, nrow(top_sum))
save_top_many <- max(save_top_many, display_top_many)

# ----------------------------
# Top CpGs from saved rank file
# ----------------------------
save_sum0 <- top_sum %>%
  dplyr::arrange(rank) %>%
  dplyr::slice(1:save_top_many) %>%
  dplyr::mutate(
    cpg_name = cpg_names[index]
  )

# ----------------------------
# CpG label map with rank
# ----------------------------
label_map <- save_sum0 %>%
  dplyr::mutate(
    cpg_label = paste0(cpg_name, " (rank ", rank, ")")
  ) %>%
  dplyr::select(cpg_name, cpg_label)

# ----------------------------
# Load posterior draws for save_top_many ranks
# ----------------------------
draw_cols <- paste0("beta_draw_top_rank_", save_sum0$rank)
missing_cols <- setdiff(draw_cols, names(top_draws))

if (length(missing_cols) > 0) {
  stop("Missing draw columns in topK draws file: ",
       paste(missing_cols, collapse = ", "))
}

save_draws <- top_draws %>%
  dplyr::select(all_of(draw_cols)) %>%
  stats::setNames(as.character(save_sum0$cpg_name)) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "cpg_name",
    values_to = "beta_draw"
  )

# ----------------------------
# Convert draws to raw scale if requested
# ----------------------------
if (use_scale == "raw") {
  scale_map <- save_sum0 %>%
    dplyr::mutate(
      scale_factor = ifelse(abs(beta_mean) > 0,
                            beta_mean_raw / beta_mean,
                            NA_real_)
    ) %>%
    dplyr::select(cpg_name, scale_factor)
  
  idx_bad <- which(is.na(scale_map$scale_factor) | !is.finite(scale_map$scale_factor))
  
  if (length(idx_bad) > 0) {
    alt <- save_sum0 %>%
      dplyr::mutate(
        scale_factor_alt = ifelse(abs(beta_median) > 0,
                                  beta_median_raw / beta_median,
                                  NA_real_)
      ) %>%
      dplyr::pull(scale_factor_alt)
    scale_map$scale_factor[idx_bad] <- alt[idx_bad]
  }
  
  save_draws <- save_draws %>%
    dplyr::left_join(scale_map, by = "cpg_name") %>%
    dplyr::mutate(beta_draw = beta_draw * scale_factor) %>%
    dplyr::select(cpg_name, beta_draw)
}

# ----------------------------
# Apply obs_cutoff if recompute_metrics == 1
# ----------------------------
if (recompute_metrics == 1) {
  save_draws_used <- save_draws %>%
    dplyr::filter(!is.na(beta_draw), abs(beta_draw) >= obs_cutoff)
} else {
  save_draws_used <- save_draws %>%
    dplyr::filter(!is.na(beta_draw))
}

# ----------------------------
# Summaries from posterior draws
# ----------------------------
draw_summary <- save_draws_used %>%
  dplyr::group_by(cpg_name) %>%
  dplyr::summarise(
    n_draw_used = dplyr::n(),
    beta_mean = mean(beta_draw),
    beta_median = median(beta_draw),
    beta_sd = ifelse(dplyr::n() > 1, sd(beta_draw), 0),
    ci_2p5 = as.numeric(quantile(beta_draw, 0.025, names = FALSE, type = 7)),
    ci_97p5 = as.numeric(quantile(beta_draw, 0.975, names = FALSE, type = 7)),
    post_prob_abs_gt_thresh = mean(abs(beta_draw) > effect_thresh),
    post_prob_positive = mean(beta_draw > 0),
    sign_mean = dplyr::case_when(
      beta_mean > 0 ~ "Positive",
      beta_mean < 0 ~ "Negative",
      TRUE ~ "Near zero"
    ),
    .groups = "drop"
  )

plot_dat <- save_sum0 %>%
  dplyr::select(rank, index, cpg_name) %>%
  dplyr::left_join(draw_summary, by = "cpg_name") %>%
  dplyr::mutate(
    cpg_name = factor(as.character(cpg_name), levels = rev(as.character(cpg_name)))
  )

plot_dat$sign_mean <- factor(
  plot_dat$sign_mean,
  levels = c("Positive", "Negative", "Near zero")
)

# save summary table
readr::write_csv(
  plot_dat %>%
    dplyr::mutate(cpg_name = as.character(cpg_name)),
  outfile_summary
)

# ============================================================
# 1. Selection-strength plot
# ============================================================

p_select <- ggplot(
  plot_dat,
  aes(x = post_prob_abs_gt_thresh, y = cpg_name, color = sign_mean)
) +
  geom_segment(
    aes(x = 0.5, xend = post_prob_abs_gt_thresh, yend = cpg_name),
    linewidth = 1.0,
    color = "grey60"
  ) +
  geom_point(
    size = 3.5,
    alpha = 0.95
  ) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey65"
  ) +
  scale_color_manual(
    values = c(
      "Positive" = "#C65A46",
      "Negative" = "#3B76AF",
      "Near zero" = "#9A9A9A"
    ),
    name = "Posterior\nmean sign"
  ) +
  scale_x_continuous(
    limits = c(.5, 1),
    breaks = c(0.5, 0.6, 0.7, 0.8, 0.9 , 1),
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_discrete(
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = bquote(P( "|" * beta * "|" > .(effect_thresh) )),
    y = "CpG"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.4),
    panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_text(size = 13, color = "black", face = "bold"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.85),
      color = "black",        # adds border
      linewidth = 0.6         # border thickness
    ),
    plot.margin = margin(10, 14, 8, 8)
  )

print(p_select)

if (save_plot) {
  ggsave(
    filename = outfile_select,
    plot = p_select,
    width = 7.9,
    height = 6.4, #max(6.4, 0.52 * save_top_many + 1.2),
    dpi = 500,
    bg = "white"
  )
}

# ============================================================
# 2. Threshold sensitivity path plot
#    P(|beta| > t) as threshold t varies
# ============================================================

draw_list <- save_draws_used %>%
  dplyr::group_by(cpg_name) %>%
  dplyr::summarise(draws = list(beta_draw), .groups = "drop")

threshold_path <- tidyr::crossing(
  cpg_name = unique(draw_list$cpg_name),
  threshold = threshold_grid
) %>%
  dplyr::left_join(draw_list, by = "cpg_name") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    prob_selected = mean(abs(unlist(draws)) > threshold, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-draws) %>%
  dplyr::left_join(
    draw_summary %>% dplyr::select(cpg_name, beta_mean),
    by = "cpg_name"
  ) %>%
  dplyr::left_join(label_map, by = "cpg_name")

threshold_path$cpg_label <- factor(
  threshold_path$cpg_label,
  levels = label_map$cpg_label
)

p_path <- ggplot(
  threshold_path,
  aes(x = threshold, y = prob_selected, group = cpg_label, color = cpg_label)
) +
  geom_line(linewidth = 1.2, alpha = 0.95) +
  scale_color_manual(
    values = scales::hue_pal()(length(unique(threshold_path$cpg_label))),
    name = "Top 10 CpGs"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_continuous(
    labels = scales::label_number(accuracy = 0.001)
  ) +
  labs(
    x = "Selection threshold",
    y = expression(paste("Posterior probability (| ", beta, "| > threshold)"))
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey93", linewidth = 0.4),
    panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
    axis.line = element_blank(),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    
    # move legend to top-right inside
    legend.position = c(0.99, 0.96),
    legend.justification = c(1, 1),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.85),
      color = "black",        # adds border
      linewidth = 0.6         # border thickness
    ),
    
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    
    plot.margin = margin(10, 14, 8, 8)
  )

print(p_path)

if (save_plot) {
  ggsave(
    filename = outfile_path,
    plot = p_path,
    width = 7.9,
    height = 6.4,
    dpi = 500,
    bg = "white"
  )
}

message("Saved selection-strength plot to: ", outfile_select)
message("Saved threshold-sensitivity plot to: ", outfile_path)
message("Saved summary table to: ", outfile_summary)