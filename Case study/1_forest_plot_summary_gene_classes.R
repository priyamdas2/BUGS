# ============================================================
# Forest plot with posterior density curves for top CpGs
# + annotation table (gene / location) for displayed/saved CpGs
# Robust EPIC annotation loading
# Optional recomputation of summaries from posterior draws
# ============================================================

setwd("U:/BUGS/Case study")
rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(tibble)

# ----------------------------
# User settings
# ----------------------------
display_top_many <- 10
save_top_many <- 10
use_scale <- "standardized"   # "standardized" or "raw"
save_plot <- TRUE
save_annotation_table <- TRUE

recompute_metrics <- 1        # 1 = recompute summaries from posterior draws after cutoff
obs_cutoff <- 5e-3            # draws with abs(beta_draw) < obs_cutoff are excluded if recompute_metrics == 1
effect_thresh <- 0.01         # threshold used for post_prob_abs_gt_thresh and selection flags

# ----------------------------
# Paths
# ----------------------------
outdir <- "Output"
datadir <- "Age methylation data"

file_summary <- file.path(outdir, "Output_FULL_topK_beta_summary.csv")
file_draws   <- file.path(outdir, "Output_FULL_topK_beta_draws.csv")
file_names   <- file.path(datadir, "X_cpg_names.rds")

outfile_plot  <- file.path(
  outdir,
  paste0(
    "Output_R_forest_density_top", display_top_many,
    "_save", save_top_many,
    "_recompute", recompute_metrics,
    "_cut", format(obs_cutoff, scientific = FALSE),
    "_CpGs.png"
  )
)

outfile_annot <- file.path(
  outdir,
  paste0(
    "Output_R_top", save_top_many,
    "_display", display_top_many,
    "_recompute", recompute_metrics,
    "_cut", format(obs_cutoff, scientific = FALSE),
    "_CpGs_annotation.csv"
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
# Keep ranking from saved file
# ----------------------------
save_sum0 <- top_sum %>%
  dplyr::arrange(rank) %>%
  dplyr::slice(1:save_top_many) %>%
  dplyr::mutate(
    cpg_name = cpg_names[index],
    row_id = dplyr::row_number()
  )

# ----------------------------
# Load posterior draws for save_top_many ranks
# ----------------------------
draw_cols_save <- paste0("beta_draw_top_rank_", save_sum0$rank)
missing_cols <- setdiff(draw_cols_save, names(top_draws))

if (length(missing_cols) > 0) {
  stop("Missing draw columns in topK draws file: ",
       paste(missing_cols, collapse = ", "))
}

save_draws <- top_draws %>%
  dplyr::select(all_of(draw_cols_save)) %>%
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
  scale_map_save <- save_sum0 %>%
    dplyr::mutate(
      scale_factor = ifelse(abs(beta_mean) > 0,
                            beta_mean_raw / beta_mean,
                            NA_real_)
    ) %>%
    dplyr::select(cpg_name, scale_factor)
  
  idx_bad <- which(is.na(scale_map_save$scale_factor) | !is.finite(scale_map_save$scale_factor))
  
  if (length(idx_bad) > 0) {
    alt <- save_sum0 %>%
      dplyr::mutate(
        scale_factor_alt = ifelse(abs(beta_median) > 0,
                                  beta_median_raw / beta_median,
                                  NA_real_)
      ) %>%
      dplyr::pull(scale_factor_alt)
    scale_map_save$scale_factor[idx_bad] <- alt[idx_bad]
  }
  
  save_draws <- save_draws %>%
    dplyr::left_join(scale_map_save, by = "cpg_name") %>%
    dplyr::mutate(beta_draw = beta_draw * scale_factor) %>%
    dplyr::select(cpg_name, beta_draw)
}

# ----------------------------
# Build summaries
# ----------------------------
if (recompute_metrics == 1) {
  
  save_draws_used <- save_draws %>%
    dplyr::filter(!is.na(beta_draw), abs(beta_draw) >= obs_cutoff)
  
  recomputed_sum <- save_draws_used %>%
    dplyr::group_by(cpg_name) %>%
    dplyr::summarise(
      n_draw_used = dplyr::n(),
      beta_mean = mean(beta_draw),
      beta_median = stats::median(beta_draw),
      beta_sd = ifelse(dplyr::n() > 1, stats::sd(beta_draw), 0),
      ci_2p5 = as.numeric(stats::quantile(beta_draw, 0.025, names = FALSE, type = 7)),
      ci_97p5 = as.numeric(stats::quantile(beta_draw, 0.975, names = FALSE, type = 7)),
      post_prob_abs_gt_thresh = mean(abs(beta_draw) > effect_thresh),
      selected_mean_thresh = as.integer(abs(mean(beta_draw)) > effect_thresh),
      selected_median_thresh = as.integer(abs(stats::median(beta_draw)) > effect_thresh),
      selected_ci_excludes_zero = as.integer(
        (as.numeric(stats::quantile(beta_draw, 0.025, names = FALSE, type = 7)) > 0) |
          (as.numeric(stats::quantile(beta_draw, 0.975, names = FALSE, type = 7)) < 0)
      ),
      sign_beta_mean = dplyr::case_when(
        mean(beta_draw) > 0 ~ "positive",
        mean(beta_draw) < 0 ~ "negative",
        TRUE ~ "zero"
      ),
      sign_beta_median = dplyr::case_when(
        stats::median(beta_draw) > 0 ~ "positive",
        stats::median(beta_draw) < 0 ~ "negative",
        TRUE ~ "zero"
      ),
      .groups = "drop"
    )
  
  save_sum <- save_sum0 %>%
    dplyr::mutate(cpg_name = as.character(cpg_name)) %>%
    dplyr::left_join(recomputed_sum, by = "cpg_name", suffix = c(".orig", ""))
  
  # x-scale for plot from recomputed summaries
  save_sum <- save_sum %>%
    dplyr::mutate(
      x_mean = beta_mean,
      x_low  = ci_2p5,
      x_high = ci_97p5
    )
  
} else {
  
  save_draws_used <- save_draws
  
  save_sum <- save_sum0
  
  if (use_scale == "standardized") {
    save_sum <- save_sum %>%
      dplyr::mutate(
        x_mean = beta_mean,
        x_low  = ci_2p5,
        x_high = ci_97p5
      )
  } else if (use_scale == "raw") {
    save_sum <- save_sum %>%
      dplyr::mutate(
        x_mean = beta_mean_raw,
        x_low  = ci_2p5_raw,
        x_high = ci_97p5_raw
      )
  }
}

# ----------------------------
# Plot subset
# ----------------------------
plot_sum <- save_sum %>%
  dplyr::slice(1:display_top_many)

plot_draws <- save_draws_used %>%
  dplyr::filter(cpg_name %in% as.character(plot_sum$cpg_name))

# ----------------------------
# Preserve row order on y-axis
# ----------------------------
plot_sum <- plot_sum %>%
  dplyr::mutate(cpg_name = factor(cpg_name, levels = rev(as.character(cpg_name))))

plot_draws <- plot_draws %>%
  dplyr::mutate(cpg_name = factor(cpg_name, levels = levels(plot_sum$cpg_name)))

xlab_txt <- if (use_scale == "standardized") {
  "Posterior coefficient (standardized scale)"
} else {
  "Posterior coefficient (raw scale)"
}

# ----------------------------
# Robust helper to fetch EPIC annotation object
# ----------------------------
load_epic_anno_object <- function() {
  pkg <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed.")
  }
  if (!requireNamespace("minfi", quietly = TRUE)) {
    stop("Package 'minfi' is not installed.")
  }
  
  ns <- asNamespace(pkg)
  obj_name <- pkg
  
  if (exists(obj_name, envir = ns, inherits = FALSE)) {
    return(get(obj_name, envir = ns, inherits = FALSE))
  }
  
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
  if (exists(obj_name, inherits = TRUE)) {
    return(get(obj_name, inherits = TRUE))
  }
  
  try(data(list = obj_name, package = pkg, envir = environment()), silent = TRUE)
  if (exists(obj_name, inherits = FALSE)) {
    return(get(obj_name, inherits = FALSE))
  }
  
  stop("Could not locate annotation object '", obj_name,
       "' inside package '", pkg, "'.")
}

# ----------------------------
# Annotation helper
# ----------------------------
get_epic_annotation <- function(cpg_ids) {
  out <- tibble(
    cpg_name = cpg_ids,
    chr = NA_character_,
    pos = NA_real_,
    strand = NA_character_,
    gene = NA_character_,
    gene_group = NA_character_,
    island_relation = NA_character_
  )
  
  if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
    message("EPIC annotation package not found. Returning NA annotation columns.")
    return(out)
  }
  
  if (!requireNamespace("minfi", quietly = TRUE)) {
    message("Package 'minfi' not found. Returning NA annotation columns.")
    return(out)
  }
  
  ann_pkg_obj <- tryCatch(
    load_epic_anno_object(),
    error = function(e) {
      message("Could not load EPIC annotation object: ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (is.null(ann_pkg_obj)) {
    return(out)
  }
  
  ann_obj <- tryCatch(
    minfi::getAnnotation(ann_pkg_obj),
    error = function(e) {
      message("minfi::getAnnotation failed: ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (is.null(ann_obj)) {
    return(out)
  }
  
  ann_df <- as.data.frame(ann_obj)
  ann_df$Name <- rownames(ann_df)
  
  keep <- ann_df %>%
    dplyr::filter(Name %in% cpg_ids)
  
  if (nrow(keep) == 0) {
    message("No CpG IDs matched the EPIC annotation object.")
    return(out)
  }
  
  chr_col <- intersect(c("chr", "CHR"), names(keep))
  pos_col <- intersect(c("pos", "MAPINFO"), names(keep))
  strand_col <- intersect(c("Strand", "strand"), names(keep))
  gene_col <- intersect(c("UCSC_RefGene_Name"), names(keep))
  gene_group_col <- intersect(c("UCSC_RefGene_Group"), names(keep))
  island_col <- intersect(c("Relation_to_Island", "Relation_to_UCSC_CpG_Island"), names(keep))
  
  keep2 <- tibble(
    cpg_name = keep$Name,
    chr = if (length(chr_col) > 0) as.character(keep[[chr_col[1]]]) else NA_character_,
    pos = if (length(pos_col) > 0) suppressWarnings(as.numeric(keep[[pos_col[1]]])) else NA_real_,
    strand = if (length(strand_col) > 0) as.character(keep[[strand_col[1]]]) else NA_character_,
    gene = if (length(gene_col) > 0) as.character(keep[[gene_col[1]]]) else NA_character_,
    gene_group = if (length(gene_group_col) > 0) as.character(keep[[gene_group_col[1]]]) else NA_character_,
    island_relation = if (length(island_col) > 0) as.character(keep[[island_col[1]]]) else NA_character_
  )
  
  out %>%
    dplyr::left_join(keep2, by = "cpg_name", suffix = c("", ".ann")) %>%
    dplyr::transmute(
      cpg_name = cpg_name,
      chr = dplyr::coalesce(chr.ann, chr),
      pos = dplyr::coalesce(pos.ann, pos),
      strand = dplyr::coalesce(strand.ann, strand),
      gene = dplyr::coalesce(gene.ann, gene),
      gene_group = dplyr::coalesce(gene_group.ann, gene_group),
      island_relation = dplyr::coalesce(island_relation.ann, island_relation)
    )
}

# ----------------------------
# Build annotation table for save_top_many
# ----------------------------
annot_tbl <- get_epic_annotation(as.character(save_sum$cpg_name))

if (recompute_metrics == 1) {
  display_tbl <- save_sum %>%
    dplyr::mutate(cpg_name = as.character(cpg_name)) %>%
    dplyr::left_join(annot_tbl, by = "cpg_name") %>%
    dplyr::select(
      rank,
      index,
      cpg_name,
      gene,
      gene_group,
      chr,
      pos,
      strand,
      island_relation,
      beta_mean,
      beta_median,
      beta_sd,
      ci_2p5,
      ci_97p5,
      post_prob_abs_gt_thresh,
      selected_mean_thresh,
      selected_median_thresh,
      selected_ci_excludes_zero,
      sign_beta_mean,
      sign_beta_median,
      n_draw_used
    )
} else {
  display_tbl <- save_sum %>%
    dplyr::mutate(cpg_name = as.character(cpg_name)) %>%
    dplyr::left_join(annot_tbl, by = "cpg_name") %>%
    dplyr::select(
      rank,
      index,
      cpg_name,
      gene,
      gene_group,
      chr,
      pos,
      strand,
      island_relation,
      beta_mean,
      beta_median,
      beta_sd,
      ci_2p5,
      ci_97p5,
      beta_mean_raw,
      beta_median_raw,
      ci_2p5_raw,
      ci_97p5_raw,
      post_prob_abs_gt_thresh,
      selected_mean_thresh,
      selected_median_thresh,
      selected_ci_excludes_zero,
      sign_beta_mean,
      sign_beta_median
    )
}

if (save_annotation_table) {
  readr::write_csv(display_tbl, outfile_annot)
}

# ----------------------------
# Plot: publication-style with meaningful color
# ----------------------------

# Ensure plot_sum has numeric x_mean / x_low / x_high already defined
plot_sum_color <- plot_sum %>%
  dplyr::mutate(
    cpg_name_chr = as.character(cpg_name)
  ) %>%
  dplyr::select(cpg_name_chr, x_mean)

plot_draws2 <- plot_draws %>%
  dplyr::mutate(cpg_name_chr = as.character(cpg_name)) %>%
  dplyr::left_join(plot_sum_color, by = "cpg_name_chr")

x_min <- min(
  plot_draws2$beta_draw,
  plot_sum$x_low,
  na.rm = TRUE
)

x_max <- max(
  plot_draws2$beta_draw,
  plot_sum$x_high,
  na.rm = TRUE
)

x_pad <- 0.06 * (x_max - x_min)
x_limits <- c(x_min - x_pad, x_max + x_pad)

fill_lim <- max(abs(plot_sum$x_mean), na.rm = TRUE)

p <- ggplot() +
  geom_density_ridges(
    data = plot_draws2,
    aes(x = beta_draw, y = cpg_name, fill = x_mean),
    scale = 0.68,
    rel_min_height = 0.008,
    color = "grey35",
    linewidth = 0.45,
    alpha = 0.9
  ) +
  geom_segment(
    data = plot_sum,
    aes(x = x_low, xend = x_high, y = cpg_name, yend = cpg_name),
    linewidth = 1.0,
    color = "black",
    lineend = "butt"
  ) +
  geom_point(
    data = plot_sum,
    aes(x = x_mean, y = cpg_name),
    shape = 21,
    size = 2.9,
    stroke = 0.8,
    fill = "white",
    color = "black"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey55"
  ) +
  scale_fill_gradient2(
    low = "#3B76AF",
    mid = "#E9E9E9",
    high = "#C65A46",
    midpoint = 0,
    limits = c(-fill_lim, fill_lim),
    name = "Posterior\nmean"
  ) +
  scale_x_continuous(
    limits = x_limits,
    expand = expansion(mult = c(0, 0))
  )  +
  scale_y_discrete(
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(
    x = "Posterior coefficient (standardized scale)",
    y = "CpG"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.4),
    
    # FULL SQUARE BORDER
    panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
    
    # remove separate axis lines (avoid double borders)
    axis.line = element_blank(),
    
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    
    axis.text.y = element_text(size = 13, color = "black", face = "bold"),
    axis.text.x = element_text(size = 11, color = "black"),
    
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    
    # LEGEND IMPROVEMENTS
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.height = unit(1.2, "cm"),
    legend.key.width = unit(0.6, "cm"),
    
    plot.margin = margin(10, 14, 8, 8)
  )
print(p)

if (save_plot) {
  ggsave(
    filename = outfile_plot,
    plot = p,
    width = 7.9,
    height = 6.4, #max(5.9, 0.78 * display_top_many + 1.3),
    dpi = 500,
    bg = "white"
  )
}

message("Saved plot to: ", outfile_plot)
message("Saved annotation table to: ", outfile_annot)
print(display_tbl)