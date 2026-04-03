# ============================================================
# Create top-CpG summary table with coefficients on original scale
#
# Inputs:
#   1. Final prediction plot data/rawY_mean_sd.csv
#   2. Final prediction plot data/rawX_sd.csv
#   3. Output/Output_R_top10_display10_recompute1_cut0.005_CpGs_annotation.csv
#
# Output:
#   Output/Output_R_CpG_gene_map_summary.csv
#
# Notes:
#   - rawX_sd.csv has no header and contains one SD per CpG
#   - annot$index is assumed to be the 1-based row index into rawX_sd.csv
#   - beta and CI are back-transformed from standardized scale to raw scale
# ============================================================

# -----------------------------
# File paths
# -----------------------------
pred_dir <- "Final prediction plot data"
out_dir  <- "Output"

file_yinfo <- file.path(pred_dir, "rawY_mean_sd.csv")
file_xsd   <- file.path(pred_dir, "rawX_sd.csv")
file_annot <- file.path(out_dir,  "Output_R_top10_display10_recompute1_cut0.005_CpGs_annotation.csv")

outfile <- file.path(out_dir, "Output_R_CpG_gene_map_summary.csv")

# -----------------------------
# Read input files
# -----------------------------
y_info <- read.csv(file_yinfo, stringsAsFactors = FALSE)
annot  <- read.csv(file_annot, stringsAsFactors = FALSE)

# rawX_sd.csv has NO header
x_sd_df <- read.csv(file_xsd, header = FALSE, stringsAsFactors = FALSE)

# -----------------------------
# Basic checks
# -----------------------------
if (!file.exists(file_yinfo)) {
  stop("File not found: ", file_yinfo)
}
if (!file.exists(file_xsd)) {
  stop("File not found: ", file_xsd)
}
if (!file.exists(file_annot)) {
  stop("File not found: ", file_annot)
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

if (!all(c("rawY_mean", "rawY_sd") %in% names(y_info))) {
  stop("rawY_mean_sd.csv must contain columns 'rawY_mean' and 'rawY_sd'.")
}

required_cols <- c(
  "index", "cpg_name", "gene", "gene_group", "chr", "pos", "island_relation",
  "beta_mean", "ci_2p5", "ci_97p5"
)

missing_cols <- setdiff(required_cols, names(annot))
if (length(missing_cols) > 0) {
  stop(
    paste0(
      "Annotation file is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

if (ncol(x_sd_df) != 1) {
  stop("rawX_sd.csv should contain exactly one column and no header.")
}

# -----------------------------
# Extract SD information
# -----------------------------
raw_y_sd <- y_info$rawY_sd[1]
raw_x_sd <- x_sd_df[[1]]

if (!is.numeric(raw_y_sd) || length(raw_y_sd) != 1 || is.na(raw_y_sd) || raw_y_sd <= 0) {
  stop("rawY_sd must be a single positive numeric value.")
}

if (any(is.na(raw_x_sd)) || any(raw_x_sd <= 0)) {
  stop("rawX_sd.csv contains missing or non-positive SD values.")
}

# -----------------------------
# Match CpG indices to raw X SDs
# Assumes annot$index is 1-based row number in rawX_sd.csv
# -----------------------------
idx <- annot$index

if (any(is.na(idx))) {
  stop("Missing values found in the 'index' column of the annotation file.")
}

if (!is.numeric(idx)) {
  idx <- as.numeric(idx)
}

if (any(is.na(idx))) {
  stop("Could not coerce annot$index to numeric values.")
}

if (any(idx < 1 | idx > length(raw_x_sd))) {
  stop("Some annotation indices are out of bounds for rawX_sd.csv.")
}

x_sd_top <- raw_x_sd[idx]

# -----------------------------
# Back-transform standardized coefficients to raw scale
# beta_raw = beta_std * (sd_y / sd_xj)
# Same transformation applies to CI endpoints
# -----------------------------
scale_factor <- raw_y_sd / x_sd_top

beta_raw    <- annot$beta_mean * scale_factor
ci_2p5_raw  <- annot$ci_2p5    * scale_factor
ci_97p5_raw <- annot$ci_97p5   * scale_factor

# -----------------------------
# Create reporting fields
# -----------------------------
chr_pos   <- paste0(annot$chr, ":", annot$pos)
ci_raw_str <- sprintf("[%.1f, %.1f]", ci_2p5_raw, ci_97p5_raw)

# -----------------------------
# Final summary table
# -----------------------------
out_tab <- data.frame(
  CpG = annot$cpg_name,
  `Gene(s)` = annot$gene,
  Region = annot$gene_group,
  `Chr:Pos` = chr_pos,
  Island = annot$island_relation,
  beta_hat_actual_scale = round(beta_raw, 1),
  ci_2p5_actual_scale = round(ci_2p5_raw, 1),
  ci_97p5_actual_scale = round(ci_97p5_raw, 1),
  `95% CI` = ci_raw_str,
  stringsAsFactors = FALSE
)

# -----------------------------
# Save CSV
# -----------------------------
write.csv(out_tab, outfile, row.names = FALSE)

cat("Saved original-scale CpG summary table to:\n", outfile, "\n\n")
print(out_tab)