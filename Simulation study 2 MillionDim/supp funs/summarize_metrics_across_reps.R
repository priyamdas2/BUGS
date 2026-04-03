summarize_metrics_across_reps <- function(
    scenario,
    n,
    p,
    method_name,
    reps_to_summarize,
    input_dir = ".",
    output_dir = input_dir
) {
  # -----------------------------
  # Standard metric names
  # -----------------------------
  standard_metric_names <- c(
    "L2_beta", "L1_beta", "RMSE_beta", "PredMSE",
    "TPR", "FPR", "MCC", "FDR"
  )
  
  # -----------------------------
  # Helper: safely read one file
  # Keep only columns up to CompTime
  # -----------------------------
  read_one_file <- function(file_path) {
    dat <- tryCatch(
      read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(dat)) {
      warning(sprintf("Could not read file: %s", file_path))
      return(NULL)
    }
    
    if (ncol(dat) < 9) {
      warning(sprintf("File has fewer than 9 columns: %s", file_path))
      return(NULL)
    }
    
    # Standardize first 8 metric column names
    colnames(dat)[1:8] <- standard_metric_names
    
    # Standardize 9th column as CompTime
    colnames(dat)[9] <- "CompTime"
    
    # Keep only first 9 columns
    dat <- dat[, 1:9, drop = FALSE]
    
    return(dat)
  }
  
  # -----------------------------
  # Build exact file name
  # -----------------------------
  build_file_name <- function(rep_id) {
    sprintf(
      "Output_scen_%d_n_%s_p_%s_%s_rep_%s.csv",
      scenario, n, p, method_name, rep_id
    )
  }
  
  # -----------------------------
  # Locate files
  # -----------------------------
  file_list <- character(0)
  rep_found <- integer(0)
  
  for (rep_id in reps_to_summarize) {
    fname <- file.path(input_dir, build_file_name(rep_id))
    
    if (!file.exists(fname)) {
      warning(sprintf("Missing file: %s", fname))
    } else {
      file_list <- c(file_list, fname)
      rep_found <- c(rep_found, rep_id)
    }
  }
  
  if (length(file_list) == 0) {
    stop("No matching files were found.")
  }
  
  # -----------------------------
  # Read files
  # -----------------------------
  data_list <- lapply(file_list, read_one_file)
  keep_idx <- !vapply(data_list, is.null, logical(1))
  
  data_list <- data_list[keep_idx]
  file_list <- file_list[keep_idx]
  rep_found <- rep_found[keep_idx]
  
  if (length(data_list) == 0) {
    stop("No readable files were found.")
  }
  
  # -----------------------------
  # Check row counts match
  # -----------------------------
  row_counts <- vapply(data_list, nrow, integer(1))
  unique_row_counts <- unique(row_counts)
  
  if (length(unique_row_counts) != 1) {
    stop(sprintf(
      "Files do not all have the same number of rows. Row counts found: %s",
      paste(unique_row_counts, collapse = ", ")
    ))
  }
  
  n_rows_each <- unique_row_counts[1]
  
  # -----------------------------
  # Common summary columns
  # -----------------------------
  common_cols <- Reduce(intersect, lapply(data_list, colnames))
  
  preferred_cols <- c(
    "L2_beta", "L1_beta", "RMSE_beta", "PredMSE",
    "TPR", "FPR", "MCC", "FDR", "CompTime"
  )
  
  summary_cols <- intersect(preferred_cols, common_cols)
  
  if (length(summary_cols) == 0) {
    stop("No common summary columns found across files.")
  }
  
  # Convert summary columns to numeric
  for (j in seq_along(data_list)) {
    for (nm in summary_cols) {
      data_list[[j]][[nm]] <- suppressWarnings(as.numeric(data_list[[j]][[nm]]))
    }
  }
  
  # -----------------------------
  # Standard error helper
  # -----------------------------
  se_fun <- function(x) {
    x <- x[is.finite(x)]
    m <- length(x)
    if (m == 0) return(NA_real_)
    if (m == 1) return(NA_real_)
    stats::sd(x) / sqrt(m)
  }
  
  # -----------------------------
  # Compute row-wise Mean and SE
  # -----------------------------
  summary_list <- vector("list", n_rows_each)
  
  for (r in seq_len(n_rows_each)) {
    mean_row <- data.frame(matrix(nrow = 1, ncol = 0))
    se_row   <- data.frame(matrix(nrow = 1, ncol = 0))
    
    for (nm in summary_cols) {
      vals <- vapply(data_list, function(d) d[r, nm], numeric(1))
      mean_row[[nm]] <- mean(vals, na.rm = TRUE)
      se_row[[nm]]   <- se_fun(vals)
    }
    
    mean_row$SummaryType <- "Mean"
    se_row$SummaryType   <- "SE"
    
    tmp <- rbind(mean_row, se_row)
    tmp <- tmp[, c("SummaryType", summary_cols), drop = FALSE]
    
    summary_list[[r]] <- tmp
  }
  
  # -----------------------------
  # Write output files
  # -----------------------------
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  output_paths <- character(length(summary_list))
  
  if (n_rows_each == 1) {
    out_name <- sprintf(
      "Output_Summary_ALL_scen_%d_n_%s_p_%s_%s.csv",
      scenario, n, p, method_name
    )
    out_path <- file.path(output_dir, out_name)
    write.csv(summary_list[[1]], out_path, row.names = FALSE)
    output_paths[1] <- out_path
  } else {
    for (r in seq_len(n_rows_each)) {
      out_name <- sprintf(
        "Output_Summary_ALL_scen_%d_n_%s_p_%s_%s%s.csv",
        scenario, n, p, method_name, r
      )
      out_path <- file.path(output_dir, out_name)
      write.csv(summary_list[[r]], out_path, row.names = FALSE)
      output_paths[r] <- out_path
    }
  }
  
  # -----------------------------
  # Return summary object
  # -----------------------------
  invisible(list(
    files_used = basename(file_list),
    reps_used = rep_found,
    n_rows_per_file = n_rows_each,
    summary_tables = summary_list,
    output_files = output_paths
  ))
}