make_combined_summary_table <- function(
    scenario,
    method_names,
    n,
    p,
    input_dir = ".",
    output_dir = input_dir,
    output_file = NULL,
    metric_cols = c("L2_beta", "L1_beta", "RMSE_beta", "PredMSE",
                    "TPR", "FPR", "MCC", "FDR", "CompTime"),
    mean_digits = NULL,
    se_digits = NULL,
    method_labels = NULL
) {
  if (is.null(mean_digits)) {
    mean_digits <- c(
      L2_beta   = 3,
      L1_beta   = 3,
      RMSE_beta = 3,
      PredMSE   = 3,
      TPR       = 3,
      FPR       = 3,
      MCC       = 3,
      FDR       = 3,
      CompTime  = 2
    )
  }
  
  if (is.null(se_digits)) {
    se_digits <- c(
      L2_beta   = 4,
      L1_beta   = 4,
      RMSE_beta = 4,
      PredMSE   = 4,
      TPR       = 4,
      FPR       = 4,
      MCC       = 4,
      FDR       = 4,
      CompTime  = 3
    )
  }
  
  if (is.null(method_labels)) {
    method_labels <- method_names
  }
  
  if (length(method_labels) != length(method_names)) {
    stop("method_labels and method_names must have the same length.")
  }
  
  fmt_num <- function(x, digits) {
    if (is.na(x)) return(NA_character_)
    formatC(x, format = "f", digits = digits)
  }
  
  fmt_mean_se <- function(mean_val, se_val, mean_d, se_d) {
    paste0(fmt_num(mean_val, mean_d), " (", fmt_num(se_val, se_d), ")")
  }
  
  standardize_metric_names <- function(dat) {
    nm <- names(dat)
    
    first8_std <- c("L2_beta", "L1_beta", "RMSE_beta", "PredMSE",
                    "TPR", "FPR", "MCC", "FDR")
    
    if ("SummaryType" %in% nm) {
      if (length(nm) >= 9) {
        nm[2:9] <- first8_std
      }
      if (length(nm) >= 10) {
        nm[10] <- "CompTime"
      }
    } else {
      if (length(nm) >= 8) {
        nm[1:8] <- first8_std
      }
      if (length(nm) >= 9) {
        nm[9] <- "CompTime"
      }
    }
    
    names(dat) <- nm
    dat
  }
  
  build_file <- function(method) {
    file.path(
      input_dir,
      sprintf("Output_Summary_ALL_scen_%d_n_%s_p_%s_%s.csv",
              scenario, n, p, method)
    )
  }
  
  out_list <- vector("list", length(method_names))
  
  for (i in seq_along(method_names)) {
    method <- method_names[i]
    method_label <- method_labels[i]
    file_path <- build_file(method)
    
    if (!file.exists(file_path)) {
      warning(sprintf("File not found: %s", basename(file_path)))
      next
    }
    
    dat <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
    dat <- standardize_metric_names(dat)
    
    if (nrow(dat) < 2) {
      warning(sprintf("File has fewer than 2 rows: %s", basename(file_path)))
      next
    }
    
    if ("SummaryType" %in% names(dat)) {
      mean_row <- dat[dat$SummaryType == "Mean", , drop = FALSE]
      se_row   <- dat[dat$SummaryType == "SE", , drop = FALSE]
      
      if (nrow(mean_row) == 0) mean_row <- dat[1, , drop = FALSE]
      if (nrow(se_row) == 0)   se_row   <- dat[2, , drop = FALSE]
    } else {
      mean_row <- dat[1, , drop = FALSE]
      se_row   <- dat[2, , drop = FALSE]
    }
    
    one_row <- data.frame(Method = method_label, stringsAsFactors = FALSE)
    
    for (metric in metric_cols) {
      if (!(metric %in% names(mean_row)) || !(metric %in% names(se_row))) {
        warning(sprintf("Metric %s not found in file %s", metric, basename(file_path)))
        one_row[[metric]] <- NA_character_
        next
      }
      
      mean_val <- suppressWarnings(as.numeric(mean_row[[metric]][1]))
      se_val   <- suppressWarnings(as.numeric(se_row[[metric]][1]))
      
      md <- if (!is.null(mean_digits[metric])) mean_digits[metric] else 3
      sd <- if (!is.null(se_digits[metric])) se_digits[metric] else 4
      
      one_row[[metric]] <- fmt_mean_se(mean_val, se_val, md, sd)
    }
    
    out_list[[i]] <- one_row
  }
  
  out_list <- Filter(Negate(is.null), out_list)
  
  if (length(out_list) == 0) {
    stop("No valid summary files were read.")
  }
  
  final_tab <- do.call(rbind, out_list)
  
  if (is.null(output_file)) {
    output_file <- sprintf("ALL_Combined_Summary_scen_%d_n_%s_p_%s.csv",
                           scenario, n, p)
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  write.csv(final_tab, file.path(output_dir, output_file), row.names = FALSE)
  
  return(final_tab)
}