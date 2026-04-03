make_combined_summary_table_sensitivity <- function(
    sensitivity_casenums = 1:12,
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
    method_labels = NULL,
    case_labels = NULL
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
      CompTime  = 3
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
      CompTime  = 4
    )
  }
  
  if (is.null(method_labels)) {
    method_labels <- method_names
  }
  
  if (length(method_labels) != length(method_names)) {
    stop("method_labels and method_names must have the same length.")
  }
  
  if (is.null(case_labels)) {
    case_labels <- paste("Case", sensitivity_casenums)
  }
  
  if (length(case_labels) != length(sensitivity_casenums)) {
    stop("case_labels and sensitivity_casenums must have the same length.")
  }
  
  fmt_num <- function(x, digits) {
    if (is.na(x) || !is.finite(x)) return(NA_character_)
    formatC(x, format = "f", digits = digits)
  }
  
  fmt_mean_se <- function(mean_val, se_val, mean_d, se_d) {
    mean_txt <- fmt_num(mean_val, mean_d)
    se_txt   <- fmt_num(se_val, se_d)
    
    if (is.na(mean_txt) && is.na(se_txt)) return(NA_character_)
    if (is.na(mean_txt)) return(NA_character_)
    if (is.na(se_txt)) return(mean_txt)
    
    paste0(mean_txt, " (", se_txt, ")")
  }
  
  standardize_metric_names <- function(dat) {
    nm <- names(dat)
    
    preferred_cols <- c(
      "sensitivity_casenum", "Method", "SummaryType",
      "L2_beta", "L1_beta", "RMSE_beta", "PredMSE",
      "TPR", "FPR", "MCC", "FDR", "CompTime"
    )
    
    # If columns already look fine, leave as is
    if (all(intersect(c("SummaryType", "L2_beta", "CompTime"), preferred_cols) %in% c(nm, preferred_cols))) {
      return(dat)
    }
    
    first8_std <- c("L2_beta", "L1_beta", "RMSE_beta", "PredMSE",
                    "TPR", "FPR", "MCC", "FDR")
    
    if ("SummaryType" %in% nm) {
      idx_start <- match("SummaryType", nm) + 1
      idx_end <- min(idx_start + 7, length(nm))
      if ((idx_end - idx_start + 1) == 8) {
        nm[idx_start:idx_end] <- first8_std
      }
      if (length(nm) >= idx_start + 8) {
        nm[idx_start + 8] <- "CompTime"
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
  
  build_file <- function(casenum, method) {
    file.path(
      input_dir,
      sprintf("Output_Summary_ALL_casenum_%d_n_%s_p_%s_%s.csv",
              casenum, n, p, method)
    )
  }
  
  out_list <- list()
  idx_out <- 1
  
  for (j in seq_along(sensitivity_casenums)) {
    casenum <- sensitivity_casenums[j]
    case_label <- case_labels[j]
    
    for (i in seq_along(method_names)) {
      method <- method_names[i]
      method_label <- method_labels[i]
      file_path <- build_file(casenum, method)
      
      if (!file.exists(file_path)) {
        warning(sprintf("Missing file for %s, %s: %s",
                        case_label, method_label, basename(file_path)))
        next
      }
      
      dat <- tryCatch(
        read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      
      if (is.null(dat)) {
        warning(sprintf("Could not read file: %s", basename(file_path)))
        next
      }
      
      dat <- standardize_metric_names(dat)
      
      if (nrow(dat) < 1) {
        warning(sprintf("Empty file: %s", basename(file_path)))
        next
      }
      
      if ("SummaryType" %in% names(dat)) {
        mean_row <- dat[dat$SummaryType == "Mean", , drop = FALSE]
        se_row   <- dat[dat$SummaryType == "SE", , drop = FALSE]
        
        if (nrow(mean_row) == 0 && nrow(dat) >= 1) {
          warning(sprintf("Mean row not found in %s; using first row.", basename(file_path)))
          mean_row <- dat[1, , drop = FALSE]
        }
        
        if (nrow(se_row) == 0) {
          if (nrow(dat) >= 2) {
            warning(sprintf("SE row not found in %s; using second row.", basename(file_path)))
            se_row <- dat[2, , drop = FALSE]
          } else {
            warning(sprintf("SE row not found in %s; SE will be missing.", basename(file_path)))
            se_row <- mean_row
            for (metric in metric_cols) {
              if (metric %in% names(se_row)) {
                se_row[[metric]] <- NA_real_
              }
            }
          }
        }
      } else {
        if (nrow(dat) >= 2) {
          mean_row <- dat[1, , drop = FALSE]
          se_row   <- dat[2, , drop = FALSE]
        } else {
          warning(sprintf("Only one row found in %s; SE will be missing.", basename(file_path)))
          mean_row <- dat[1, , drop = FALSE]
          se_row   <- mean_row
          for (metric in metric_cols) {
            if (metric %in% names(se_row)) {
              se_row[[metric]] <- NA_real_
            }
          }
        }
      }
      
      one_row <- data.frame(
        Case = case_label,
        sensitivity_casenum = casenum,
        Method = method_label,
        stringsAsFactors = FALSE
      )
      
      for (metric in metric_cols) {
        if (!(metric %in% names(mean_row)) || !(metric %in% names(se_row))) {
          warning(sprintf("Metric %s not found in file %s", metric, basename(file_path)))
          one_row[[metric]] <- NA_character_
          next
        }
        
        mean_val <- suppressWarnings(as.numeric(mean_row[[metric]][1]))
        se_val   <- suppressWarnings(as.numeric(se_row[[metric]][1]))
        
        md <- if (metric %in% names(mean_digits)) mean_digits[[metric]] else 3
        sd <- if (metric %in% names(se_digits)) se_digits[[metric]] else 4
        
        one_row[[metric]] <- fmt_mean_se(mean_val, se_val, md, sd)
      }
      
      out_list[[idx_out]] <- one_row
      idx_out <- idx_out + 1
    }
  }
  
  out_list <- Filter(Negate(is.null), out_list)
  
  if (length(out_list) == 0) {
    warning("No valid summary files were found. Returning empty data frame.")
    final_tab <- data.frame(
      Case = character(0),
      sensitivity_casenum = integer(0),
      Method = character(0),
      stringsAsFactors = FALSE
    )
    
    if (is.null(output_file)) {
      output_file <- sprintf("ALL_Combined_Summary_casenums_n_%s_p_%s.csv", n, p)
    }
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    write.csv(final_tab, file.path(output_dir, output_file), row.names = FALSE)
    return(final_tab)
  }
  
  final_tab <- do.call(rbind, out_list)
  
  # Order by case then method
  final_tab <- final_tab[order(final_tab$sensitivity_casenum, final_tab$Method), , drop = FALSE]
  rownames(final_tab) <- NULL
  
  if (is.null(output_file)) {
    if (length(sensitivity_casenums) == 1) {
      output_file <- sprintf(
        "ALL_Combined_Summary_casenum_%d_n_%s_p_%s.csv",
        sensitivity_casenums, n, p
      )
    } else {
      output_file <- sprintf(
        "ALL_Combined_Summary_casenums_%d_to_%d_n_%s_p_%s.csv",
        min(sensitivity_casenums), max(sensitivity_casenums), n, p
      )
    }
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  write.csv(final_tab, file.path(output_dir, output_file), row.names = FALSE)
  
  return(final_tab)
}