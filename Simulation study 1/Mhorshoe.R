# ============================================================
# Simulation study with Mhorseshoe
# Reads data from Data/, saves metrics to Results/
# ============================================================
setwd("U:/BUGS/Simulation study 1")
source("supp funs/summarize_metrics_across_reps.R")
library(Mhorseshoe)
reps <- 1:10
scenario <- 1

n <- 100
p <- 200

burn <- 1000
iter <- 5000

# -----------------------------
# Toggle which methods to run
# -----------------------------
run_exact  <- 1   # 0 or 1
run_approx <- 1   # 0 or 1

if (run_exact == 0 && run_approx == 0) {
  stop("Set at least one of run_exact or run_approx to 1.")
}

# Selection rule:
tol <- 0.01   

data_dir    <- "Data"
results_dir <- "Results"

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

for(rep_id in reps) {
  
  set.seed(rep_id)
  cat(sprintf("performing data rep: %d | scenario: %d | n: %d | p: %d\n",
              rep_id, scenario, n, p))
  
  
  
  # -----------------------------
  # Read data
  # -----------------------------
  fileX <- sprintf("%s/Data_Scenario_%d_X_n_%d_p_%d_rep_%d.csv",
                   data_dir, scenario, n, p, rep_id)
  fileY <- sprintf("%s/Data_Scenario_%d_Y_n_%d_p_%d_rep_%d.csv",
                   data_dir, scenario, n, p, rep_id)
  fileB <- sprintf("%s/Data_Scenario_%d_beta_n_%d_p_%d.csv",
                   data_dir, scenario, n, p)
  
  X <- as.matrix(read.csv(fileX, header = FALSE))
  y_raw <- as.numeric(read.csv(fileY, header = FALSE)[, 1])
  beta_true <- as.numeric(read.csv(fileB, header = FALSE)[, 1])
  
  # -----------------------------
  # Standardize y and beta_true
  # -----------------------------
  y_raw <- y_raw - mean(y_raw)
  sy <- sd(y_raw)
  
  if (sy == 0) stop("Standard deviation of y_raw is zero.")
  
  y <- y_raw / sy
  beta_true_std <- beta_true / sy
  
  # -----------------------------
  # Helper: evaluate fit
  # -----------------------------
  evaluate_fit <- function(beta_est, beta_ref, X, comp_time = NA_real_,
                           use_topk = FALSE, tol = 0.01) {
    beta_est <- as.numeric(beta_est)
    beta_ref <- as.numeric(beta_ref)
    
    # Estimation metrics
    l2_beta   <- sqrt(sum((beta_est - beta_ref)^2))
    l1_beta   <- sum(abs(beta_est - beta_ref))
    rmse_beta <- sqrt(mean((beta_est - beta_ref)^2))
    
    # Prediction metric
    pred_mse <- mean((as.vector(X %*% beta_ref) - as.vector(X %*% beta_est))^2)
    
    # Support recovery
    true_supp <- beta_ref != 0
    p <- length(beta_ref)
    
    if (use_topk) {
      k <- sum(true_supp)
      est_supp <- rep(FALSE, p)
      ord <- order(abs(beta_est), decreasing = TRUE)
      est_supp[ord[1:k]] <- TRUE
    } else {
      est_supp <- abs(beta_est) > tol
    }
    
    TP <- sum(true_supp & est_supp)
    FP <- sum(!true_supp & est_supp)
    TN <- sum(!true_supp & !est_supp)
    FN <- sum(true_supp & !est_supp)
    
    TPR <- TP / max(TP + FN, 1)
    FPR <- FP / max(FP + TN, 1)
    FDR <- FP / max(TP + FP, 1)
    
    # Avoid integer overflow in MCC
    mcc_num <- as.numeric(TP) * as.numeric(TN) - as.numeric(FP) * as.numeric(FN)
    
    a <- as.numeric(TP + FP)
    b <- as.numeric(TP + FN)
    c <- as.numeric(TN + FP)
    d <- as.numeric(TN + FN)
    
    mcc_den <- sqrt(a * b * c * d)
    
    if (is.na(mcc_den) || mcc_den <= 0) {
      MCC <- NA_real_
    } else {
      MCC <- mcc_num / mcc_den
    }
    
    data.frame(
      L2       = l2_beta,
      L1       = l1_beta,
      RMSE     = rmse_beta,
      PredMSE  = pred_mse,
      TPR      = TPR,
      FPR      = FPR,
      MCC      = MCC,
      FDR      = FDR,
      CompTime = comp_time,
      TP       = TP,
      FP       = FP,
      TN       = TN,
      FN       = FN,
      sel_thresh = ifelse(use_topk, NA, tol)
    )
  }
  # -----------------------------
  # Store outputs
  # -----------------------------
  summary_list <- list()
  
  # ============================================================
  # Run approximate horseshoe
  # ============================================================
  if (run_approx == 1) {
    t_approx <- system.time({
      fit_approx <- approx_horseshoe(
        y = y,
        X = X,
        burn = burn,
        iter = iter
      )
    })
    
    beta_est_approx <- as.numeric(fit_approx$BetaHat)
    
    metrics_approx <- evaluate_fit(
      beta_est = beta_est_approx,
      beta_ref = beta_true_std,
      X = X,
      comp_time = as.numeric(t_approx["elapsed"]),
      tol = tol
    )
    
    summary_list[[length(summary_list) + 1]] <- metrics_approx
    
    outfile_approx <- sprintf(
      "%s/Output_scen_%d_n_%d_p_%d_MhorseshoeApprox_rep_%d.csv",
      results_dir, scenario, n, p, rep_id
    )
    write.csv(metrics_approx, outfile_approx, row.names = FALSE)
  }
  
  # ============================================================
  # Run exact horseshoe
  # ============================================================
  if (run_exact == 1) {
    t_exact <- system.time({
      fit_exact <- exact_horseshoe(
        y = y,
        X = X,
        burn = burn,
        iter = iter
      )
    })
    
    beta_est_exact <- as.numeric(fit_exact$BetaHat)
    
    metrics_exact <- evaluate_fit(
      beta_est = beta_est_exact,
      beta_ref = beta_true_std,
      X = X,
      comp_time = as.numeric(t_exact["elapsed"]),
      tol = tol
    )
    
    
    summary_list[[length(summary_list) + 1]] <- metrics_exact
    
    outfile_exact <- sprintf(
      "%s/Output_scen_%d_n_%d_p_%d_MhorseshoeExact_rep_%d.csv",
      results_dir, scenario, n, p, rep_id
    )
    write.csv(metrics_exact, outfile_exact, row.names = FALSE)
  }
  
  # ============================================================
  # Save combined summary for this rep
  # ============================================================
  if (length(summary_list) > 0) {
    summary_tab <- do.call(rbind, summary_list)
    
    outfile_summary <- sprintf(
      "%s/Output_scen_%d_n_%d_p_%d_MhorseshoeALL_rep_%d.csv",
      results_dir, scenario, n, p, rep_id
    )
    write.csv(summary_tab, outfile_summary, row.names = FALSE)
    
    print(summary_tab)
  }
  
}

res <- summarize_metrics_across_reps(
  scenario = scenario,
  n = n,
  p = p,
  method_name = "MhorseshoeExact",
  reps_to_summarize = reps,
  input_dir = "Results",
  output_dir = "Results"
)

res <- summarize_metrics_across_reps(
  scenario = scenario,
  n = n,
  p = p,
  method_name = "MhorseshoeApprox",
  reps_to_summarize = reps,
  input_dir = "Results",
  output_dir = "Results"
)