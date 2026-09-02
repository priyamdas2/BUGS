setwd("U:/BUGS/Simulation study 1 Highdim")
source("supp funs/summarize_metrics_across_reps.R")

# ============================================================
# Simulation study with SIS-LASSO and SIS-SCAD
# Reads data from Data/, saves metrics to Results/
# ============================================================

library(glmnet)
library(ncvreg)
library(matrixStats)

# -----------------------------
# User settings
# -----------------------------
scenario <- 1

n <- 200
p <- 10000

num_reps <- 10

# Selection rule:
tol <- 0.01   # only used if use_topk = FALSE

data_dir <- "Data"
results_dir <- "Results"

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# -----------------------------
# Helper: evaluate fit
# -----------------------------
evaluate_fit <- function(beta_est, beta_ref, X, comp_time = NA_real_,
                         use_topk = FALSE, tol = 0.01,
                         intercept_est = 0, intercept_ref = 0) {
  
  beta_est <- as.numeric(beta_est)
  beta_ref <- as.numeric(beta_ref)
  
  # Estimation metrics
  l2_beta   <- sqrt(sum((beta_est - beta_ref)^2))
  l1_beta   <- sum(abs(beta_est - beta_ref))
  rmse_beta <- sqrt(mean((beta_est - beta_ref)^2))
  
  # Prediction metric
  pred_ref <- as.vector(intercept_ref + X %*% beta_ref)
  pred_est <- as.vector(intercept_est + X %*% beta_est)
  pred_mse <- mean((pred_ref - pred_est)^2)
  
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
  mcc_num <- as.numeric(TP) * as.numeric(TN) -
    as.numeric(FP) * as.numeric(FN)
  
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
    L2         = l2_beta,
    L1         = l1_beta,
    RMSE       = rmse_beta,
    PredMSE    = pred_mse,
    TPR        = TPR,
    FPR        = FPR,
    MCC        = MCC,
    FDR        = FDR,
    CompTime   = comp_time,
    TP         = TP,
    FP         = FP,
    TN         = TN,
    FN         = FN,
    sel_thresh = ifelse(use_topk, NA, tol)
  )
}

# ============================================================
# Loop over reps
# ============================================================

for (rep_id in 1:num_reps) {
  
  cat(sprintf(
    "\nperforming data rep: %d | scenario: %d | n: %d | p: %d\n",
    rep_id, scenario, n, p
  ))
  
  set.seed(rep_id)
  
  # -----------------------------
  # Read data
  # -----------------------------
  fileX <- sprintf(
    "%s/Data_Scenario_%d_X_n_%d_p_%d_rep_%d.csv",
    data_dir, scenario, n, p, rep_id
  )
  
  fileY <- sprintf(
    "%s/Data_Scenario_%d_Y_n_%d_p_%d_rep_%d.csv",
    data_dir, scenario, n, p, rep_id
  )
  
  fileB <- sprintf(
    "%s/Data_Scenario_%d_beta_n_%d_p_%d.csv",
    data_dir, scenario, n, p
  )
  
  X <- as.matrix(read.csv(fileX, header = FALSE))
  y_raw <- as.numeric(read.csv(fileY, header = FALSE)[, 1])
  beta_true <- as.numeric(read.csv(fileB, header = FALSE)[, 1])
  
  # -----------------------------
  # Standardize y and beta_true
  # -----------------------------
  y_raw <- y_raw - mean(y_raw)
  sy <- sd(y_raw)
  
  if (sy == 0) {
    stop("Standard deviation of y_raw is zero.")
  }
  
  y <- y_raw / sy
  beta_true_std <- beta_true / sy
  
  
  # ============================================================
  # Step 1: Sure Independence Screening
  # ============================================================
  
  t_sis <- system.time({
    
    # Fan and Lv (2008):
    # retain d = floor(n / log(n))
    d_sis <- min(p, floor(n / log(n)))
    
    # Fast predictor standard deviations
    x_sd <- matrixStats::colSds(X)
    
    # Absolute marginal correlation scores
    #
    # Since y is centered, centering each X_j does not affect
    # X_j' y. Scaling by sd(X_j) gives the SIS ranking.
    sis_score <- abs(as.numeric(crossprod(X, y))) / x_sd
    
    # Handle zero-variance or otherwise invalid predictors
    sis_score[!is.finite(sis_score)] <- -Inf
    
    # Select predictors with largest SIS scores
    sis_idx <- order(
      sis_score,
      decreasing = TRUE
    )[1:d_sis]
    
    # Reduced design matrix
    X_sis <- X[, sis_idx, drop = FALSE]
  })
  
  sis_time <- as.numeric(t_sis["elapsed"])
  
  cat(sprintf(
    "  SIS retained %d of %d predictors | time: %.3f sec\n",
    d_sis, p, sis_time
  ))
  
  
  # ============================================================
  # Step 2A: SIS-LASSO
  # ============================================================
  
  t_lasso <- system.time({
    
    fit_sis_lasso <- cv.glmnet(
      x = X_sis,
      y = y,
      alpha = 1
    )
  })
  
  lasso_time <- as.numeric(t_lasso["elapsed"])
  
  # Extract coefficients at lambda.min
  coef_lasso <- as.matrix(
    coef(fit_sis_lasso, s = "lambda.min")
  )
  
  intercept_sis_lasso <- as.numeric(coef_lasso[1, 1])
  beta_est_lasso_reduced <- as.numeric(coef_lasso[-1, 1])
  
  # Expand back to original p-dimensional space
  beta_est_sis_lasso <- rep(0, p)
  beta_est_sis_lasso[sis_idx] <- beta_est_lasso_reduced
  
  # Total computational time = SIS + LASSO
  total_time_sis_lasso <- sis_time + lasso_time
  
  metrics_sis_lasso <- evaluate_fit(
    beta_est = beta_est_sis_lasso,
    beta_ref = beta_true_std,
    X = X,
    comp_time = total_time_sis_lasso,
    tol = tol,
    intercept_est = intercept_sis_lasso,
    intercept_ref = 0
  )
  
  # Save SIS-LASSO results separately
  outfile_sis_lasso <- sprintf(
    "%s/Output_scen_%d_n_%d_p_%d_SIS_LASSO_rep_%d.csv",
    results_dir, scenario, n, p, rep_id
  )
  
  write.csv(
    metrics_sis_lasso,
    outfile_sis_lasso,
    row.names = FALSE
  )
  
  cat(sprintf(
    "  SIS-LASSO fit time: %.3f sec | total: %.3f sec\n",
    lasso_time, total_time_sis_lasso
  ))
  
  print(metrics_sis_lasso)
  
  
  # ============================================================
  # Step 2B: SIS-SCAD
  # ============================================================
  
  t_scad <- system.time({
    
    fit_sis_scad <- cv.ncvreg(
      X = X_sis,
      y = y,
      family = "gaussian",
      penalty = "SCAD",
      gamma = 3.7
    )
  })
  
  scad_time <- as.numeric(t_scad["elapsed"])
  
  # coef(cv.ncvreg object) returns coefficients corresponding
  # to lambda.min
  coef_scad <- as.numeric(coef(fit_sis_scad))
  
  intercept_sis_scad <- coef_scad[1]
  beta_est_scad_reduced <- coef_scad[-1]
  
  # Expand back to original p-dimensional space
  beta_est_sis_scad <- rep(0, p)
  beta_est_sis_scad[sis_idx] <- beta_est_scad_reduced
  
  # Total computational time = SIS + SCAD
  total_time_sis_scad <- sis_time + scad_time
  
  metrics_sis_scad <- evaluate_fit(
    beta_est = beta_est_sis_scad,
    beta_ref = beta_true_std,
    X = X,
    comp_time = total_time_sis_scad,
    tol = tol,
    intercept_est = intercept_sis_scad,
    intercept_ref = 0
  )
  
  # Save SIS-SCAD results separately
  outfile_sis_scad <- sprintf(
    "%s/Output_scen_%d_n_%d_p_%d_SIS_SCAD_rep_%d.csv",
    results_dir, scenario, n, p, rep_id
  )
  
  write.csv(
    metrics_sis_scad,
    outfile_sis_scad,
    row.names = FALSE
  )
  
  cat(sprintf(
    "  SIS-SCAD fit time: %.3f sec | total: %.3f sec\n",
    scad_time, total_time_sis_scad
  ))
  
  print(metrics_sis_scad)
}


# ============================================================
# Summarize SIS-LASSO across reps
# ============================================================

res_sis_lasso <- summarize_metrics_across_reps(
  scenario = scenario,
  n = n,
  p = p,
  method_name = "SIS_LASSO",
  reps_to_summarize = 1:num_reps,
  input_dir = "Results",
  output_dir = "Results"
)


# ============================================================
# Summarize SIS-SCAD across reps
# ============================================================

res_sis_scad <- summarize_metrics_across_reps(
  scenario = scenario,
  n = n,
  p = p,
  method_name = "SIS_SCAD",
  reps_to_summarize = 1:num_reps,
  input_dir = "Results",
  output_dir = "Results"
)