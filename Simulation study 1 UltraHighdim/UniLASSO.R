setwd("U:/BUGS/Simulation study 1 UltraHighdim")
source("supp funs/summarize_metrics_across_reps_Highdim.R")
library(uniLasso)

# ============================================================
# Simulation study with uniLasso
# Reads data from Data/, saves metrics to Results/
# ============================================================

# -----------------------------
# User settings
# -----------------------------
scenario <- 1

n <- 200
p <- 100000

reps <- 1:10

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

# ============================================================
# Loop over reps
# ============================================================
for (rep_id in reps) {
  
  cat(sprintf("performing data rep: %d | scenario: %d | n: %d | p: %d\n",
              rep_id, scenario, n, p))
  
  set.seed(rep_id)
  
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
  
  # ============================================================
  # Run uniLasso (package defaults)
  # ============================================================
  t_unilasso <- system.time({
    fit_unilasso <- cv.uniLasso(
      x = X,
      y = y
    )
  })
  
  # Extract coefficients at lambda.min
  coef_mat <- as.matrix(coef(fit_unilasso, s = "lambda.min"))
  intercept_unilasso <- as.numeric(coef_mat[1, 1])
  beta_est_unilasso <- as.numeric(coef_mat[-1, 1])
  
  metrics_unilasso <- evaluate_fit(
    beta_est = beta_est_unilasso,
    beta_ref = beta_true_std,
    X = X,
    comp_time = as.numeric(t_unilasso["elapsed"]),
    tol = tol,
    intercept_est = intercept_unilasso,
    intercept_ref = 0
  )
  
  outfile_unilasso <- sprintf(
    "%s/Output_scen_%d_n_%d_p_%d_UniLasso_rep_%d.csv",
    results_dir, scenario, n, p, rep_id
  )
  write.csv(metrics_unilasso, outfile_unilasso, row.names = FALSE)
  
  print(metrics_unilasso)
}

# ============================================================
# Summarize across reps
# ============================================================
res <- summarize_metrics_across_reps_Highdim(
  scenario = scenario,
  n = n,
  p = p,
  method_name = "UniLasso",
  reps_to_summarize = reps,
  input_dir = "Results",
  output_dir = "Results"
)