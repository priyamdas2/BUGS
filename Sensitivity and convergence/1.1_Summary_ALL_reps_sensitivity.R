rm(list=ls())
setwd("U:/BUGS/Sensitivity and convergence")

source("supp funs/summarize_metrics_across_reps_sensitivity.R")

# -----------------------------
# Settings
# -----------------------------

n_here <- 200
p_here <- 1000




for (sensitivity_casenum_here in 1:12) {   # <-- change this (1 to 12)



# -----------------------------
# Run summary
# -----------------------------
res <- summarize_metrics_across_reps_sensitivity(
  sensitivity_casenum = sensitivity_casenum_here,
  n = n_here,
  p = p_here,
  method_name = "BUGS",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

}