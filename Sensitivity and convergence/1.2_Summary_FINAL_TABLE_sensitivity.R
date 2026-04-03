rm(list=ls())
setwd("U:/BUGS/Sensitivity and convergence")
source("supp funs/make_combined_summary_table_sensitivity.R")

n_here <- 200
p_here <- 1000

tab_sens <- make_combined_summary_table_sensitivity(
  sensitivity_casenums = 1:12,
  method_names = "BUGS",
  n = n_here,
  p = p_here,
  input_dir = "Results",
  output_dir = "Results"
)

tab_sens