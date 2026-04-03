rm(list=ls())
setwd("U:/BUGS/Simulation study 2 UltraHighdim")
source("supp funs/summarize_metrics_across_reps_Highdim.R")

scenario_here <- 2
n_here <- 200
p_here <- 100000


res <- summarize_metrics_across_reps_Highdim(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "BUGSactive",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)


res <- summarize_metrics_across_reps_Highdim(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "LASSO",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

res <- summarize_metrics_across_reps_Highdim(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "UniLasso",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

