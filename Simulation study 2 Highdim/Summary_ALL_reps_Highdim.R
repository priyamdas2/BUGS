rm(list=ls())
setwd("U:/BUGS/Simulation study 2 Highdim")
source("supp funs/summarize_metrics_across_reps.R")

scenario_here <- 2
n_here <- 200
p_here <- 10000

res <- summarize_metrics_across_reps(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "BUGS",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

res <- summarize_metrics_across_reps(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "BUGSactive",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

res <- summarize_metrics_across_reps(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "BayesRegALL",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)


res <- summarize_metrics_across_reps(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "SSLASSO",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

res <- summarize_metrics_across_reps(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "LASSO",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

res <- summarize_metrics_across_reps(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "UniLasso",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

res <- summarize_metrics_across_reps(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "SIS_LASSO",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

res <- summarize_metrics_across_reps(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "SIS_SCAD",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

res <- summarize_metrics_across_reps(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "SIS_RHS",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

