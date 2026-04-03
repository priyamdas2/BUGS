rm(list=ls())
setwd("U:/BUGS/Simulation study 1")
source("supp funs/summarize_metrics_across_reps.R")

scenario_here <- 1
n_here <- 100
p_here <- 200

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
  method_name = "BUGSRegHS",
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
  method_name = "MhorseshoeExact",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

res <- summarize_metrics_across_reps(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "MhorseshoeApprox",
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
  method_name = "R2D2",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

res <- summarize_metrics_across_reps(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "R2D2DL",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)
