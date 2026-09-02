rm(list=ls())
setwd("U:/BUGS/Simulation study 2 MillionDim")
source("supp funs/summarize_metrics_across_reps_Highdim.R")

scenario_here <- 2
n_here <- 500
p_here <- 1000000

res <- summarize_metrics_across_reps_Highdim(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "SIS_LASSO",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

res <- summarize_metrics_across_reps_Highdim(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "SIS_SCAD",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

res <- summarize_metrics_across_reps_Highdim(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "SIS_RHS",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)


res <- summarize_metrics_across_reps_Highdim(
  scenario = scenario_here,
  n = n_here,
  p = p_here,
  method_name = "BUGSactive",
  reps_to_summarize = 1:10,
  input_dir = "Results",
  output_dir = "Results"
)

