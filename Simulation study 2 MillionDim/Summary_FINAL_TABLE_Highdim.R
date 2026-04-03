rm(list=ls())
setwd("U:/BUGS/Simulation study 2 MillionDim")
source("supp funs/make_combined_summary_table_Highdim.R")

scenario_here <- 2
n_here <- 500
p_here <- 1000000

methods <- c(
  "BUGSactive"
)


labels <- c(
  "BUGS-Active"
)





tab <- make_combined_summary_table_Highdim(
  scenario = scenario_here,
  method_names = methods,
  n = n_here,
  p = p_here,
  input_dir = "Results",
  output_dir = "Results",
  method_labels = labels
)

print(tab)