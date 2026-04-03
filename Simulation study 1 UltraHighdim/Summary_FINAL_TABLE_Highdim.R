rm(list=ls())
setwd("U:/BUGS/Simulation study 1 UltraHighdim")
source("supp funs/make_combined_summary_table_Highdim.R")

scenario_here <- 1
n_here <- 200
p_here <- 100000

methods <- c(
  "LASSO",
  "UniLasso",
  "BUGSactive"
)


labels <- c(
  "LASSO",
  "UniLASSO",
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