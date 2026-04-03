rm(list=ls())
setwd("U:/BUGS/Simulation study 1 Highdim")
source("supp funs/make_combined_summary_table.R")

scenario_here <- 1
n_here <- 200
p_here <- 10000

methods <- c(
  "LASSO",
  "UniLasso",
  "BayesRegALL1",
  "SSLASSO",
  "BayesRegALL2",
  "BayesRegALL3",
  "BUGS",
  "BUGSactive"
)


labels <- c(
  "LASSO",
  "UniLASSO",
  "BayesLASSO",
  "SSLASSO",
  "Horseshoe",
  "Horseshoe+",
  "BUGS",
  "BUGS-Active"
)





tab <- make_combined_summary_table(
  scenario = scenario_here,
  method_names = methods,
  n = n_here,
  p = p_here,
  input_dir = "Results",
  output_dir = "Results",
  method_labels = labels
)

print(tab)