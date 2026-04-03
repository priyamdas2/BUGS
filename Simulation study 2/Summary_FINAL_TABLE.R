rm(list=ls())
setwd("U:/BUGS/Simulation study 2")
source("supp funs/make_combined_summary_table.R")

scenario_here <- 2
n_here <- 150
p_here <- 500

methods <- c(
  "LASSO",
  "UniLasso",
  "BayesRegALL1",
  "R2D2DL",
  "R2D2",
  "SSLASSO",
  "BayesRegALL2",
  "BayesRegALL3",
  "MhorseshoeExact",
  "BUGS",
  "BUGSactive"
)


labels <- c(
  "LASSO",
  "UniLASSO",
  "BayesLASSO",
  "Dirich-Laplace",
  "R2D2",
  "SSLASSO",
  "Horseshoe",
  "Horseshoe+",
  "Mhorseshoe",
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