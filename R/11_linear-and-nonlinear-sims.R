#' ---
#' title: "Linear and nonlinear simulations"
#' output: github_document
#' ---
#' 
library(knitr)
opts_chunk$set(warning = FALSE, message = FALSE, cache = FALSE, fig.width = 7, fig.height = 7)

#'
library(dplyr)
library(clustermq)
library(here)
library(purrr)
library(readr)
library(glue)
library(ggplot2)
library(longsurr)

fs::dir_create(here('results/tmp'))

sim_parameters <- expand.grid(
  run = 1:1000,
  n = c(250, 500, 1000),
  n_i = c(5, 10, 25),
  m = c('linear', 'nonlinear'),
  s_y = 1,
  s_x = 1,
  delta = c(1, 5, 15, 25),
  B = c(0, 1000)
) %>%
  filter(B == 0 | n == 250)

#' We used a server to run these simulations. If running these locally, you may want to reduce the number of runs.

sim_res <- Q(longsurr:::lsa_sim, 
             n = sim_parameters$n,
             n_i = sim_parameters$n_i,
             m = sim_parameters$m,
             s_y = sim_parameters$s_y,
             s_x = sim_parameters$s_x,
             B = sim_parameters$B,
             delta = sim_parameters$delta,
             run = sim_parameters$run,
             const = list(tmpdir = here('results/tmp')),
             n_jobs = 500,
             # memory = 2000,
             fail_on_error = FALSE
)
saveRDS(sim_res, 
        here('results/linear_and_nonlinear_sims.rds'))
fs::dir_delete(here('results/tmp'))
