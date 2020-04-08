#' ---
#' title: "Nonsmooth simulations"
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
  B = c(0, 1000),
  n = c(250, 500, 1000),
  n_i = c(5, 10, 25),
  k = 1,
  delta = c(0.1, 0.25, 0.5),
  run = 1:1000
) %>%
  filter(B == 0 | n == 250)

#' We used a server to run these simulations. If running these locally, you may want to reduce the number of runs.

sim_res <- Q_rows(sim_parameters, 
                  fun = longsurr:::dc_sim,
                  const = list(tmpdir = tmpdir),
                  n_jobs = 500,
                  fail_on_error = FALSE)
saveRDS(sim_res, 
        here('results/nonsmooth_sims.rds'))
fs::dir_delete(here('results/tmp'))
