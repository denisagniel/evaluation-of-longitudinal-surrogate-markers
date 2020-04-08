#' ---
#' title: "ACTG-based simulations"
#' output: github_document
#' ---
#' 
library(knitr)
opts_chunk$set(warning = FALSE, message = FALSE, cache = FALSE, fig.width = 7, fig.height = 7)

#'

library(tidyverse)
library(here)
library(zeallot)
library(longsurr)
library(refund)
library(fda.usc)
library(Rsurrogate)
library(clustermq)

analysis_data <- read_csv(here('data/hiv-analysis-data.csv'))
smoothed_data <- read_rds(here('results/hiv-smoothed-data.rds'))
sim_list <- read_rds(here('results/hiv-sim-tools.rds'))
final_list <- map(1:(length(sim_list) + 2), function(i) {
  if (i %in% 1:length(sim_list)) return(sim_list[[i]])
  else if (i == length(sim_list) + 1) smoothed_data
  else analysis_data
})
names(final_list) <- c(names(sim_list), 'smoothed_data', 'analysis_data')
sim_pars <- expand.grid(mean_fn = c('kernel', 'gam', 'linear'),
                        s = 1:1000)

#' We used a server to run these simulations. If running these locally, you may want to reduce the number of runs.

sim_res <- Q_rows(sim_pars, 
                  longsurr:::hiv_sim_fn,
                  n_jobs = 100,
                  export=final_list,
                  pkgs=c('tidyverse',
                         'here',
                         'zeallot',
                         'longsurr',
                         'refund',
                         'fda.usc',
                         'Rsurrogate'))
write_rds(sim_res, here('results/hiv-full-sim-reanalysis-results.rds'))  

