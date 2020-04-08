#' ---
#' title: "Repeatedly subsample to find best time grid"
#' output: github_document
#' ---
#' 
library(knitr)
opts_chunk$set(warning = FALSE, message = FALSE, cache = FALSE, fig.width = 7, fig.height = 7)

#'
#'
library(tidyverse)
library(here)
library(zeallot)
library(longsurr)
library(refund)
library(fda.usc)
library(Rsurrogate)
library(clustermq)

select <- dplyr::select
analysis_data <- read_csv(here('data/hiv-analysis-data.csv')) %>%
  arrange(id, tt)
all_ids <- analysis_data %>%
  select(id) %>% 
  unique

c(trt_xhat_wide, ctrl_xhat_wide) %<-%
  presmooth_data(obs_data = analysis_data, 
                 options = 
                   list(plot = TRUE, 
                        # methodBwCov = 'GCV',
                        methodBwMu = 'CV',
                        methodSelectK = 'AIC',
                        useBinnedCov = FALSE,
                        verbose = TRUE,
                        usergrid = FALSE,
                        nRegGrid = 51))

smoothed_data <- as_tibble(trt_xhat_wide, rownames = 'id') %>%
  mutate(a = 1) %>%
  full_join(
    as_tibble(ctrl_xhat_wide, rownames = 'id') %>%
      mutate(a = 0)
  ) %>%
  gather(tt, x, -id, -a) %>%
  mutate(tt = as.numeric(tt))

ttl <- read_rds(here('data/time-grids-for-cv.rds'))

fs::dir_create(here('results/tmp'))
sim_pars <- tidyr::expand_grid(s = 1:500,
                               i = 1:length(ttl))

#' We ran the following on a high-performance cluster. If running locally, you may want to select fewer time grids and/or fewer subsamples.
sim_res <- Q_rows(sim_pars,
                  longsurr:::hiv_cv,
                  n_jobs = 500,
                  const = list(
                    time_list = ttl,
                    all_ids = all_ids,
                    analysis_data = analysis_data,
                    smoothed_data = smoothed_data,
                    trt_xhat_wide = trt_xhat_wide,
                    tmpdir = here('results/tmp')
                  ),
                  pkgs=c('tidyverse',
                         'here',
                         'zeallot',
                         'longsurr',
                         'refund',
                         'fda.usc',
                         'Rsurrogate'),
                  fail_on_error = FALSE)
write_rds(sim_res, here('results/hiv-full-cv-results.rds'))  
fs::dir_delete(here('results/tmp'))
