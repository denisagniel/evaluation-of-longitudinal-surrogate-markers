#' ---
#' title: "Using CV to choose best model and time grid"
#' output: github_document
#' ---
#' 
library(knitr)
opts_chunk$set(warning = FALSE, message = FALSE, cache = FALSE, fig.width = 7, fig.height = 7)

#'
#'
library(tidyverse)
library(glue)
library(here)


cv_res <- read_rds(here('results/hiv-full-cv-results.rds')) %>%
  bind_rows %>%
  mutate(s = sim_pars$s,
         t_i = sim_pars$i)
ttl <- read_rds(here('data/time-grids-for-cv.rds'))
sim_pars <- tidyr::expand_grid(s = 1:500,
                               i = 1:length(ttl))
tt_info <- map(ttl, function(t) {
  tibble(min_t = t[1],
         max_t = t[length(t)],
         length_t = length(t))
}) %>%
  bind_rows(.id = 't_i') %>%
  mutate(t_i = as.numeric(t_i))

cv_sum <- cv_res %>%
  inner_join(tt_info) %>%
  group_by(t_i, min_t, max_t, length_t) %>%
  summarise_all(mean) %>%
  ungroup %>%
  transmute(
    t_i, min_t, max_t, length_t,
    g_mse,
    g_pct_min = g_mse/min(g_mse),
    k3_mse,
    k3_pct_min = k3_mse/min(k3_mse),
    k4_mse,
    k4_pct_min = k4_mse/min(k4_mse),
    l_mse,
    l_pct_min = l_mse/min(l_mse)
  )
write_csv(cv_sum, here('results/cv-summarized-results.csv'))
