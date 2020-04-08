#' ---
#' title: "ACTG 175 analysis"
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
select <- dplyr::select
analysis_data <- read_csv(
  here('data/hiv-analysis-data.csv'))
cv_results <- read_csv(here('results/cv-summarized-results.csv'))
set.seed(0)
trt_ds <- analysis_data  %>%
  filter(a == 1)
ctrl_ds <- analysis_data %>%
  filter(a == 0)

n_trt <- trt_ds %>%
  summarise(n_trt = length(unique(id))) %>%
  pull(n_trt)
y_t <- trt_ds %>%
  select(id, y) %>%
  unique %>%
  pull(y)

n_ctrl <- ctrl_ds %>%
  summarise(n_ctrl = length(unique(id))) %>%
  pull(n_ctrl)
y_c <- ctrl_ds %>%
  select(id, y) %>%
  unique %>%
  pull(y)


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

overall_treatment_effect <- mean(y_t) - mean(y_c)

#' Get the grids for each treatment.

k_res <- cv_results %>%
  arrange(k4_mse) %>%
  head(1) 
k_grid <- seq(k_res$min_t, k_res$max_t, length = k_res$length_t) %>% round
k_xt <- trt_xhat_wide[,k_grid]
k_xc <- ctrl_xhat_wide[,k_grid]

g_res <- cv_results %>%
  arrange(g_mse) %>%
  head(1) 
g_grid <- seq(g_res$min_t, g_res$max_t, length = g_res$length_t) %>% round
g_xt <- trt_xhat_wide[,g_grid]
g_xc <- ctrl_xhat_wide[,g_grid]

l_res <- cv_results %>%
  arrange(l_mse) %>%
  head(1) 
l_grid <- seq(l_res$min_t, l_res$max_t, length = l_res$length_t) %>% round
l_xt <- trt_xhat_wide[,l_grid]
l_xc <- ctrl_xhat_wide[,l_grid]

#'
#' ### Results on smoothed data
linear_res <- estimate_surrogate_value(y_t = y_t,
                                       y_c = y_c,
                                       X_t = l_xt,
                                       X_c = l_xc,
                                       method = 'linear',
                                       bootstrap_samples = 250)
gam_res <- estimate_surrogate_value(y_t = y_t,
                                    y_c = y_c,
                                    X_t = g_xt,
                                    X_c = g_xc,
                                    method = 'gam',
                                    bootstrap_samples = 250)
k_res <- estimate_surrogate_value(y_t = y_t,
                                  y_c = y_c,
                                  X_t = k_xt,
                                  X_c = k_xc,
                                  method = 'kernel',
                                  bootstrap_samples = 250,
                                  k = 4)


delta_res <- 
  full_join(linear_res %>% mutate(method = 'linear'), 
            gam_res %>% mutate(method = 'gam')) %>%
  full_join(k_res %>% mutate(method = 'kernel')) %>%
  kable(digits = 3)


delta_res
write_csv(delta_res, here('results/hiv-results.csv'))

#' Save some stuff for simulation.
smoothed_data <- as_tibble(trt_xhat_wide, rownames = 'id') %>%
  mutate(a = 1) %>%
  full_join(
    as_tibble(ctrl_xhat_wide, rownames = 'id') %>%
      mutate(a = 0)
  ) %>%
  gather(tt, xh, -id, -a) %>%
  mutate(tt = as.numeric(tt),
         id = as.numeric(id)) %>%
  full_join(analysis_data %>%
              dplyr::select(id, tt, x, a)) %>%
  inner_join(analysis_data %>%
               dplyr::select(id, y) %>%
               unique)
write_rds(smoothed_data, here('results/hiv-smoothed-data.rds'))

trt_guys <- analysis_data %>%
  filter(a == 1) %>%
  dplyr::select(id, y) %>%
  unique
y_t <- trt_guys %>%
  pull(y)
control_guys <- analysis_data %>%
  filter(a == 0) %>%
  dplyr::select(id, y) %>%
  unique
y_c <- control_guys %>%
  pull(y)
kernel_fit_1 <- fregre.np.cv(fdataobj = fdata(k_xt), y = y_t,
                             metric = longsurr:::make_semimetric_pca(4))
kernel_fit_0 <- fregre.np.cv(fdataobj = fdata(k_xc), y = y_c,
                             metric = longsurr:::make_semimetric_pca(4))
lin_1 <- fgam(y_t ~ lf(l_xt))
lin_0 <- fgam(y_c ~ lf(l_xc))
fgam_fit_1 <- fgam(y_t ~ af(g_xt))
fgam_fit_0 <- fgam(y_c ~ af(g_xc))
k_sigma_0 <- kernel_fit_0$sr2 %>% sqrt
k_sigma_1 <- kernel_fit_1$sr2 %>% sqrt
l_sigma_0 <- lin_0$sig2 %>% sqrt
l_sigma_1 <- lin_1$sig2 %>% sqrt
g_sigma_0 <- fgam_fit_0$sig2 %>% sqrt
g_sigma_1 <- fgam_fit_1$sig2 %>% sqrt

k_trt_guys <- trt_guys %>%
  mutate(y_1 = predict(kernel_fit_1),
         y_0 = predict(kernel_fit_0, fdata(k_xt)))
k_control_guys <- control_guys %>%
  mutate(y_0 = predict(kernel_fit_0),
         y_1 = predict(kernel_fit_1, fdata(k_xc)))

k_sim_pool <- smoothed_data %>%
  full_join(k_trt_guys %>%
              full_join(k_control_guys))

k_sim_id_data <- k_sim_pool %>%
  dplyr::select(id, a, y_1, y_0) %>%
  unique
k_sim_facts <- k_sim_id_data %>%
  summarise(Delta = mean(y_1[a == 1]) - mean(y_0[a==0]),
            Delta_S = mean(y_1[a==0]) - mean(y_0[a==0]),
            R = 1 - Delta_S/Delta)

g_trt_guys <- trt_guys %>%
  mutate(y_1 = predict(fgam_fit_1),
         y_0 = predict(fgam_fit_0, newdata = list(g_xc = g_xt), type = 'response'))
g_control_guys <- control_guys %>%
  mutate(y_0 = predict(fgam_fit_0),
         y_1 = predict(fgam_fit_1, newdata = list(g_xt = g_xc), type = 'response'))

g_sim_pool <- smoothed_data %>%
  full_join(g_trt_guys %>%
              full_join(g_control_guys))
g_sim_id_data <- g_sim_pool %>%
  dplyr::select(id, a, y_1, y_0) %>%
  unique
g_sim_facts <- g_sim_id_data %>%
  summarise(Delta = mean(y_1[a == 1]) - mean(y_0[a==0]),
            Delta_S = mean(y_1[a==0]) - mean(y_0[a==0]),
            R = 1 - Delta_S/Delta)

l_trt_guys <- trt_guys %>%
  mutate(y_1 = predict(lin_1),
         y_0 = predict(lin_0, newdata = list(l_xc = l_xt), type = 'response'))
l_control_guys <- control_guys %>%
  mutate(y_0 = predict(lin_0),
         y_1 = predict(lin_1, newdata = list(l_xt = l_xc), type = 'response'))

l_sim_pool <- smoothed_data %>%
  full_join(l_trt_guys %>%
              full_join(l_control_guys))
l_sim_id_data <- l_sim_pool %>%
  dplyr::select(id, a, y_1, y_0) %>%
  unique
l_sim_facts <- l_sim_id_data %>%
  summarise(Delta = mean(y_1[a == 1]) - mean(y_0[a==0]),
            Delta_S = mean(y_1[a==0]) - mean(y_0[a==0]),
            R = 1 - Delta_S/Delta)

sim_list <- list(k_sigma_1 = k_sigma_1,
                 k_sigma_0 = k_sigma_0,
                 k_sim_id_data = k_sim_id_data,
                 k_grid = k_grid,
                 k_sim_facts = k_sim_facts,
                 g_sigma_1 = g_sigma_1,
                 g_sigma_0 = g_sigma_0,
                 g_sim_id_data = g_sim_id_data,
                 g_grid = g_grid,
                 g_sim_facts = g_sim_facts,
                 l_sigma_1 = l_sigma_1,
                 l_sigma_0 = l_sigma_0,
                 l_sim_id_data = l_sim_id_data,
                 l_grid = l_grid,
                 l_sim_facts = l_sim_facts)
write_rds(sim_list, here('results/hiv-sim-tools.rds'))
