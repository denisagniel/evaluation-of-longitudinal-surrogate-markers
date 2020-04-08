#' ---
#' title: "Collect possible time grids"
#' output: github_document
#' ---
#' 
library(knitr)
opts_chunk$set(warning = FALSE, message = FALSE, cache = FALSE, fig.width = 7, fig.height = 7)

#'
#'
library(tidyverse)
library(here)
tst <- tidyr::expand_grid(
  left_trim = 1:5,
  right_trim = 1:20,
  size = seq(11, 51, by = 2)
)

tstl <- apply(tst, 1, function(x) {
  lt <- x[1]
  rt <- x[2]
  sz <- x[3]
  out <- seq(0 + lt, 52 - rt, length = sz) %>% round %>% unique
  out
}) %>% unique

write_rds(tstl, here('data/time-grids-for-cv.rds'))
