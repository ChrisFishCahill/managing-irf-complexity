# develop harvest control rules for Alberta Walleye lakes
# Cahill & Walters 2022 
# load packages

library(tidyverse)
library(tidybayes)
library(purrr)

# install.packages("devtools")
# devtools::install_github("seananderson/ggsidekick")

# list the fits
list.files("fits/", full.names = TRUE)

# extract all the saved .stan fit names
paths <- dir("fits/", pattern = "\\.rds$")
paths <- paste0(getwd(), "/fits/", paths)
fits <- map(paths, readRDS) %>%
  set_names(paths)

my_string <- names(fits)[10] #lac ste. anne bh cr = 6
fit <- fits[[which(names(fits)==my_string)]]

#----------------------------------------------------------------------
# ***N.B.***
# We need to preserve the historical frequency of weak and strong 
# year classes in our recruitment time series for our retrospective 
# analysis
# Thus, we will select 1990-2015 for these lakes for a 25 yr ref period
#----------------------------------------------------------------------

# Extract things from some (identical) rows of the posterior to 
# preserve correlation structure

n_draws <- 5 

devs <- fit %>% 
  spread_draws(R0, ar, br, sbr0_kick) %>%
  sample_draws(n_draws)

# index which draws were selected to ensure same draws for each simulation
draw_idx <- unique(devs$.draw)

w_devs <- fit %>%
  spread_draws(w[year]) %>%
  filter(.draw %in% draw_idx) # force same .draw to be taken as above

v_devs <- fit %>%
  spread_draws(v[period]) %>%
  filter(.draw %in% draw_idx)

F_devs <- fit %>%
  spread_draws(F_vec[year]) %>%
  filter(.draw %in% draw_idx)

# extract nta from stan
nta_stan <- fit %>%
  spread_draws(Nat_array[age, year]) %>%
  pivot_wider(names_from = age, values_from = Nat_array) %>%
  filter(.draw %in% draw_idx) 
