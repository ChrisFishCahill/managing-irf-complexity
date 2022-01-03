#----------------------------------------------------------------------
# develop harvest control rules for Alberta Walleye lakes
# Cahill & Walters 3 Jan 2022 
#----------------------------------------------------------------------

# load packages
library(tidyverse)
library(tidybayes)
library(purrr)

# install.packages("devtools")
# devtools::install_github("seananderson/ggsidekick")

# declare some variables 
ages <- 2:20
t <- 2000 # first survey year
max_a <- max(ages) # max age
rec_a <- min(ages) # age at recruitment
initial_yr <- t - max_a + rec_a - 2 
initial_yr_minus_one <- initial_yr - 1

# list the fits
list.files("fits/", full.names = TRUE)

# extract all the saved .stan fit names
paths <- dir("fits/", pattern = "\\.rds$")
paths <- paste0(getwd(), "/fits/", paths)
fits <- map(paths, readRDS) %>%
  set_names(paths)

# pick one
my_string <- names(fits)[1] 
fit <- fits[[which(names(fits)==my_string)]]

#----------------------------------------------------------------------
# extract things from some (identical) rows of the posterior 
n_draws <- 5 

devs <- fit %>% 
  spread_draws(R0, ar, br, sbr0_kick) %>%
  sample_draws(n_draws)

# index which draws were selected to preserve correlation structure
draw_idx <- unique(devs$.draw)

w_devs <- fit %>%
  spread_draws(w[year]) %>%
  filter(.draw %in% draw_idx) %>% # force same .draw to be taken as above
  mutate(year = year + initial_yr_minus_one) #make years 1980-2028

v_devs <- fit %>%
  spread_draws(v[period]) %>%
  filter(.draw %in% draw_idx)

F_devs <- fit %>%
  spread_draws(F_vec[year]) %>%
  filter(.draw %in% draw_idx) %>%
  mutate(year = year + initial_yr_minus_one)

# extract nta from stan
nta_stan <- fit %>%
  spread_draws(Nat_array[age, year]) %>%
  pivot_wider(names_from = age, values_from = Nat_array) %>%
  filter(.draw %in% draw_idx) %>%
  mutate(year = year + initial_yr_minus_one)

# which columns have ages: 
age_cols <- which(!is.na(str_extract(string = colnames(nta_stan), 
                         pattern = "[0-9]|10[0-9]") ))

# rename age columns to correct ages 2-20
colnames(nta_stan)[age_cols] = ages 

#----------------------------------------------------------------------
# ***N.B.***
# We need to preserve the historical frequency of weak and strong 
# year classes in our recruitment time series for the retrospective 
# analysis
# 
# Thus, we will select 1990-2015 for these lakes for a 25 yr recruitment 
# reference period as most FWIN surveys have information on recruitment
# back to 1990 
#----------------------------------------------------------------------

# initialize the simulation using BERTA posterior 

retro_initial_yr <- 1990 
retro_terminal_yr <- 2015
n_repeats <- 3 # 3 recruitment series repeats 

# 3 recruitment sequence repeats: 
n_sim_years <- length(retro_initial_yr:retro_terminal_yr)*n_repeats 

# extract the age structure corresponding to 1990   
nta_init <- nta_stan %>% filter(year == retro_initial_yr)

# extract w estimates from 1990-2015
w_devs <- w_devs %>%
  filter(year %in% retro_initial_yr:retro_terminal_yr)

# extract F estimates from 1990-2015
F_devs <- F_devs %>% 
  filter(year %in% retro_initial_yr:retro_terminal_yr)



