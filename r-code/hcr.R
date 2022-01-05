#----------------------------------------------------------------------
# Harvest control rules for Alberta Walleye lakes
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
fit <- fits[[which(names(fits) == my_string)]]

#----------------------------------------------------------------------
# extract things from some (identical) rows of the posterior
n_draws <- 1

devs <- fit %>%
  spread_draws(Ro, ar, br, sbro_report) %>%
  sample_draws(n_draws)

# index which draws were selected to preserve correlation structure
draw_idx <- unique(devs$.draw)

w_devs <- fit %>%
  spread_draws(w[year]) %>%
  filter(.draw %in% draw_idx) %>% # force same .draw to be taken as above
  mutate(year = year + initial_yr_minus_one) # make years 1980-2028

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
age_cols <- which(!is.na(str_extract(
  string = colnames(nta_stan),
  pattern = "[0-9]|10[0-9]"
)))

# rename age columns to correct ages 2-20
colnames(nta_stan)[age_cols] <- ages

#----------------------------------------------------------------------
#                             ***N.B.***
# We need to preserve the historical frequency of weak and strong
# year classes in our recruitment time series for the retrospective
# analysis
#
# Thus, we will select 1990-2015 for these lakes for a 26 yr recruitment
# reference period as most FWIN surveys have information on recruitment
# to approximately 1990
#----------------------------------------------------------------------

# initialize the retrospective simulation using BERTA posterior

retro_initial_yr <- 1990
retro_terminal_yr <- 2015

# extract the age structure corresponding to 1990
nta_init <- nta_stan %>% filter(year == retro_initial_yr)

# extract w estimates from 1990-2015
w_devs <- w_devs %>%
  filter(year %in% retro_initial_yr:retro_terminal_yr)

# extract F estimates from 1990-2015
F_devs <- F_devs %>%
  filter(year %in% retro_initial_yr:retro_terminal_yr)

# join all those devs into one big tibble called "sampled_post"
sampled_post <- left_join(w_devs, nta_stan,
  join_by = c(.draw, year)
)

sampled_post <- left_join(sampled_post, F_devs,
  join_by = c(.draw, year)
)

sampled_post <- left_join(sampled_post, devs,
  join_by = .draw
) %>% arrange(.draw)

#----------------------------------------------------------------------
# set up cslope, bmin sequences and performance metric output matrices

# Extract MAP estimate of Ro from posterior (for indexing bmin)
Ro_summary <- rstan::summary(fit, pars = "Ro")$summary
Ro_map <- Ro_summary[, "mean"]
Ro_map 

#---------------------------------------------------------
# extrac leading parameters from stan to get vbro
leading_pars <- fit %>%
  spread_draws(Lo_report[age], 
               l_a_report[age], 
               v_a_report[age],
               v_f_a_report[age], 
               f_a_report[age], 
               w_a_report[age]) %>%
  filter(.draw %in% draw_idx)

leading_pars$age <- ages # correct ages

w_a <- leading_pars$w_a_report
f_a <- leading_pars$f_a_report
v_survey <- leading_pars$v_a_report
v_fish <- leading_pars$v_f_a_report
Lo <- leading_pars$Lo_report

vbro <- sum(Lo * v_survey * w_a)
sbro <- sum(f_a*Lo)

vbro

# recruitment sequence repeats:
n_repeats <- 3 
n_sim_years <- length(retro_initial_yr:retro_terminal_yr) * n_repeats

# c_slope and bmin ranges for harvest control rule
c_slope_seq <- seq(from = 0.05, to = 1.0, by = 0.05)
bmin_seq <- seq(from = 0, to = 1.0 * Ro_map * vbro, length.out = length(c_slope_seq))
tot_y <- tot_u <- prop_below <- TAC_zero <- matrix(0, nrow = length(c_slope_seq), ncol = length(bmin_seq))
yield_array <- array(0, dim = c(length(c_slope_seq), length(bmin_seq), n_sim_years))

#----------------------------------------------------------------------
k = sampled_post$.draw[1] # pick a draw
sub_post <- subset(sampled_post, sampled_post$.draw == k)
wt <- sub_post$w

rec_var <- 1.0 # 1.2 might be fun to try
wt <- c(wt, wt * rec_var, wt * rec_var)
df <- data.frame(wt = wt, sim_yrs = 1:n_sim_years)

df %>%
ggplot(aes(x=sim_yrs, y=wt)) + 
  geom_point() + 
  geom_line() + 
  xlab("Year of Simulation") + 
  ylab("Recruitment Anomaly ln(wt)") + 
  ggsidekick::theme_sleek()


#run the hcr psedo-code
#for each cslope i 
# for each bmin j 
#  for each draw k 
#   for each year in n_sim_years l 
#    [run model, record performance]
# end i, j k, l
