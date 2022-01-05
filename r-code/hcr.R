#----------------------------------------------------------------------
# Harvest control rules for Alberta Walleye lakes
# Cahill & Walters 3 Jan 2022
# TODO: assessing every few years code needs updating
#----------------------------------------------------------------------

# load packages
library(tidyverse)
library(tidybayes)
library(purrr)

# install.packages("devtools")
# devtools::install_github("seananderson/ggsidekick")

# declare some variables
ages <- 2:20
n_ages <- length(ages)
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

# set rec_ctl 
rec_ctl <- ifelse(grep("bh", my_string), "bh", "ricker")

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
               w_a_report[age], 
               M_a_report[age]) %>%
  filter(.draw %in% draw_idx)

leading_pars$age <- ages # correct ages

w_a <- leading_pars$w_a_report
f_a <- leading_pars$f_a_report
v_survey <- leading_pars$v_a_report
v_fish <- leading_pars$v_f_a_report
Lo <- leading_pars$Lo_report
M_a <- leading_pars$M_a_report

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

# hcr psedo-code
# for each cslope i 
#  for each bmin j 
#   for each draw k 
#    for each year in n_sim_years l 
#     [run model, record performance]
# end i, j k, l

# Run retrospective simulation for each cslope, bmin, draw, and simulation year
for (i in seq_along(c_slope_seq)) {
  c_slope <- c_slope_seq[i]
  for (j in seq_along(bmin_seq)) {
    b_lrp <- bmin_seq[j]
    set.seed(83) # challenge each bmin, cslope combo with same set of unpredictable rec seqs
    for (k in seq_len(n_draws)) {
      # k = 1
      # pick a single draw
      sub_post <- subset(sampled_post, sampled_post$.draw == unique(sampled_post$.draw)[k])
      
      # set leading parameters from sampled draw
      rec_a <- sub_post$ar[1]
      rec_b <- sub_post$br[1]
      Ro <- sub_post$Ro[1]
      sbo <- Ro * sbro
      vbo <- Ro * vbro
      
      wt <- sub_post$w
      rec_var <- 1.0 # 1.2 might be fun to try
      wt <- c(wt, wt * rec_var, wt * rec_var)
      
      # extract the initial age structure
      Ninit <- sub_post[
        which(sub_post$year == retro_initial_yr),
        which(colnames(sub_post) %in% ages)
      ] %>%
        slice() %>%
        unlist(., use.names = FALSE)
      
      # nta matrix
      nta <- matrix(NA, nrow = length(wt), ncol = length(ages))
      
      # SSB, Rpred vectors
      SSB <- Rpred <- vB_fish <- vB_survey <- rep(0, length(wt))
      
      nta[1, ] <- Ninit # initialize from posterior for retro_initial_yr
      
      # run age-structured model for sim years
      for (t in seq_len(n_sim_years)[-n_sim_years]) { # years 1 to (n_sim_year-1)
        SSB[t] <- sum(nta[t, ] * f_a * w_a)
        vB_fish[t] <- sum(nta[t, ] * v_fish * w_a)
        vB_survey[t] <- sum(nta[t, ] * v_survey * w_a)
        
        # set observed vb from "true" 
        obs_sd <- 0.1
        q_survey <- 1.0 #q_survey assumed to be 1.0 in Cahill et al. 2021
        vB_obs <- q_survey * vB_survey[t] * exp(obs_sd * (rnorm(1)) - 0.5 * (obs_sd)^2) #-0.5*(0.1)^2 corrects exponential effect on mean observation
        
        # Set up hcr
        TAC <- c_slope * (vB_obs - b_lrp)
        if (TAC < 0) {
          TAC <- 0
          Ut <- 0
        } else {
          Ut <- ifelse((TAC / vB_fish[t]) < 0.9, (TAC / vB_fish[t]), 0.9)
        }
        
        # stock-recruitment 
        if (rec_ctl == "ricker") { 
          Rpred[t] <- rec_a * SSB[t] * exp(-rec_b * SSB[t] + wt[t])
        }
        if (rec_ctl == "bh") { 
          Rpred[t] <- rec_a * SSB[t] * exp(wt[t]) / (1 + rec_b * SSB[t])
        }
        
        # update the age structure
        for (a in seq_len(n_ages)[-n_ages]) { #ages 2-19 
          nta[t + 1, a + 1] <- nta[t, a] * exp(-M_a[a])* (1 - Ut * v_fish[a])
        }
        
        # record performance metrics
        yield <- Ut * vB_fish[t]
        yield_array[i, j, t] <- yield_array[i, j, t] + yield
        tot_y[i, j] <- tot_y[i, j] + yield
        tot_u[i, j] <- tot_u[i, j] + yield^0.3
        prop_below[i, j] <- prop_below[i, j] + ifelse(SSB[t] < 0.1 * sbo, 1, 0) # SSBs falling below 0.1*sbo
        TAC_zero[i, j] <- TAC_zero[i, j] + ifelse(TAC == 0, 1, 0)
        
        # set rec value for next t
        nta[t + 1, 1] <- Rpred[t]
      }
    }
  }
}
  