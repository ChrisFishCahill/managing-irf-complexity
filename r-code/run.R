#----------------------------------------------------------------------
# single lake version of BERTA run.R for AEP contract 
# Cahill 2021
#----------------------------------------------------------------------

# libraries
library(tidyverse)
library(rstan)

# show stan all cores
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)

# compile the model
m <- rstan::stan_model("src/BERTA_single_lake.stan", verbose = F)

# read in the data from ALL lakes used in Cahill et al. 2021
data <- readRDS("data/BERTA-wide-0-25.rds")
stocking <- readRDS("data/stocking_matrix_ha.rds")

# create function to run .stan model

get_fit <- function(which_lake = "pigeon lake",
                    rec_model = c("bev-holt", "ricker"),
                    cr_prior = c(6, 12),
                    n_iter = n_iter, n_chains = n_chains,
                    n_warmup = n_iter / 2,
                    data = data,
                    stocking = stocking, 
                    ...) {
  rec_model <- match.arg(rec_model)
  cat(
    crayon::green(
      clisymbols::symbol$tick
    ),
    fitted = "model fitted = ", which_lake, rec_model, 
    sep = " "
  )
  cat("\n")
  
  #filter the run data from all data, re-order it 
  run_data <- data %>% filter(name %in% which_lake)
  run_data <-
    within(run_data, lake <-
      as.numeric(interaction(
        run_data$WBID,
        drop = TRUE, lex.order = F
      )))
  run_data <- run_data[order(run_data$lake), ]
  
  # stocking stuff was run in different versions, now just for plotting:
  run_stocking <- stocking[which(rownames(stocking) %in% which_lake), ]
  
  # Add ten years of zero for short term projections
  proj_stock <- rep(0, 10) 
  run_stocking <- round(c(run_stocking, proj_stock)) #add to stocking data (for plots)
  
  # Set up the Rbar years
  suppressMessages(
    survey_yrs <- run_data %>%
      group_by(lake) %>%
      summarise(
        min_yr = min(year) + length(initial_yr:(t - 1)),
        max_yr = max(year) + length(initial_yr:(t - 1))
      ) %>%
      as.numeric()
  )
 survey_yrs <- survey_yrs[2:3]
  # summarize the life history relationships
  suppressMessages(
    life_hist <- run_data %>%
      group_by(lake) %>%
      summarize(
        a50 = unique(a50),
        vbk = unique(vbk),
        linf = unique(linf),
        wl_beta = unique(beta_wl)
      )
  )
  
  Fseq <- seq(from = 0.01, to = 1.0, by = 0.01)
  
  # declare the tagged data list for stan
  stan_data <- list(
    n_surveys = nrow(run_data),
    n_ages = length(Ages),
    n_obs = nrow(run_data) * length(Ages),
    n_years = length(initial_yr:2028),
    n_lakes = length(unique(run_data$lake)),
    caa = run_data[, which(colnames(run_data) %in% Ages)],
    prop_aged = run_data$p_aged,
    effort = run_data$effort,
    lake = run_data$lake,
    year = run_data$year + length(initial_yr:(t - 1)),
    ages = Ages,
    survey_yrs = survey_yrs,
    which_year = 1996 - initial_yr + 2, # which integer corresponds to year = 1997
    v_prior_early = 0.3,
    v_prior_late = 0.1,
    prior_sigma_v = c(0.1, 0.5),
    R0_mean = log(6),
    R0_sd = log(3),
    ar_sd = 0.1,
    prior_mean_w = 0,
    prior_sigma_w = 1.2,
    vbk = life_hist$vbk,
    linf = life_hist$linf,
    a50 = life_hist$a50,
    wl_beta = life_hist$wl_beta,
    lbar = 57.57, # From cahill et al. 2020
    M = 0.1,
    theta = 0.85, # Lorenzen M exponent
    phi = 2.02, # vulnerability parameter (nets)
    psi = 2, # vulnerability parameter (angling)
    G_bound = c(0, Inf),
    get_SSB_obs = 1L,
    obs_cv_prior = 0.15,
    SSB_penalty = 0,
    prior_sigma_G = 1,
    Rinit_ctl = 0,
    length_Fseq = length(Fseq),
    Fseq = Fseq,
    rec_model = ifelse(rec_model == "ricker", 0, 1),
    cr_prior = cr_prior
  )

  # create a function for start values
  vk <- c(0.3, 0.3)
  inits <- function() {
    list(
      v = jitter(vk, amount = 0.1),
      R0 = jitter(15, amount = 2),
      G = jitter(1, amount = 0.1),
      w = jitter(rep(0,
                     stan_data$n_years - 2), 
                 amount = 0.1),
      sigma_w = jitter(0.5, amount = 0.05),
      ar = jitter(0.5, amount = 0.01)
    )
  }

  #run the model
  fit <-
    rstan::sampling(
      m,
      data = stan_data,
      pars =
        c(
          "ar_mean_kick", "F_ratio", "Fmsy", "MSY",
          "G", "cr", "ar", "SPR", "br",
          "SBR", "sbr0_kick", "R0", "v", "SSB",
          "R2", "SSB_obs", "caa_pred", "b_ratio", "w"
        ),
      iter = n_iter,
      warmup = n_warmup,
      chains = n_chains,
      init = inits,
      control = list(
        adapt_delta = 0.999,
        max_treedepth = 15
      )
    )
  fit
}

# declare some model index values 
Ages <- 2:20
t <- 2000 # first survey year
max_a <- max(Ages)
rec_a <- min(Ages)
initial_yr <- t - max_a + rec_a - 2 
add_year <- initial_yr - 1

# declare HMC run parameters 
n_iter = 100
n_chains = 3
n_warmup = n_iter/2
names <- unique(data$name)

#----------------------------------------------------------------------
# test with a single model 

fit <- get_fit(which_lake = "buck lake", 
               rec_model = "bev-holt", 
               cr_prior = 6, 
               n_iter = n_iter, n_chains = n_chains, 
               n_warmup = n_warmup, 
               data=data, stocking=stocking)

#----------------------------------------------------------------------

#TODO: Run with all six lakes using some fancy functional programming stuff
