#----------------------------------------------------------------------
# single lake version of BERTA run.R for AEP contract 
# Cahill 2021
#----------------------------------------------------------------------

# libraries
library(tidyverse)
library(rstan)
library(purrr)
library(furrr)
library(future) 
library(tidybayes)

# detect number of CPUs on current host 
options(mc.cores = parallel::detectCores())

# eliminate redundant compilations
rstan::rstan_options(auto_write = TRUE)

# compile the model
m <- rstan::stan_model("src/BERTA_single_lake.stan", verbose = F)

# read in the data from all lakes used in Cahill et al. 2021
data <- readRDS("data/BERTA-wide-0-25.rds")
stocking <- readRDS("data/stocking_matrix_ha.rds")

#----------------------------------------------------------------------
# create function to run .stan model

get_fit <- function(which_lake = "pigeon lake",
                    rec_ctl = c("bev-holt", "ricker"),
                    cr_prior = c(6, 12),
                    n_iter = n_iter, n_chains = n_chains,
                    n_warmup = n_iter / 2,
                    ...) {
  rec_ctl <- match.arg(rec_ctl)
  cat(
    crayon::green(
      clisymbols::symbol$tick
    ),
    fitted = "model fitted = ", which_lake, rec_ctl, 
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
    n_ages = length(ages),
    n_obs = nrow(run_data) * length(ages),
    n_years = length(initial_yr:2028),
    n_lakes = length(unique(run_data$lake)),
    caa = run_data[, which(colnames(run_data) %in% ages)],
    prop_aged = run_data$p_aged,
    effort = run_data$effort,
    lake = run_data$lake,
    year = run_data$year + length(initial_yr:(t - 1)),
    ages = ages,
    survey_yrs = survey_yrs,
    which_year = 1996 - initial_yr + 2, # which integer corresponds to year = 1997
    v_prior_early = 0.3,
    v_prior_late = 0.1,
    prior_sigma_v = c(0.1, 0.5),
    Ro_mean = log(6),
    Ro_sd = log(3),
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
    rec_ctl = ifelse(rec_ctl == "ricker", 0, 1),
    cr_prior = cr_prior
  )

  # create a function for start values
  vk <- c(0.3, 0.3)
  inits <- function() {
    list(
      v = jitter(vk, amount = 0.1),
      Ro = jitter(15, amount = 2),
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
          "F_ratio", "Fmsy", "MSY",
          "G", "cr", "ar", "SPR", "br", "Nat_array",
          "SBR", "Ro", "v", "F_vec", "SSB",
          "R2", "SSB_obs", "caa_pred", "b_ratio", "w", 
          # report stuff: 
          "sbro_report", "ar_mean_report", "l_a_report", 
          "Lo_report", "v_a_report", "v_f_a_report", 
          "f_a_report", "w_a_report"
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
  # create name and save .rds files for each run 
  if(rec_ctl=="ricker"){
    my_name <- paste0(which_lake, "_ricker.rds")
  }
  if(rec_ctl=="bev-holt"){
    my_name <- paste0(which_lake, "_bh.rds")
  }
  stan_file <- "fits/"
  stan_file <- str_c(stan_file, my_name)
  stan_file <- stan_file %>% gsub(" ", "_", .)
  if(cr_prior == 6){
    stan_file <- stan_file %>% gsub(".rds", "_cr_6.rds", .)
  }
  if(cr_prior == 12){
    stan_file <- stan_file %>% gsub(".rds", "_cr_12.rds", .)
  }
  if (file.exists(stan_file)) {
    return(NULL)
  } else {
    saveRDS(fit, file = stan_file)
  }
}

#----------------------------------------------------------------------
# declare some indeces for stan model  
ages <- 2:20
t <- 2000 # first survey year
max_a <- max(ages)
rec_a <- min(ages)
initial_yr <- t - max_a + rec_a - 2 
initial_yr_minus_one <- initial_yr - 1

# declare HMC run parameters 
n_iter = 2000
n_chains = 4
n_warmup = n_iter/2

#----------------------------------------------------------------------
# test with a single lake / stock-recruitment function  

fit <- get_fit(which_lake = "lac ste. anne",
               rec_ctl = "bev-holt",
               cr_prior = 6,
               n_iter = n_iter, n_chains = n_chains,
               n_warmup = n_warmup)

#----------------------------------------------------------------------
# naughty functional programming black magjicks  

contract_lakes <- c("lac ste. anne", "baptiste lake", 
                    "pigeon lake", "calling lake", 
                    "moose lake", "lake newell"
)

to_fit <- tibble(which_lake = contract_lakes)
to_fit$n_iter <- n_iter
to_fit$n_chains <- n_chains
to_fit$n_warmup <- n_warmup
to_fit$rec_ctl <- "ricker"
to_fit$cr_prior <- 6

to_fit2 <- to_fit
to_fit2$rec_ctl <- "bev-holt"
to_fit <- rbind(to_fit, to_fit2)

to_fit3 <- to_fit
to_fit3$cr_prior <- 12
to_fit <- rbind(to_fit, to_fit3)

# set up parallel processing plan 
future::plan(multisession)

# Run models and save fits -- 5.5 hours on my laptop
system.time({ 
  future_pwalk(to_fit, get_fit, 
               .options = furrr_options(seed = TRUE)
  ) 
})

# future_pwalk looks like magic, is from the purrr/furrr family of libraries
# see https://adv-r.hadley.nz/fp.html sections on functional programming 

# pseudo-code for another more familiar way:
# for each row of to_fit
# run get_fit using parameters from single row of to fit
# save results
# end for-loop

# remove all fits? 
# do.call(file.remove, list(list.files("fits/", full.names = TRUE)))

#----------------------------------------------------------------------
# some diagnostics -- can do this for any fit

which_file <- list.files("fits/", full.names = TRUE)[1]
fit <- readRDS(which_file)
fit 
shinystan::launch_shinystan(fit)

#-----------------------------------------------------------------------
# end 
