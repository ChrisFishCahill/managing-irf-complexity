#----------------------------------------------------------------------
# Harvest control rules for Alberta Walleye lakes
# Cahill & Walters 3 Jan 2022
#----------------------------------------------------------------------
# load packages
library(tidyverse)
library(tidybayes)
library(purrr)
library(future)
library(furrr)
# install.packages("devtools")
# devtools::install_github("seananderson/ggsidekick")

#--------------------------------------------------------------------
# write a function to take BERTA posteriors and run MSE

get_hcr <- function(which_lake = "lac_ste._anne", ass_int = 1,
                    sd_survey = 0.05, d_mort = 0.3,
                    rule = c("linear", "rectilinear")) {
  #--------------------------------------------------------------------
  # initialize stuff
  #--------------------------------------------------------------------
  n_draws <- hcr_pars$n_draws
  n_repeats <- hcr_pars$n_repeats
  retro_initial_yr <- hcr_pars$retro_initial_yr
  retro_terminal_yr <- hcr_pars$retro_terminal_yr
  n_sim_yrs <- hcr_pars$n_sim_yrs
  ages <- hcr_pars$ages
  n_ages <- length(ages)
  initial_yr <- hcr_pars$initial_yr
  initial_yr_minus_one <- hcr_pars$initial_yr_minus_one
  ret_a <- hcr_pars$ret_a
  Ut_overall <- hcr_pars$Ut_overall
  sbo_prop <- hcr_pars$sbo_prop
  rho <- hcr_pars$rho
  sd_wt <- hcr_pars$sd_wt
  psi_wt <- hcr_pars$psi_wt
  n_historical_yrs <- length(retro_initial_yr:retro_terminal_yr)
  grid_size <- hcr_pars$grid_size

  #--------------------------------------------------------------------
  # subset lake-specific posterior from all fits
  #--------------------------------------------------------------------
  lake_str <- gsub(" ", "_", which_lake)
  fit_idx <- grep(lake_str, names(fits))
  fit <- fits[[fit_idx]]
  rec_ctl <- ifelse(grepl("bh", names(fits)[fit_idx]), "bh", "ricker")

  #--------------------------------------------------------------------
  # Get Bmay, Umay for rectilinear HCR (based on BH CR 6 for now)
  #--------------------------------------------------------------------
  if (rule == "rectilinear") {
    B_may <- may_data %>%
      filter(lake == gsub("_", " ", which_lake)) %>%
      select(c("BMAY"))
    B_may <- B_may$BMAY

    U_may <- may_data %>%
      filter(lake == gsub("_", " ", which_lake)) %>%
      select(c("UMAY"))
    U_may <- U_may$UMAY
  }

  # extract estimated and derived parameters from BERTA
  devs <- fit %>%
    spread_draws(Ro, ln_ar, br) %>%
    sample_draws(n_draws)

  draw_idx <- unique(devs$.draw) # which posterior rows did sample_draws() take?

  w_devs <- fit %>%
    spread_draws(w[year]) %>%
    filter(.draw %in% draw_idx) %>% # force same .draw to be taken as above
    mutate(year = year + initial_yr_minus_one) # make years 1980-2028

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

  # extract MAP estimate of Ro from posterior (for indexing bmin)
  Ro_summary <- rstan::summary(fit, pars = "Ro")$summary
  Ro_map <- Ro_summary[, "mean"]

  # extract the age structure corresponding to 1990
  nta_init <- nta_stan %>% filter(year == retro_initial_yr)

  # extract w estimates from 1990-2015
  w_devs <- w_devs %>%
    filter(year %in% retro_initial_yr:retro_terminal_yr)

  # join all those devs into one big tibble called "post"
  suppressMessages(
    post <- left_join(w_devs, nta_stan,
      join_by = c(.draw, year)
    )
  )

  suppressMessages(
    post <- left_join(post, devs,
      join_by = .draw
    ) %>% arrange(.draw)
  )

  # extract leading parameters from stan to get vbro
  leading_pars <- fit %>%
    spread_draws(
      Lo_report[age],
      l_a_report[age],
      v_a_report[age],
      v_f_a_report[age],
      f_a_report[age],
      w_a_report[age],
      M_a_report[age]
    ) %>%
    filter(.draw %in% draw_idx[1]) # indexed on 1 bc not changing
  leading_pars$age <- ages # correct ages
  w_a <- leading_pars$w_a_report
  f_a <- leading_pars$f_a_report

  v_survey <- leading_pars$v_a_report
  v_fish <- leading_pars$v_f_a_report
  Lo <- leading_pars$Lo_report
  M_a <- leading_pars$M_a_report
  vbro <- sum(Lo * v_survey * w_a)
  sbro <- sum(f_a * Lo)

  # put some extra stuff in post
  post$sbro <- sbro
  post$Ro_map <- Ro_map
  post$vbro <- vbro # biomass vulnerable per recruit *surveys*
  post$vbo <- sum(Ro_map * Lo * v_fish * w_a) # biomass vulnerable to *fishing*

  #----------------------------------------------------------------------
  # set up cslope, bmin sequences and performance metric output matrices
  #----------------------------------------------------------------------
  bmin_max_value <- Ro_map * sum(Lo * v_fish * w_a) # ifelse(Ro_map * vbro < 20, Ro_map * vbro, 20)
  if (which_lake == "calling lake") {
    bmin_max_value <- 100
  } else {
    bmin_max_value <- ifelse(bmin_max_value > 75, 75, bmin_max_value)
  }
  c_slope_seq <- seq(from = 0.01, to = 1.0, length.out = round(grid_size / 2))
  bmin_seq_low <- seq(from = 0, to = 20, length.out = round(grid_size / 2))
  bmin_seq_high <- seq(from = 20.5, to = bmin_max_value, length.out = grid_size - length(bmin_seq_low))
  bmin_seq <- c(bmin_seq_low, bmin_seq_high)

  if (rule == "rectilinear") {
    # only considering one rectilinear rule, so set loops to one
    c_slope_seq <- 1
    bmin_seq <- 1
  }

  tot_y <- tot_u <- prop_below <- TAC_zero <-
    matrix(0, nrow = length(c_slope_seq), ncol = length(bmin_seq))
  yield_array <- vB_fish_array <-
    array(0, dim = c(length(c_slope_seq), length(bmin_seq), n_sim_yrs))
  wt_seqs <- matrix(NA, nrow = n_sim_yrs, ncol = n_draws)

  #----------------------------------------------------------------------
  # run retrospective simulation for each cslope, bmin, draw, and sim yr
  #----------------------------------------------------------------------
  for (i in seq_along(c_slope_seq)) {
    c_slope <- c_slope_seq[i]
    for (j in seq_along(bmin_seq)) {
      b_lrp <- bmin_seq[j]
      set.seed(24) # challenge each bmin, cslope combo with same set of rec seqs
      wt_re_mat <- matrix(rnorm(n = n_sim_yrs * n_draws, mean = 0, sd = sd_wt), nrow = n_sim_yrs, ncol = n_draws) # generate random deviates
      for (k in seq_len(n_draws)) {
        # pick a single draw
        sub_post <- subset(post, post$.draw == unique(post$.draw)[k])

        # set leading parameters from sampled draw
        rec_a <- exp(sub_post$ln_ar[1])
        rec_b <- sub_post$br[1]
        Ro <- sub_post$Ro[1]
        sbo <- Ro * sbro
        vbo <- Ro * vbro

        # repeat the historical recruitment series
        wt_historical <- sub_post$w
        wt_bar <- mean(wt_historical)
        wt_historical <- wt_historical - wt_bar # correct nonzero wt over initialization period
        wt_historical <- rep(wt_historical, n_repeats) # 8 x 26 year sequence of values
 
        # set the random effect vector
        wt_re <- wt_re_mat[, k]

        # generate auto-correlated w(t)'s
        wt_sim <- wt <- rep(NA, length(wt_historical))
        wt_sim[1] <- wt_re[1] # initialize the process for t = 1
        
        # set wt = BERTA estimated values for yrs 1-26
        wt[1:n_historical_yrs] <- wt_historical[1:n_historical_yrs]

        # create autoregressive wt_sim[t]
        for (t in seq_len(n_sim_yrs)[-n_sim_yrs]) { # t = 1 to 207
          wt_sim[t + 1] <- rho * wt_sim[t] + wt_re[t + 1]
        }

        # calculate wt differently for yrs 26 +
        for (t in n_historical_yrs:n_sim_yrs) { # t = 26 to 208
          wt[t] <- psi_wt * wt_historical[t] + (1 - psi_wt) * wt_sim[t]
        }

        # extract the initial age structure
        Ninit <- sub_post[
          which(sub_post$year == retro_initial_yr),
          which(colnames(sub_post) %in% ages)
        ] %>%
          slice() %>%
          unlist(., use.names = FALSE)

        Rinit_yr_2 <- sub_post[
          which(sub_post$year == retro_initial_yr + 1),
          which(colnames(sub_post) == 2)
        ] %>%
          slice() %>%
          unlist(., use.names = FALSE)

        # nta matrix
        nta <- matrix(NA, nrow = length(wt) + 2, ncol = length(ages)) # recruit @ age 2

        # SSB, Rpred, vulnerable biomass vectors
        SSB <- Rpred <- vB_fish <- vB_survey <- rep(0, length(wt))
        nta[1, ] <- Ninit # initialize from posterior for retro_initial_yr
        nta[2, 1] <- Rinit_yr_2
        t_last_ass <- 1 - ass_int # when was last survey / assessment (initialize)

        # run age-structured model for sim years
        for (t in seq_len(n_sim_yrs)[-n_sim_yrs]) { # years 1 to (n_sim_year-1)
          SSB[t] <- sum(nta[t, ] * f_a * w_a)
          vB_fish[t] <- sum(nta[t, ] * v_fish * w_a * ret_a) # without ret_a, overestimate yield
          vB_survey[t] <- sum(nta[t, ] * v_survey * w_a)
          if (t - t_last_ass == ass_int) { # assess every ass_int yrs
            t_last_ass <- t
            # note -0.5*(0.1)^2 corrects exponential effect on mean observation:
            vB_obs <- vB_fish[t] * exp(sd_survey * (rnorm(1)) - 0.5 * (sd_survey)^2)

            if (rule == "linear") {
              TAC <- c_slope * (vB_obs - b_lrp)
            }

            if (rule == "rectilinear") {
              b_lrp <- 0.4 * B_may
              u_lrp <- 0.8 * B_may
              Ut <- U_may * (vB_obs - b_lrp) / (u_lrp - b_lrp)
              if (Ut < 0) {
                Ut <- 0
              }
              if (Ut > U_may) {
                Ut <- U_may
              }
              TAC <- Ut * vB_obs
            }

            if (TAC < 0) {
              TAC <- 0
            }
            Ut <- ifelse((TAC / vB_fish[t]) < Ut_overall, (TAC / vB_fish[t]), Ut_overall)
          }
          rett <- ifelse(Ut / Ut_overall <= 1.0, Ut / Ut_overall, 1.0) # cap rett annual retention proportion at 1.0
          if (any(rett * ret_a > 1)) {
            message("rett*ret_a yielded values > 1.0! \ncalculations cannot be trusted!")
            break
          }

          # stock-recruitment
          if (rec_ctl == "ricker") {
            Rpred[t] <- rec_a * SSB[t] * exp(-rec_b * SSB[t] + wt[t])
          }
          if (rec_ctl == "bh") {
            Rpred[t] <- rec_a * SSB[t] * exp(wt[t]) / (1 + rec_b * SSB[t])
          }

          # update the age structure
          for (a in seq_len(n_ages)[-n_ages]) { # a = 2-19
            nta[t + 1, a + 1] <- nta[t, a] * exp(-M_a[a]) *
              (1 - Ut_overall * v_fish[a] * (ret_a[a] * rett + (1 - ret_a[a] * rett) * d_mort))
          }

          # set rec value for t + 2
          nta[t + 2, 1] <- Rpred[t]

          # record performance metrics
          yield <- Ut_overall * rett * vB_fish[t]
          yield_array[i, j, t] <- yield_array[i, j, t] + yield
          vB_fish_array[i, j, t] <- vB_fish_array[i, j, t] + vB_fish[t]
          tot_y[i, j] <- tot_y[i, j] + yield
          tot_u[i, j] <- tot_u[i, j] + yield^0.3 #pp = exponent term 
          prop_below[i, j] <- prop_below[i, j] + ifelse(SSB[t] < sbo_prop * sbo, 1, 0)
          TAC_zero[i, j] <- TAC_zero[i, j] + ifelse(rett == 1, 1, 0)
        }
        wt_seqs[, k] <- wt
      }
    }
  }
  #--------------------------------------------------------------------
  # process simulation output - get annual values & return a list
  #--------------------------------------------------------------------
  tot_y <- tot_y / (n_draws * n_sim_yrs)
  tot_u <- tot_u / (n_draws * n_sim_yrs)
  prop_below <- prop_below / (n_draws * n_sim_yrs)
  TAC_zero <- TAC_zero / (n_draws * n_sim_yrs)

  row.names(tot_y) <- row.names(tot_u) <-
    row.names(prop_below) <- row.names(TAC_zero) <- c_slope_seq
  colnames(tot_y) <- colnames(tot_u) <-
    colnames(prop_below) <- colnames(TAC_zero) <- round(bmin_seq, 2)

  yield_array <- yield_array / n_draws # expected yield seqs for each cslope, bmin
  vB_fish_array <- vB_fish_array / n_draws # expected vul bio seq for each i,j

  row_idx <- which(tot_y == max(tot_y), arr.ind = TRUE)[1]
  col_idx <- which(tot_y == max(tot_y), arr.ind = TRUE)[2]
  MSY_yields <- yield_array[row_idx, col_idx, ] # best MSY yields
  MSY_vB_fish <- vB_fish_array[row_idx, col_idx, ] # vB for MSY yields

  row_idx <- which(tot_u == max(tot_u), arr.ind = TRUE)[1]
  col_idx <- which(tot_u == max(tot_u), arr.ind = TRUE)[2]
  HARA_yields <- yield_array[row_idx, col_idx, ] # best HARA yields
  HARA_vB_fish <- vB_fish_array[row_idx, col_idx, ]

  hcr_sim_list <- list(
    "tot_y" = tot_y, "tot_u" = tot_u,
    "prop_below" = prop_below, "TAC_zero" = TAC_zero,
    "yield_array" = yield_array, "MSY_yields" = MSY_yields,
    "HARA_yields" = HARA_yields, "vB_fish_array" = vB_fish_array,
    "MSY_vB_fish" = MSY_vB_fish, "HARA_vB_fish" = HARA_vB_fish,
    "post" = post, "leading_pars" = leading_pars, "wt_seqs" = wt_seqs
  )
  # create name and save .rds files for each run
  file_name <- str_extract(
    string = names(fits)[fit_idx],
    pattern = "(?<=fits/).*(?=.rds)"
  )
  file_name <- paste0(
    "sims/", file_name, "_hcr", "_ass_int_", ass_int,
    "_sd_", sd_survey, "_d_mort_", d_mort,
    "_rule_", rule, ".rds"
  )
  if (file.exists(file_name)) {
    return(NULL)
  } else {
    saveRDS(hcr_sim_list, file = file_name)
  }
}

#----------------------------------------------------------------------
# declare some values for the simulations
#----------------------------------------------------------------------

n_draws <- 30
n_repeats <- 8 # recruitment repeats
retro_initial_yr <- 1990 # initial year for retrospective analysis
retro_terminal_yr <- 2015
n_sim_yrs <- length(retro_initial_yr:retro_terminal_yr) * n_repeats
ages <- 2:20
t <- 2000 # first survey year
max_a <- max(ages) # max age
recruit_a <- min(ages) # age at recruitment
initial_yr <- t - max_a + recruit_a - 2 # 1980 = BERTA initial yr
initial_yr_minus_one <- initial_yr - 1
ah_ret <- 5
sd_ret <- 1
ret_a <- 1 / (1 + exp(-(ages - ah_ret) / sd_ret)) # retention by age vector
Ut_overall <- 0.5 # annual capture rate of fully vulnerable fish
sbo_prop <- 0.1 # performance measure value to see if SSB falls below sbo_prop*sbo
rho <- 0.6 # correlation for recruitment terms
sd_wt <- 1.1 # std. dev w(t)'s
psi_wt <- 0.5 # weighting multiplier for wt_historical vs. wt_sim, aka "wthistory" in spreadsheets
grid_size <- 75 # how many bmins, round(grid_size/2) gives how many cslope

# put it all in a tagged list
hcr_pars <- list(
  "n_draws" = n_draws,
  "n_repeats" = n_repeats,
  "retro_initial_yr" = retro_initial_yr,
  "retro_terminal_yr" = retro_terminal_yr,
  "n_sim_yrs" = n_sim_yrs,
  "ages" = ages,
  "initial_yr" = initial_yr,
  "initial_yr_minus_one" = initial_yr_minus_one,
  "ret_a" = ret_a,
  "Ut_overall" = Ut_overall,
  "sbo_prop" = sbo_prop,
  "rho" = rho,
  "sd_wt" = sd_wt,
  "psi_wt" = psi_wt,
  "grid_size" = grid_size
)

#----------------------------------------------------------------------
# read in the stan fits and run retrospective mse
#----------------------------------------------------------------------

# extract some saved .stan fit names
paths <- dir("fits/", pattern = "\\.rds$")
paths <- paths[grep("bh_cr_6", paths)]
paths <- paste0(getwd(), "/fits/", paths)

fits <- map(paths, readRDS) %>%
  set_names(paths)

which_lakes <- str_extract(
  string = paths,
  "(?<=fits/).*(?=_bh|_ricker)"
)

# things to simulate hcrs across
ass_ints <- c(1, 3, 5, 10) # how often to assess / run FWIN
sd_surveys <- c(0.05, 0.4) # survey sd
d_morts <- c(0.03, 0.15, 0.3) # discard mortalities

# use expand.grid() and distinct() to get all possible combinations
to_sim <- expand.grid(which_lakes, ass_ints, sd_surveys, d_morts)
names(to_sim) <- c("which_lake", "ass_int", "sd_survey", "d_mort")
to_sim <- to_sim %>% distinct()
glimpse(to_sim)
to_sim$rule <- "linear"

#----------------------------------------------------------------------
# one lake at a time
# [you must understand the Tao of Programming before transcending structure]
#
# contract_lakes <- c(
#   "lac ste. anne", "baptiste lake",
#   "pigeon lake", "calling lake",
#   "moose lake", "lake newell"
# )
get_hcr(which_lake = "lac ste. anne", ass_int = 1,
        sd_survey = 0.05, d_mort = 0.03, rule = "linear")
# purrr
# to_sim <- tibble(which_lake = "lac ste. anne", ass_int = 1,
#                  sd_survey = 0.05, d_mort = 0.03, rule = "linear")
# pwalk(to_sim, get_hcr)
#----------------------------------------------------------------------
# run the linear policies 

# if you run all these your gov computer may explode
if (FALSE) {
  options(future.globals.maxSize = 8000 * 1024^2) # 8 GB
  future::plan(multisession)
  system.time({
    future_pwalk(to_sim, get_hcr,
      .options = furrr_options(seed = TRUE)
    )
  })
}

# remove all fits?
# do.call(file.remove, list(list.files("sims/", full.names = TRUE)))
#----------------------------------------------------------------------
# Try to run the rectilinear rules
may_data <- read.csv("data/may_yield_calcs.csv")

# run <- get_hcr(which_lake = "pigeon lake", ass_int = 1,
#                 sd_survey = 0.4, d_mort = 0.3, rule = "rectilinear")

to_sim <- data.frame(
  which_lake = which_lakes, ass_int = 1,
  sd_survey = 0.4, d_mort = 0.3, rule = "rectilinear"
)
pwalk(to_sim, get_hcr)

# remove rectilinear sims?
do.call(file.remove, list(list.files("sims/", pattern = "rectilinear", full.names = TRUE)))
#----------------------------------------------------------------------
# end
