#----------------------------------------------------------------------
# stripped down hcr function to play around with results faster
#----------------------------------------------------------------------
# load packages
library(tidyverse)
library(tidybayes)
library(purrr)
library(future)
library(furrr)
# library(rstan)
# install.packages("devtools")
# devtools::install_github("seananderson/ggsidekick")

#--------------------------------------------------------------------

get_hcr <- function(which_lake = "lac ste. anne", ass_int = 1) {
  #--------------------------------------------------------------------
  # initialize stuff
  #--------------------------------------------------------------------
  n_draws <- hcr_pars$n_draws
  rec_var <- hcr_pars$rec_var
  n_repeats <- hcr_pars$n_repeats
  retro_initial_yr <- hcr_pars$retro_initial_yr
  retro_terminal_yr <- hcr_pars$retro_terminal_yr
  n_sim_yrs <- hcr_pars$n_sim_yrs
  ages <- hcr_pars$ages
  n_ages <- length(ages)
  initial_yr <- hcr_pars$initial_yr
  initial_yr_minus_one <- hcr_pars$initial_yr_minus_one
  d_mort <- hcr_pars$d_mort
  ret_a <- hcr_pars$ret_a
  Ut_overall <- hcr_pars$Ut_overall
  sd_survey <- hcr_pars$sd_survey
  q_survey <- hcr_pars$q_survey
  sbo_prop <- hcr_pars$sbo_prop
  
  #--------------------------------------------------------------------
  # subset lake-specific posterior from all fits
  #--------------------------------------------------------------------
  lake_str <- gsub(" ", "_", which_lake)
  fit_idx <- grep(lake_str, names(fits))
  fit <- fits[[fit_idx]]
  rec_ctl <- ifelse(grepl("bh", names(fits)[fit_idx]), "bh", "ricker")

  # extract estimated and derived parameters from BERTA
  devs <- fit %>%
    spread_draws(Ro, ln_ar, br, cr) %>%
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
  #----------------------------------------------------------------------
  # set up cslope, bmin sequences and performance metric output matrices
  #----------------------------------------------------------------------
  # bmin_max_value <- Ro_map * sum(Lo * v_fish * w_a) # ifelse(Ro_map * vbro < 20, Ro_map * vbro, 20)
  # bmin_max_value <- ifelse(bmin_max_value > 100, 100, bmin_max_value)
  # c_slope_seq <- seq(from = 0.05, to = 1.0, by = 0.05)
  # bmin_seq <- seq(from = 0, to = bmin_max_value, length.out = length(c_slope_seq))
  # tot_y <- tot_u <- prop_below <- TAC_zero <-
  #   matrix(0, nrow = length(c_slope_seq), ncol = length(bmin_seq))
  # yield_array <- vB_fish_array <-
  #   array(0, dim = c(length(c_slope_seq), length(bmin_seq), n_sim_yrs))
  bmin_max_value <- 0
  c_slope_seq <- 0.01
  bmin_seq <- bmin_max_value
  tot_y <- tot_u <- prop_below <- TAC_zero <-
    matrix(0, nrow = length(c_slope_seq), ncol = length(bmin_seq))
  yield_array <- vB_fish_array <- 
    array(0, dim = c(length(c_slope_seq), length(bmin_seq), n_sim_yrs))
  #----------------------------------------------------------------------
  # run retrospective simulation for each cslope, bmin, draw, and sim yr
  #----------------------------------------------------------------------
  for (i in seq_along(c_slope_seq)) {
    c_slope <- c_slope_seq[i]
    for (j in seq_along(bmin_seq)) {
      b_lrp <- bmin_seq[j]
      set.seed(83) # challenge each bmin, cslope combo with same set of rec seqs
      for (k in seq_len(n_draws)) {
        # pick a single draw
        sub_post <- subset(post, post$.draw == unique(post$.draw)[k])
        
        # set leading parameters from sampled draw
        CR <- sub_post$cr[1]
        #rec_b <- sub_post$br[1]
        Ro <- sub_post$Ro[1]
        sbo <- Ro * sbro
        vbo <- Ro * vbro
        rec_a <- CR / sbro
        # stock-recruitment b
        if (rec_ctl == "ricker") { # Ricker
          rec_b <- log(CR) / (Ro * sbro) #rec_a <-  / sbro # exp(sub_post$ln_ar[1])
        }
        if (rec_ctl == "bh") { # Beverton-Holt
          rec_b <- (CR - 1) / sbo
        }
        
        wt <- sub_post$w
        wt_bar <- mean(wt)
        wt <- wt - wt_bar # correct nonzero wt over initialization period
        wt <- rep(wt, n_repeats)
  
        # extract the initial age structure
        Ninit <- sub_post[
          which(sub_post$year == retro_initial_yr),
          which(colnames(sub_post) %in% ages)
        ] %>%
          slice() %>%
          unlist(., use.names = FALSE)
        
        Ninit_yr_2 <- sub_post[
          which(sub_post$year == retro_initial_yr + 1),
          which(colnames(sub_post) %in% ages)
        ] %>%
          slice() %>%
          unlist(., use.names = FALSE)
        
        # nta matrix
        nta <- matrix(NA, nrow = length(wt) + 2, ncol = length(ages))
        
        # SSB, Rpred, vulnerable biomass vectors
        SSB <- Rpred <- vB_fish <- vB_survey <- rep(0, length(wt))
        nta[1, ] <- Ninit # initialize from posterior for retro_initial_yr
        nta[2, ] <- Ninit_yr_2
        t_last_ass <- 1 - ass_int # when was last survey / assessment (initialize)
        
        # run age-structured model for sim years
        for (t in seq_len(n_sim_yrs)[-n_sim_yrs]) { # years 1 to (n_sim_year-1)
          SSB[t] <- sum(nta[t, ] * f_a * w_a)
          vB_fish[t] <- sum(nta[t, ] * v_fish * w_a)
          vB_survey[t] <- sum(nta[t, ] * v_survey * w_a)
          
          if (t - t_last_ass == ass_int) { # assess every ass_int yrs
            t_last_ass <- t
            # note -0.5*(0.1)^2 corrects exponential effect on mean observation:
            # q_survey no longer used
            vB_obs <- vB_fish[t] * exp(sd_survey * (rnorm(1)) - 0.5 * (sd_survey)^2)
            TAC <- c_slope * (vB_obs - b_lrp)
            if (TAC < 0) {
              TAC <- 0
            }
            Ut <- ifelse((TAC / vB_fish[t]) < Ut_overall, (TAC / vB_fish[t]), Ut_overall)
          }
          rett <- ifelse(Ut / Ut_overall <= 1.0, Ut / Ut_overall, 1.0) #cap rett annual retention proportion at 1.0
          
          if(any(rett*ret_a > 1)){
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
          for (a in seq_len(n_ages)[-n_ages]) { # ages 2-19
            nta[t + 1, a + 1] <- nta[t, a] * exp(-M_a[a]) *
              (1 - Ut_overall * v_fish[a] * (ret_a[a] * rett + (1 - ret_a[a] * rett) * d_mort))       
          }
          
          # set rec value for next t
          nta[t + 2, 1] <- Rpred[t]
          
          # record performance metrics
          yield <- Ut_overall * rett * vB_fish[t]
          yield_array[i, j, t] <- yield_array[i, j, t] + yield
          vB_fish_array[i, j, t] <- vB_fish_array[i, j, t] + vB_fish[t]
          tot_y[i, j] <- tot_y[i, j] + yield
          tot_u[i, j] <- tot_u[i, j] + yield^0.3
          prop_below[i, j] <- prop_below[i, j] + ifelse(SSB[t] < sbo_prop * sbo, 1, 0)
          TAC_zero[i, j] <- TAC_zero[i, j] + ifelse(rett == 1, 1, 0)
        }
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
    "post" = post, "leading_pars" = leading_pars 
  )
  # create name and save .rds files for each run
 return(hcr_sim_list)
}

n_draws <- 30
rec_var <- 1.5 # variability of recruitment seqs after first seq
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
d_mort <- 0.3 # discard mortality
ah_ret <- 5
sd_ret <- 1
ret_a <- 1 / (1 + exp(-(ages - ah_ret) / sd_ret)) # retention by age vector
Ut_overall <- 1.0 # max U that fishermen can exert
sd_survey <- 0.4 # survey observation error
q_survey <- 1.0 # Cahill et al. 2021 assumed q_survey = 1.0
sbo_prop <- 0.1 # performance measure value to see if SSB falls below sbo_prop*sbo

# put it all in a tagged list
hcr_pars <- list(
  "n_draws" = n_draws,
  "rec_var" = rec_var,
  "n_repeats" = n_repeats,
  "retro_initial_yr" = retro_initial_yr,
  "retro_terminal_yr" = retro_terminal_yr,
  "n_sim_yrs" = n_sim_yrs,
  "ages" = ages,
  "initial_yr" = initial_yr,
  "initial_yr_minus_one" = initial_yr_minus_one,
  "d_mort" = d_mort,
  "ret_a" = ret_a,
  "Ut_overall" = Ut_overall,
  "sd_survey" = sd_survey,
  "q_survey" = q_survey,
  "sbo_prop" = sbo_prop
)

run <- get_hcr(which_lake = "lake newell", ass_int = 1)

my_df <- tibble(MSY_vB = run$MSY_vB_fish,
                MSY_yields = run$MSY_yields, 
                year = 1:hcr_pars$n_sim_yrs)

run$tot_y


long_data <- 
  lake_data %>%
  as.data.frame.table(., responseName = "value", dnn = c("cslope", "bmin")) %>%
  rename(
    "cslope" = "Var1",
    "bmin" = "Var2"
  ) %>%
  mutate(
    cslope = as.numeric(as.character(cslope)),
    bmin = as.numeric(as.character(bmin)), 
    lake = gsub("_", " ", i)
  )
if(names(hara_list)[1]==i){
  my_data <- long_data
} else {
  my_data <- rbind(my_data, long_data)
}



p <-
  my_df %>%
  ggplot(aes(x = year, y = MSY_vB)) +
  geom_line(colour = "#80b1d3", size = 0.5, alpha = 1) +
  ylab("Biomass vulnerable to fishing (blue) or yield (orange) ") +
  xlab("Year") +
  ggsidekick::theme_sleek() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) + 
  geom_line(aes(x=year, y=MSY_yields), colour="darkorange2", size=0.5)
p

#ggsave("plots/newell2.pdf")

# ass_int <- 10
# t_last_ass <- 1 - ass_int
# 
# for (t in 1:100) { 
#  print(t)
#   if (t - t_last_ass == ass_int) { # assess every ass_int yrs
#     t_last_ass <- t
#     message("assessment happened")
#   }
# }
