#----------------------------------------------------------------------
# Harvest control rules for Alberta Walleye lakes
# Cahill & Walters 3 Jan 2022
# TODO:
# 1) save output somehow
# 2) domestic netting
#----------------------------------------------------------------------

# load packages
library(tidyverse)
library(tidybayes)
library(purrr)

# install.packages("devtools")
# devtools::install_github("seananderson/ggsidekick")

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


get_hcr <- function(which_lake = "lac ste. anne", hcr_pars = hcr_pars) {
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
  ass_int <- hcr_pars$ass_int
  obs_sd <- hcr_pars$obs_sd
  q_survey <- hcr_pars$q_survey
  Ut_limit <- hcr_pars$Ut_limit
  sbo_prop <- hcr_pars$sbo_prop

  #--------------------------------------------------------------------
  # subset lake-specific posterior from all fits
  #--------------------------------------------------------------------
  fit_idx <- which(grepl(gsub(" ", "_", which_lake), names(fits)))
  fit <- fits[[fit_idx]]
  rec_ctl <- ifelse(grep("bh", names(fits)[fit_idx]), "bh", "ricker")

  # extract estimated and derived parameters from BERTA
  devs <- fit %>%
    spread_draws(Ro, ar, br) %>%
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
  c_slope_seq <- seq(from = 0.05, to = 1.0, by = 0.05)
  bmin_seq <- seq(from = 0, to = 1.0 * Ro_map * vbro, length.out = length(c_slope_seq))
  tot_y <- tot_u <-
    prop_below <- TAC_zero <- matrix(0, nrow = length(c_slope_seq), ncol = length(bmin_seq))
  yield_array <- array(0, dim = c(length(c_slope_seq), length(bmin_seq), n_sim_yrs))

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
        rec_a <- sub_post$ar[1]
        rec_b <- sub_post$br[1]
        Ro <- sub_post$Ro[1]
        sbo <- Ro * sbro
        vbo <- Ro * vbro

        wt <- sub_post$w
        wt_bar <- mean(wt)
        wt <- wt - wt_bar # correct nonzero wt over initialization period
        wt <- rep(wt, n_repeats)

        # multiply recruitment anomalies after initial period by rec_var (ugly)
        wt[(length(retro_initial_yr:retro_terminal_yr) + 1):length(wt)] <-
          wt[(length(retro_initial_yr:retro_terminal_yr) + 1):length(wt)] * rec_var

        # extract the initial age structure
        Ninit <- sub_post[
          which(sub_post$year == retro_initial_yr),
          which(colnames(sub_post) %in% ages)
        ] %>%
          slice() %>%
          unlist(., use.names = FALSE)

        # nta matrix
        nta <- matrix(NA, nrow = length(wt), ncol = length(ages))

        # SSB, Rpred, vulnerable biomass vectors
        SSB <- Rpred <- vB_fish <- vB_survey <- rep(0, length(wt))
        nta[1, ] <- Ninit # initialize from posterior for retro_initial_yr
        t_last_ass <- 1 - ass_int # when was last survey / assessment (initialize)

        # run age-structured model for sim years
        for (t in seq_len(n_sim_yrs)[-n_sim_yrs]) { # years 1 to (n_sim_year-1)
          SSB[t] <- sum(nta[t, ] * f_a * w_a)
          vB_fish[t] <- sum(nta[t, ] * v_fish * w_a)
          vB_survey[t] <- sum(nta[t, ] * v_survey * w_a)

          if (t - t_last_ass == ass_int) { # assess every ass_int yrs
            t_last_ass <- t
            # note -0.5*(0.1)^2 corrects exponential effect on mean observation:
            vB_obs <- q_survey * vB_survey[t] * exp(obs_sd * (rnorm(1)) - 0.5 * (obs_sd)^2)
            TAC <- c_slope * (vB_obs - b_lrp)
            if (TAC < 0) {
              TAC <- 0
            }
            Ut <- ifelse((TAC / vB_fish[t]) < Ut_limit, (TAC / vB_fish[t]), Ut_limit)
          }
          rett <- Ut / Ut_overall # rett = annual retention proportion

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
          nta[t + 1, 1] <- Rpred[t]

          # record performance metrics
          yield <- Ut_overall * rett * vB_fish[t]
          yield_array[i, j, t] <- yield_array[i, j, t] + yield
          tot_y[i, j] <- tot_y[i, j] + yield
          tot_u[i, j] <- tot_u[i, j] + yield^0.3
          prop_below[i, j] <- prop_below[i, j] + ifelse(SSB[t] < sbo_prop * sbo, 1, 0)
          TAC_zero[i, j] <- TAC_zero[i, j] + ifelse(TAC == 0, 1, 0)
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
  row_idx <- which(tot_y == max(tot_y), arr.ind = TRUE)[1]
  col_idx <- which(tot_y == max(tot_y), arr.ind = TRUE)[2]
  MSY_yields <- yield_array[row_idx, col_idx, ] # best MSY yields

  row_idx <- which(tot_u == max(tot_u), arr.ind = TRUE)[1]
  col_idx <- which(tot_u == max(tot_u), arr.ind = TRUE)[2]
  HARA_yields <- yield_array[row_idx, col_idx, ] # best HARA yields

  hcr_sim_list <- list(
    "tot_y" = tot_y, "tot_u" = tot_u,
    "prop_below" = prop_below, "TAC_zero" = TAC_zero,
    "yield_array" = yield_array, "MSY_yields" = MSY_yields,
    "HARA_yields" = HARA_yields
  )
  hcr_sim_list
}

#----------------------------------------------------------------------
# declare some values for the simulations
n_draws <- 1
rec_var <- 1.0 # variability of recruitment seqs after first seq
n_repeats <- 8 # recruitment repeats
retro_initial_yr <- 1990 # initial year for retrospective analysis
retro_terminal_yr <- 2015
n_sim_yrs <- length(retro_initial_yr:retro_terminal_yr) * n_repeats
ages <- 2:20
t <- 2000 # first survey year
max_a <- max(ages) # max age
rec_a <- min(ages) # age at recruitment
initial_yr <- t - max_a + rec_a - 2 # 1980 = BERTA initial yr
initial_yr_minus_one <- initial_yr - 1
d_mort <- 0.3 # discard mortality
ah_ret <- 5
sd_ret <- 1
ret_a <- 1 / (1 + exp(-(ages - ah_ret) / sd_ret)) # retention by age vector
Ut_overall <- 0.5
ass_int <- 3 # how often to assess / run FWIN
obs_sd <- 0.1
q_survey <- 1.0 # Cahill et al. 2021 assumed q_survey = 1.0
Ut_limit <- 0.9 # limit TAC mortality to < this value
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
  "ass_int" = ass_int,
  "obs_sd" = obs_sd,
  "q_survey" = q_survey,
  "Ut_limit" = Ut_limit,
  "sbo_prop" = sbo_prop
)

#----------------------------------------------------------------------
# read in the data and run 
#----------------------------------------------------------------------

# extract some saved .stan fit names
paths <- dir("fits/", pattern = "\\.rds$")
paths <- paste0(getwd(), "/fits/", paths)
paths <- paths[grep("bh_cr_6", paths)]
fits <- map(paths, readRDS) %>%
  set_names(paths)

run <- get_hcr(which_lake = "pigeon lake", hcr_pars = hcr_pars)


#----------------------------------------------------------------------
# plots
#----------------------------------------------------------------------

# yield plot
tot_y2 <-
  as.data.frame.table(tot_y, responseName = "value", dnn = c("cslope", "bmin")) %>%
  rename(
    "cslope" = "Var1",
    "bmin" = "Var2"
  ) %>%
  mutate(
    cslope = as.numeric(as.character(cslope)),
    bmin = as.numeric(as.character(bmin))
  )

highlight <- tot_y2 %>%
  filter(value == max(value))

p1 <- tot_y2 %>%
  ggplot(aes(bmin, cslope, z = value)) +
  geom_contour_filled(bins = 15) +
  ggsidekick::theme_sleek() +
  labs(fill = "Yield") +
  geom_point(data = highlight, aes(x = bmin, y = cslope), size = 1.75, show.legend = F) +
  scale_color_manual(values = c(NA, "black")) +
  xlab("Limit reference biomass (kg)")
p1

tot_u2 <-
  as.data.frame.table(tot_u, responseName = "value", dnn = c("cslope", "bmin")) %>%
  rename(
    "cslope" = "Var1",
    "bmin" = "Var2"
  ) %>%
  mutate(
    cslope = as.numeric(as.character(cslope)),
    bmin = as.numeric(as.character(bmin))
  )
highlight <- tot_u2 %>%
  filter(value == max(value))

# utility plot
p2 <- tot_u2 %>%
  ggplot(aes(bmin, cslope, z = value)) +
  geom_contour_filled(bins = 15) +
  ggsidekick::theme_sleek() +
  labs(fill = "Utility") +
  geom_point(data = highlight, aes(x = bmin, y = cslope), size = 1.75, show.legend = F) +
  scale_color_manual(values = c(NA, "black")) +
  xlab("Limit reference biomass (kg)")
p2

prop_below2 <-
  as.data.frame.table(prop_below, responseName = "value", dnn = c("cslope", "bmin")) %>%
  rename(
    "cslope" = "Var1",
    "bmin" = "Var2"
  ) %>%
  mutate(
    cslope = as.numeric(as.character(cslope)),
    bmin = as.numeric(as.character(bmin))
  ) %>%
  mutate(color = max(value) == value)

# proportion failing plot
p3 <- prop_below2 %>%
  ggplot(aes(bmin, cslope, z = value)) +
  geom_contour_filled(bins = 15) +
  ggsidekick::theme_sleek() +
  labs(fill = "Proportion of years \nbelow 10% of average \nunfished SSB") +
  scale_fill_viridis_d(direction = -1) +
  xlab("Limit reference biomass (kg)")

p3

TAC_zero2 <-
  as.data.frame.table(TAC_zero, responseName = "value", dnn = c("cslope", "bmin")) %>%
  rename(
    "cslope" = "Var1",
    "bmin" = "Var2"
  ) %>%
  mutate(
    cslope = as.numeric(as.character(cslope)),
    bmin = as.numeric(as.character(bmin))
  ) %>%
  mutate(color = max(value) == value)

# zero catch
p4 <- TAC_zero2 %>%
  ggplot(aes(bmin, cslope, z = value)) +
  geom_contour_filled(bins = 15) +
  ggsidekick::theme_sleek() +
  labs(fill = "Proportion of years \nwith no harvest") +
  scale_fill_viridis_d(direction = -1) +
  xlab("Limit reference biomass (kg)")

p4

# make the comparison plot for policies
msys <- data.frame(
  "yield" = MSY_yields,
  "Policy" = "MSY policy",
  "year" = retro_initial_yr:(retro_initial_yr + n_sim_yrs - 1)
) %>%
  filter(year <= retro_terminal_yr)

haras <- data.frame(
  "yield" = HARA_yields,
  "Policy" = "HARA policy",
  "year" = retro_initial_yr:(retro_initial_yr + n_sim_yrs - 1)
) %>%
  filter(year <= retro_terminal_yr)

# obs_yields <- yields %>%
#   select(year, med) %>%
#   filter(year <= terminal_yr) %>%
#   filter(year >= initialization_yr) %>%
#   mutate(
#     "yield" = med,
#     "Policy" = "historical yield",
#     "year" = year
#   ) %>%
#   select("yield", "Policy", "year")

all_yields <- rbind(msys, haras) # , obs_yields
p5 <- all_yields %>%
  ggplot(aes(x = year, y = yield, linetype = Policy, color = Policy)) +
  geom_line(size = 1.5) +
  scale_linetype_manual(values = c("dotted", "solid")) + # , "solid"
  scale_color_manual(values = c("black", "grey")) + # , "black"
  xlab("Year") +
  ylab("Yield (kg)") +
  ggsidekick::theme_sleek() +
  guides(fill = guide_legend(title = ""))
p5



# plot title for area
which_x <- min(all_yields$year)
hjust <- 0
size <- 3

p5 <- p5 +
  annotate("text", which_x, Inf,
    vjust = 3, hjust = hjust,
    label = which_lake, size = size
  )

my_plot <- cowplot::plot_grid(p1, p2, p3, p4,
  nrow = 2
)

filename <- paste0("plots/", which_lake, "_cr6_hcr_plot.pdf")
filename <- gsub(" ", "_", filename)
my_tableau <- cowplot::plot_grid(p5, my_plot, nrow = 2, rel_heights = c(0.4, 0.6))
ggsave(
  filename = filename,
  width = 10, height = 11, units = "in"
)


##############################################################################################
##############################################################################################
##############################################################################################

# EXTRA stuff to get wt sequences to CJ

# extract most likely wt sequences from lakes bh cr =6
# for carl
# find the fits corresponding to beverton-holt compensation ratio = 6
#
# retro_initial_yr <- 1990
# retro_terminal_yr <- 2015
#
# stan_files <- list.files("fits/", full.names = TRUE)
# stan_files <- stan_files[grep("bh_cr_12", stan_files)]
#
# big_list <-
#   stan_files %>%
#   purrr::set_names(.) %>%
#   purrr::map(readRDS)
#
# wt <- big_list %>%
#   map_dfr(function(big_list) { # extract recruits
#     big_list %>%
#       spread_draws(w[year]) %>%
#       mutate(
#         value = w,
#         year = year + initial_yr_minus_one
#       ) %>%
#       summarise(
#         med = quantile(w, 0.5), # posterior median
#       )
#   }, .id = "stan_file") %>%
#   mutate("name" = str_extract(
#     string = stan_file,
#     pattern = "(?<=fits/).*(?=_bh|ricker)"
#   )) %>%
#   mutate(name = gsub("_", " ", name)) %>%
#   filter(year %in% retro_initial_yr:retro_terminal_yr)
#
# wts <- wt %>%
#   pivot_wider(
#     id_cols = -stan_file,
#     names_from = name,
#     values_from = med
#   )

# write.csv(wts, "data/most_likely_wts_cr_12.csv")

##############################################################################################
#----------------------------------------------------------------------
# k = sampled_post$.draw[1] # pick a draw
# sub_post <- subset(sampled_post, sampled_post$.draw == k)
# wt <- sub_post$w
#
# rec_var <- 1.0 # 1.2 might be fun to try
# wt <- rep(wt, n_repeats)
# df <- data.frame(wt = wt, sim_yrs = 1:n_sim_yrs)
#
# df %>%
# ggplot(aes(x=sim_yrs, y=wt)) +
#   geom_point() +
#   geom_line() +
#   xlab("Year of Simulation") +
#   ylab("Recruitment Anomaly ln(wt)") +
#   ggsidekick::theme_sleek()
#----------------------------------------------------------------------
