# make a multipanel trajectory plot for deliverable

# load packages
library(tidyverse)
library(tidybayes)
library(purrr)
# install.packages("devtools")
# devtools::install_github("seananderson/ggsidekick")

# declare some indeces 
t <- 2000 # first survey year
max_a <- max(ages)
rec_a <- min(ages)
initial_yr <- t - max_a + rec_a - 2 
initial_yr_minus_one <- initial_yr - 1

# find the fits corresponding to beverton-holt compensation ratio = 6
stan_files <- list.files("fits/", full.names = TRUE)
stan_files <- stan_files[grep("bh_cr_6", stan_files)]

# read fits into a big list in using map
big_list <-
  stan_files %>%
  purrr::set_names(.) %>%
  purrr::map(readRDS)

# extract r2, ssb, and ssb_c for pretty plotting 
r2 <- 
  big_list %>% map_dfr(function(big_list) { # extract recruits 
    big_list %>%
      spread_draws(R2[year]) %>%
      mutate(
        value = R2,
        year = year + initial_yr_minus_one
      ) %>%
      summarise(
        lwr = quantile(R2, 0.1),
        med = quantile(R2, 0.5),
        upr = quantile(R2, 0.9),
        lwr2 = quantile(R2, 0.25),
        upr2 = quantile(R2, 0.75),
      ) 
  }, .id = "stan_file") %>%
  mutate("name" = str_extract(string = stan_file, 
                              pattern = "(?<=fits/).*(?=_bh|ricker)")) %>%
  mutate(name = gsub("_", " ", name))

r2$which <- "recruits"

ssb <- big_list %>% map_dfr(function(big_list) { # extract estimated ssb
  big_list %>%
    spread_draws(SSB[year]) %>%
    mutate(
      value = SSB,
      year = year + initial_yr_minus_one
    ) %>%
    summarise(
      lwr = quantile(SSB, 0.1),
      med = quantile(SSB, 0.5),
      upr = quantile(SSB, 0.9),
      lwr2 = quantile(SSB, 0.25),
      upr2 = quantile(SSB, 0.75),
    )
}, .id = "stan_file") %>%
  mutate("name" = str_extract(string = stan_file, 
                              pattern = "(?<=fits/).*(?=_bh|ricker)")) %>%
  mutate(name = gsub("_", " ", name))

ssb$which <- "female biomass"

ssb_c <- big_list %>% map_dfr(function(big_list) { # extract "observed" ssb
  big_list %>%
    spread_draws(SSB_obs[counter]) %>%
    median_qi() %>%
    mutate(
      SSB_obs = SSB_obs,
      SSB_lower = SSB_obs - 0.30 * SSB_obs, # cv of 0.3 for error bars on plots
      SSB_upper = SSB_obs + 0.30 * SSB_obs
    )
}, .id = "stan_file") %>%
  mutate("name" = str_extract(string = stan_file, 
                              pattern = "(?<=fits/).*(?=_bh|ricker)")) %>%
  mutate(name = gsub("_", " ", name))

yr_dat <- data %>% 
  filter(name %in% contract_lakes) %>%
  summarise(name = name,
            year = year + t - 1) # t - 1 = 1999

ssb_c$year <- yr_dat$year

# join ssb, r2, ssb_c together into one large dataframe  
trajectory_dat <- rbind(ssb, r2)

trajectory_dat <- left_join(trajectory_dat, ssb_c,
                            by = c("name", "year")
)

# plot everything 
p <- trajectory_dat %>%
  ggplot(aes(x = year, y = med, colour = factor(which), fill = factor(which))) +
  geom_line(lwd = 0.5) +
  geom_ribbon(aes(ymin = lwr2, ymax = upr2),
              linetype = 0,
              alpha = 0.5
  ) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              linetype = 2, lwd = 0.5,
              alpha = 0.25
  ) +
  scale_fill_manual(values = c("steelblue", "darkorange2"), name = "") +
  scale_color_manual(values = c("steelblue", "darkorange2"), name = "") +
  xlab("Year") +
  ylab("Female biomass (kg/ha)") +
  xlim(1980, 2028) +
  scale_x_continuous(breaks = c(1980, 1990, 2000, 2010, 2020, 2028)) +
  ggsidekick::theme_sleek() +
  geom_point(aes(x = year, y = SSB_obs),
             size = 0.5, colour = "black",
             show.legend = FALSE
  ) +
  geom_linerange(aes(x = year, ymin = SSB_lower, ymax = SSB_upper),
                 colour = "black"
  ) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 1,
                                         name = "Age 2 recruits (N/ha)"
  )) +
  theme(
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(
      angle = 90, size = 8,
      hjust = 0.95, vjust = 0.2
    ),
    axis.text.y = element_text(size = 8)
  ) +
  facet_wrap(~name, ncol = 2, scales = "free")
p

ggsave("plots/trajectory_plot_bh6.pdf", width = 8,
       height = 5)
