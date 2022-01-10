#Misc run to get wt with Fearly = 0.5 for Carl
# declare HMC run parameters 
n_iter = 300
n_chains = 1
n_warmup = n_iter/2

contract_lakes <- c("lac ste. anne", "baptiste lake", 
                    "pigeon lake", "calling lake", 
                    "moose lake", "lake newell"
)

to_fit <- tibble(which_lake = contract_lakes)
to_fit$n_iter <- n_iter
to_fit$n_chains <- n_chains
to_fit$n_warmup <- n_warmup
to_fit$rec_ctl <- "bev-holt"
to_fit$cr_prior <- 6

# set up parallel processing plan 
future::plan(multisession)

# Run models and save fits -- 5.5 hours on my laptop
system.time({ 
  future_pwalk(to_fit, get_fit, 
               .options = furrr_options(seed = TRUE)
  ) 
})

# EXTRA stuff to get wt sequences to CJ

# extract most likely wt sequences from lakes bh cr =6
# for carl
# find the fits corresponding to beverton-holt compensation ratio = 6
#
retro_initial_yr <- 1990
retro_terminal_yr <- 2015

stan_files <- list.files("fits/", full.names = TRUE)
stan_files <- stan_files[grep("bh_cr_6", stan_files)]

big_list <-
  stan_files %>%
  purrr::set_names(.) %>%
  purrr::map(readRDS)

wt <- big_list %>%
  map_dfr(function(big_list) { # extract recruits
    big_list %>%
      spread_draws(w[year]) %>%
      mutate(
        value = w,
        year = year + initial_yr_minus_one
      ) %>%
      summarise(
        med = quantile(w, 0.5), # posterior median
      )
  }, .id = "stan_file") %>%
  mutate("name" = str_extract(
    string = stan_file,
    pattern = "(?<=fits/).*(?=_bh|ricker)"
  )) %>%
  mutate(name = gsub("_", " ", name)) %>%
  filter(year %in% retro_initial_yr:retro_terminal_yr)

wt %>% ggplot(aes(x=year, y=med, colour=name)) +
  geom_line() + 
  ylab("w(t)")

wts <- wt %>%
  pivot_wider(
    id_cols = -stan_file,
    names_from = name,
    values_from = med
  )


write.csv(wts, "data/most_likely_wts_cr_6_high_fearly.csv")

v <- big_list %>%
  map_dfr(function(big_list) { # extract recruits
    big_list %>%
      spread_draws(v[period]) %>%
      mutate(
        value = v,
      ) 
  }, .id = "stan_file") %>%
  mutate("name" = str_extract(
    string = stan_file,
    pattern = "(?<=fits/).*(?=_bh|ricker)"
  )) %>%
  mutate(name = gsub("_", " ", name)) %>%
  filter(period == 1)

v %>% ggplot(aes(x=value)) +
  geom_histogram() + 
  facet_wrap(~name) + 
  xlab(expression(F[early]))



str_extract(string = names(fits)[3], 
            pattern = "(?<=fits/).*(?=_bh|ricker)")) 

file_name <- paste0("sims/", file_name, "_hcr", "_ass_int_", ass_int, ".rds")

#------------------------------------

#debugging newell 
linf <- 74
lbar <- 57.57
M <- 0.2
theta <- 0.85
phi <- 2.02
vbk <- 0.2 
wl_beta <- 3.28
ages <- 2:20
n_ages <- length(ages)
v_a <- w_a <- l_a <- M_a <- Lo <- rep(NA, length(ages))

for(a in 1:n_ages){
  v_a[a] = ((linf/lbar)*(1 - exp(-vbk * ages[a])))^phi; 
  l_a[a] = (linf/lbar)*(1 - exp(-vbk * ages[a])); 
  M_a[a] = M/l_a[a]^theta;
  w_a[a] <- 0.00001*(linf*(1 - exp(-vbk * ages[a])))^wl_beta
  if(a == 1){ 
    Lo[a] = 1;
  } else{
    Lo[a] = Lo[a-1]*exp(-M_a[a-1]); 
  }
}

vbro <- sum(Lo * v_a * w_a)
vbro



data %>% group_by(name) %>% summarize(mean(beta_wl)) %>%
  print(n=Inf)
