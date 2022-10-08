#----------------------------------------------------------------------

library(reshape2)
library(rstan) 
library(tidyverse)

vbk <- 0.23
pbig <- 0.05 # Pr(spasmodic cohort)
surv <- 0.84
cost <- 0.05 # small cost to get solver to stop fishing after peak years

years <- 1:198
# create a vector indicating whether a big cohort happen
# carl does in excel using uniform deviate and algebra
# but I use rbinom()

set.seed(1)
R <- rbinom(n = length(years), size = 1, pbig)
obj_ctl <- 1 # utility = 1, yield = 0

m <- rstan::stan_model("om-sims/src/om.stan", verbose = F)

# declare the tagged data list for stan
stan_data <- list(
  n_year = length(years),
  vbk = vbk,
  surv = surv, 
  cost = cost,
  R = R, 
  obj_ctl = obj_ctl
)

inits <- function() {
  list(
    Ut = rep(0.1, length(years))
  )
}

#run the model
# fit <-
#   rstan::sampling(
#     m,
#     data = stan_data,
#     init = inits,
#     pars = 
#       c(
#         "Ut", "B", "Wdom", "Adom", "yield", 
#         "utility"
#       ),
#      iter = 2000,
#     # warmup = 500,
#     chains = 1
#   )

inits <- function() {
  list(
    Ut = rep(0.5, length(years))
  )
}

stan_data$obj_ctl <- 0 
#run the model
fit <-
  rstan::optimizing(
    m,
    data = stan_data,
    init = inits
  )

Ut_MLEs <- fit$par[grep("Ut", names(fit$par))]

Ut_MLEs <- Ut_MLEs[-length(Ut_MLEs)]
Ut_MLEs
plot(Ut_MLEs~years, type = "b")

inits <- function() {
  list(
    Ut = Ut_MLEs
  )
}

library(tidybayes)
library(ggqfc)

fit %>%
  spread_draws(B[year]) %>%
  mutate(
    value = B,
    year = year
  ) %>%
  summarise(
    lwr = quantile(B, 0.1),
    med = quantile(B, 0.5),
    upr = quantile(B, 0.9),
    lwr2 = quantile(B, 0.25),
    upr2 = quantile(B, 0.75),
  ) %>%
  ggplot(aes(x = year, y = med)) + 
  geom_point() + 
  geom_line() + 
  theme_qfc()


fit %>%
  spread_draws(Ut[year]) %>%
  mutate(
    value = Ut,
    year = year
  ) %>%
  summarise(
    lwr = quantile(Ut, 0.1),
    med = quantile(Ut, 0.5),
    upr = quantile(Ut, 0.9),
    lwr2 = quantile(Ut, 0.25),
    upr2 = quantile(Ut, 0.75),
  ) %>%
  ggplot(aes(x = year, y = med)) + 
  geom_point() + 
  geom_line() + 
  theme_qfc()




fit %>%
  spread_draws(obj) %>%
  mutate(
    value = obj
  ) %>%
  ggplot(aes(x = value))+
  geom_histogram() + 
  theme_qfc()





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