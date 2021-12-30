#libraries
library(tidyverse)
library(rstan)

#run it on all cores
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)

#compile the model
m <- rstan::stan_model("src/BERTA_single_lake.stan", verbose = F)

#read in the data

