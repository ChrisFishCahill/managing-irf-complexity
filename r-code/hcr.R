# develop harvest control rules for Alberta Walleye lakes
# Cahill & Walters 2022 
# load packages

library(tidyverse)
library(tidybayes)
library(purrr)

# install.packages("devtools")
# devtools::install_github("seananderson/ggsidekick")

# list the fits
list.files("fits/", full.names = TRUE)

# pick a file 
which_file <- "fits/lac_ste._anne_bh_cr_6.rds"

# read fits into a big list in using map
fit <-
  which_file %>%
  purrr::set_names(.) %>%
  purrr::map(readRDS)

#----------------------------------------------------------------------
# ***N.B.***
# We need to preserve the historical frequency of weak and strong 
# year classes in our recruitment time series for our retrospective 
# analysis
# Thus, we will select 1990-2015 for these lakes for a 25 yr ref period
#----------------------------------------------------------------------
