#----------------------------------------------------------------------
# hcr plots
#----------------------------------------------------------------------

library(paletteer)
library(viridis)
library(purrr)
library(ggplot2)
library(dplyr)
library(stringr)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(grid)

retro_initial_yr <- 1990 # initial year for retrospective analysis
retro_terminal_yr <- 2015
n_repeats <- 8

#----------------------------------------------------------------------
# read in the hcr sim files
sim_files <- list.files("sims/", full.names = TRUE)
sim_files <- sim_files[grep("bh_cr_6_hcr_ass_int_1_sd_0.4_d_mort_0.3", sim_files)]

# read in the data from all lakes used in Cahill et al. 2021
data <- readRDS("data/BERTA-wide-0-25.rds")

contract_lakes <- c(
  "lac ste. anne", "baptiste lake",
  "pigeon lake", "calling lake",
  "moose lake", "lake newell"
)

# sim_files <- sim_files[grep("ricker_cr_6", sim_files)]

names <- str_extract(
  string = sim_files,
  pattern = "(?<=sims/).*(?=_bh|_ricker)"
)

big_list <-
  sim_files %>%
  purrr::set_names(names) %>%
  purrr::map(readRDS)

all_posts <-
  big_list %>%
  purrr::map_dfr(~ .x$post, .id = "lake") %>%
  mutate(lake = gsub("_", " ", lake))

#----------------------------------------------------------------------
# make the yield isopleth plots

my_data <- NA
yield_list <- big_list %>%
  purrr::map(~ .x$tot_y)

# these loops are r-binding the total yield matrices for each file into one big dataframe
for (i in names(yield_list)) {
  lake_data <- yield_list[[i]]
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
  if (names(yield_list)[1] == i) {
    my_data <- long_data
  } else {
    my_data <- rbind(my_data, long_data)
  }
}
# rm(my_data)

plot_list <- vector("list", length(unique(my_data$lake)))
for (i in unique(my_data$lake)) {
  highlight <- my_data %>%
    filter(lake == i) %>%
    filter(value == max(value))
  p1 <- my_data %>%
    filter(lake == i) %>%
    ggplot(aes(bmin, cslope, z = value)) +
    geom_contour_filled(bins = 15) +
    ggsidekick::theme_sleek() +
    labs(fill = "Yield") +
    geom_point(data = highlight, aes(x = bmin, y = cslope), size = 1.75, show.legend = F) +
    scale_color_manual(values = c(NA, "black")) +
    xlab("Limit reference biomass (kg)") +
    ggtitle(i) +
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
  # geom_vline(xintercept=10)
  plot_list[[which(unique(my_data$lake) == i)]] <- p1
}

bigp <- gridExtra::grid.arrange(grobs = plot_list, ncol = 3)

# bigp <- gridExtra::grid.arrange(grobs = plot_list, ncol=1)
# ggsave("plots/yield_isos_ass_int_1_dmort_0.3_sd_0.4.pdf", bigp,
#        width = 8,
#        height = 5
# )
#----------------------------------------------------------------------
# calculate may_data for msy

may_data <-
  my_data %>%
  group_by(lake) %>%
  filter(value == max(value)) %>%
  mutate(
    MAY = value,
    BMAY = MAY / cslope + bmin,
    UMAY = MAY / BMAY
  ) %>%
  select(!value)

# may_data %>%
#  write.csv(.,file = "C:/Users/Chris Cahill/Documents/GitHub/managing-irf-complexity/data/may_yield_calcs.csv")

may_data


#----------------------------------------------------------------------
# calculate may_data for hara
my_data <- NA
hara_list <- big_list %>%
  purrr::map(~ .x$tot_u)

for (i in names(hara_list)) {
  lake_data <- hara_list[[i]]
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
  if (names(hara_list)[1] == i) {
    my_data <- long_data
  } else {
    my_data <- rbind(my_data, long_data)
  }
}


may_data <-
  my_data %>%
  group_by(lake) %>%
  filter(value == max(value)) %>%
  mutate(
    MAY = value,
    BMAY = MAY / cslope + bmin,
    UMAY = MAY / BMAY
  ) %>%
  select(!value)

# may_data %>%
#   write.csv(.,file = "C:/Users/Chris Cahill/Documents/GitHub/managing-irf-complexity/data/may_hara_calcs.csv")

may_data


#----------------------------------------------------------------
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
library(tidybayes)

may_data$Fmsy_equil <- may_data$msy_equil <- NA
for (i in unique(may_data$lake)) {
  lake_str <- gsub(" ", "_", i)
  fit_idx <- grep(lake_str, names(fits))
  fit <- fits[[fit_idx]]

  # extract estimated and derived parameters from BERTA
  devs <- fit %>%
    spread_draws(MSY, Fmsy) %>%
    median_qi()
  may_data$Fmsy_equil[which(may_data$lake == i)] <- devs$Fmsy
  may_data$msy_equil[which(may_data$lake == i)] <- devs$MSY
}

may_data %>%
  write.csv(., file = "C:/Users/Chris Cahill/Documents/GitHub/managing-irf-complexity/data/may_yield_calcs.csv")
