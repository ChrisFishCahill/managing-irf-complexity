#----------------------------------------------------------------------
# hcr plots
#----------------------------------------------------------------------

# packages
library(paletteer)
library(viridis)
library(purrr)
library(ggplot2)
library(dplyr)
library(stringr)
library(gridExtra)

#----------------------------------------------------------------------
# read in the hcr sim files
sim_files <- list.files("sims/", full.names = TRUE)
#sim_files <- sim_files[grep("bh_cr_6", sim_files)]
#sim_files <- sim_files[grep("ricker_cr_6", sim_files)]

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
# plot the recruitment anomaly w(t) sequences used in hcrs
p <-
  all_posts %>%
  ggplot(aes(x = year, y = w, colour = lake, group = .draw)) +
  geom_line(colour = "#80b1d3", size = 0.05, alpha = 0.35) +
  ylab(expression(ln(w[t]))) +
  xlab("Year") +
  facet_wrap(~lake) +
  ggsidekick::theme_sleek() +
  ggtitle("Recruitment anomaly sequences used for harvest control rule development") +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
p

ggsave("plots/hcr_wts.pdf",
  width = 8,
  height = 5
)


#----------------------------------------------------------------------
# plots
#----------------------------------------------------------------------
# make the yield isopleth plots
my_data <- NA
yield_list <- big_list %>%
  purrr::map(~ .x$tot_y) 

for(i in names(yield_list)){
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
  if(names(yield_list)[1]==i){
    my_data <- long_data
    } else {
      my_data <- rbind(my_data, long_data)
    }
}
# rm(my_data)

plot_list <- vector('list', length(unique(my_data$lake)))
for(i in unique(my_data$lake)){
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
  plot_list[[which(unique(my_data$lake) == i )]] <- p1 
}

bigp <- gridExtra::grid.arrange(grobs = plot_list, ncol=3)

ggsave("plots/yield_isos.pdf", bigp,
       width = 8,
       height = 5
)

#----------------------------------------------------------------------
# make the utility isopleth plots

my_data <- NA
hara_list <- big_list %>%
  purrr::map(~ .x$tot_u) 

for(i in names(hara_list)){
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
  if(names(hara_list)[1]==i){
    my_data <- long_data
  } else {
    my_data <- rbind(my_data, long_data)
  }
}
# rm(my_data)

# make the utility isopleth plots
plot_list <- vector('list', length(unique(my_data$lake)))
for(i in unique(my_data$lake)){
  highlight <- my_data %>%
    filter(lake == i) %>%
    filter(value == max(value))
  p1 <- my_data %>%
    filter(lake == i) %>%
    ggplot(aes(bmin, cslope, z = value)) +
    geom_contour_filled(bins = 15) +
    ggsidekick::theme_sleek() +
    labs(fill = "Utility") +
    geom_point(data = highlight, aes(x = bmin, y = cslope), size = 1.75, show.legend = F) +
    scale_color_manual(values = c(NA, "black")) +
    xlab("Limit reference biomass (kg)") + 
    ggtitle(i) +
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
  plot_list[[which(unique(my_data$lake) == i )]] <- p1 
}

bigp <- gridExtra::grid.arrange(grobs = plot_list, ncol=3)

ggsave("plots/hara_isos.pdf", bigp,
       width = 8,
       height = 5
)

#----------------------------------------------------------------------
# 
test <- my_data %>%
  filter(lake == "lake newell")

summary(test$value)

p2 <- my_data %>%
  filter(lake == "lake newell") %>%
  ggplot(aes(bmin, cslope, z = value)) +
  geom_contour_filled(bins = 15) +
  ggsidekick::theme_sleek() +
  labs(fill = "Utility") +
  geom_point(data = highlight, aes(x = bmin, y = cslope), size = 1.75, show.legend = F) +
  scale_color_manual(values = c(NA, "black")) +
  xlab("Limit reference biomass (kg)") + 
  ggtitle(i) +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p2






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

# make a recruitment series plot
# k = post$.draw[1] # pick a draw
# sub_post <- subset(post, sampled_post$.draw == k)
# wt <- sub_post$w
#
# rec_var <- 1.0 # 1.2 might be fun to try
# wt <- rep(wt, n_repeats)
# df <- data.frame(wt = wt, sim_yrs = 1:n_sim_yrs)
#
# df %>%
#   ggplot(aes(x=sim_yrs, y=wt)) +
#   geom_point() +
#   geom_line() +
#   xlab("Year of Simulation") +
#   ylab("Recruitment Anomaly ln(wt)") +
#   ggsidekick::theme_sleek()
#----------------------------------------------------------------------
