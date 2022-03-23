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
# plot the recruitment anomaly w(t) sequences estimated by BERTA and used
# for hcr simulation

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

# ggsave("plots/hcr_wts.pdf",
#   width = 8,
#   height = 5
# )

my_data <- NA
wt_list <- big_list %>%
  purrr::map(~ .x$wt_seqs)

for (i in names(wt_list)) {
  lake_data <- wt_list[[i]]
  long_data <- setNames(melt(lake_data), c("sim_yr", "draw", "wt"))
  long_data$lake <- gsub("_", " ", i)
  if (names(wt_list)[1] == i) {
    my_data <- long_data
  } else {
    my_data <- rbind(my_data, long_data)
  }
}
my_data$sim_yr <- my_data$sim_yr + 1989

trace_data <- my_data %>%
  filter(draw == 1) %>%
  filter(sim_yr <= retro_terminal_yr)

trace_data2 <- my_data %>%
  filter(draw == 1)

p <-
  my_data %>%
  ggplot(aes(x = sim_yr, y = wt, colour = lake, group = draw)) +
  geom_line(color = "#80b1d3", size = 0.05, alpha = 0.35) +
  ylab(expression(ln(w[t]))) +
  xlab("Year") +
  facet_wrap(~lake) +
  ggsidekick::theme_sleek() +
  ggtitle("Recruitment anomaly sequences used for harvest control rule development") +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  geom_line(
    data = trace_data2, aes(x = sim_yr, y = wt), color = "black", size = 0.3,
    alpha = 0.4
  ) +
  geom_line(data = trace_data, aes(x = sim_yr, y = wt), color = "black", size = 0.3)
p

# ggsave("plots/hcr_wts.pdf",
#   width = 8,
#   height = 5
# )

# ggsave("plots/hcr_wts.jpeg",
#   width = 8,
#   height = 5
# )

my_data %>% filter(lake == "pigeon lake") %>%
  ggplot(aes(x = sim_yr, y = wt, colour = lake, group = draw)) +
  geom_line(color = "#80b1d3", size = 0.05, alpha = 0.35)

trace_data2 <- 
  trace_data2 %>%
  filter(lake == "baptiste lake")

trace_data <- 
  trace_data %>%
  filter(lake == "baptiste lake")
p <-
  my_data %>%
  filter(lake == "baptiste lake") %>%
  ggplot(aes(x = sim_yr, y = wt, colour = lake, group = draw)) +
  geom_line(color = "#80b1d3", size = 0.05, alpha = 0.35) +
  ylab(expression(ln(w[t]))) +
  xlab("Year") +
  ggsidekick::theme_sleek() +
  ggtitle("Recruitment anomaly sequences used for Baptiste Lake") + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  geom_line(
    data = trace_data2, aes(x = sim_yr, y = wt), color = "black", size = 0.75,
    alpha = 0.4
  ) +
  geom_line(data = trace_data, aes(x = sim_yr, y = wt), color = "black", size = 0.75)
p

ggsave("plots/hcr_wts_baptiste.jpeg",
       width = 8,
       height = 5
)


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
ggsave("plots/yield_isos_ass_int_1_dmort_0.3_sd_0.4.pdf", bigp,
       width = 8,
       height = 5
)

ggsave("plots/yield_isos_ass_int_1_dmort_0.3_sd_0.4.jpeg", bigp,
       width = 7.3,
       height = 4
)

#----------------------------------------------------------------------
# make the utility isopleth plots

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
# rm(my_data)

# make the utility isopleth plots
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
  plot_list[[which(unique(my_data$lake) == i)]] <- p1
}

bigp <- gridExtra::grid.arrange(grobs = plot_list, ncol = 3)
# bigp <- gridExtra::grid.arrange(grobs = plot_list, ncol=1)

# ggsave("plots/hara_isos.pdf", bigp,
#   width = 8,
#   height = 5
# )

ggsave("plots/hara_isos_ass_int_1_dmort_0.3_sd_0.4.jpeg", bigp,
       width = 7.3,
       height = 4
)


#----------------------------------------------------------------------
# plot vulnerable biomass vs. yield for msy policy

vb_list <- big_list %>%
  purrr::map(~ .x$MSY_vB_fish)

names(vb_list) <- gsub("_", " ", names(vb_list))

my_data <- tibble()

for (i in names(vb_list)) {
  lake_data <- vb_list[[i]]
  long_data <- tibble(
    MSY_vB = lake_data,
    year = retro_initial_yr:(retro_initial_yr + n_repeats * length(retro_initial_yr:retro_terminal_yr) - 1),
    lake = gsub("_", " ", i)
  )
  if (names(vb_list)[1] == i) {
    my_data <- long_data
  } else {
    my_data <- rbind(my_data, long_data)
  }
}

# now get MSY yields
yields_list <- big_list %>%
  purrr::map(~ .x$MSY_yields)

# names(yields_list) <- gsub("_", " ", names(yields_list))

my_data2 <- tibble()

for (i in names(yields_list)) {
  lake_data <- yields_list[[i]]
  long_data <- tibble(
    MSY_yield = lake_data,
    year = retro_initial_yr:(retro_initial_yr + 8 * length(retro_initial_yr:retro_terminal_yr) - 1),
    lake = gsub("_", " ", i)
  )
  if (names(vb_list)[1] == i) {
    my_data2 <- long_data
  } else {
    my_data2 <- rbind(my_data2, long_data)
  }
}

my_data$MSY_yields <- my_data2$MSY_yield
my_data[which(my_data$year == max(my_data$year)), c("MSY_vB", "MSY_yields")] <- NA

# p <-
#   my_data %>%
#   ggplot(aes(x = year, y = MSY_vB)) +
#   geom_line(colour = "#80b1d3", size = 0.5, alpha = 1) +
#   ylab("Biomass vulnerable to fishing (blue) or yield (orange) ") +
#   xlab("Year") +
#   facet_wrap(~lake, scales="free") +
#   ggsidekick::theme_sleek() +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "none",
#     plot.title = element_text(face = "bold", hjust = 0.5)
#   ) +
#   geom_line(aes(x=year, y=MSY_yields), colour="darkorange2", size=0.5)
# p
#
# ggsave("plots/vb_yield_msy_policy.pdf",
#        width = 8,
#        height = 5
# )
#
# p <-
#   my_data %>%
#   ggplot(aes(y = MSY_yields, x = MSY_vB)) +
#   geom_point() +
#   ylab("Optimum Yield") +
#   xlab("Vulnerable Biomass") +
#   facet_wrap(~lake, scales="free") +
#   ggsidekick::theme_sleek() +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "none",
#     plot.title = element_text(face = "bold", hjust = 0.5)
#   ) +
#   ggtitle("MSY policies") +
#   geom_abline() +
#   geom_abline(slope=0.3)
# p


#----------------------------------------------------------------------
# vulnerable biomass vs. yield for hara

vb_list <- big_list %>%
  purrr::map(~ .x$HARA_vB_fish)

names(vb_list) <- gsub("_", " ", names(vb_list))

my_data <- tibble()

for (i in names(vb_list)) {
  lake_data <- vb_list[[i]]
  long_data <- tibble(
    Utility_vB = lake_data,
    year = retro_initial_yr:(retro_initial_yr + 8 * length(retro_initial_yr:retro_terminal_yr) - 1),
    lake = gsub("_", " ", i)
  )
  if (names(hara_list)[1] == i) {
    my_data <- long_data
  } else {
    my_data <- rbind(my_data, long_data)
  }
}

# now get hara yields
yields_list <- big_list %>%
  purrr::map(~ .x$HARA_yields)

# names(yields_list) <- gsub("_", " ", names(yields_list))

my_data2 <- tibble()

for (i in names(yields_list)) {
  lake_data <- yields_list[[i]]
  long_data <- tibble(
    HARA_yield = lake_data,
    year = retro_initial_yr:(retro_initial_yr + 8 * length(retro_initial_yr:retro_terminal_yr) - 1),
    lake = gsub("_", " ", i)
  )
  if (names(hara_list)[1] == i) {
    my_data2 <- long_data
  } else {
    my_data2 <- rbind(my_data2, long_data)
  }
}

my_data$HARA_yield <- my_data2$HARA_yield
my_data[which(my_data$year == max(my_data$year)), c("Utility_vB", "HARA_yield")] <- NA

# p <-
#   my_data %>%
#   ggplot(aes(x = year, y = Utility_vB)) +
#   geom_line(colour = "#80b1d3", size = 0.5, alpha = 1) +
#   ylab("Biomass vulnerable to fishing (blue) or yield (orange) ") +
#   xlab("Year") +
#   facet_wrap(~lake, scales="free") +
#   ggsidekick::theme_sleek() +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "none",
#     plot.title = element_text(face = "bold", hjust = 0.5)
#   ) +
#   geom_line(aes(x=year, y=HARA_yield), colour="darkorange2", size=0.5)
# p

ggsave("plots/vb_yield_hara_policy.pdf",
       width = 8,
       height = 5
)

# p <-
#   my_data %>%
#   ggplot(aes(y = HARA_yield, x = Utility_vB)) +
#   geom_point() +
#   ylab("Optimum Yield") +
#   xlab("Vulnerable Biomass") +
#   facet_wrap(~lake, scales="free") +
#   ggsidekick::theme_sleek() +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "none",
#     plot.title = element_text(face = "bold", hjust = 0.5)
#   ) +
#   ggtitle("HARA policies") +
#   geom_abline() +
#   geom_abline(slope=0.3) +
#   geom_abline(slope=0.2)
# p
#
#----------------------------------------------------------------------
# Let's make a plot of the frontier of tradeoffs between utility and yield
# policies for simple linear hcrs
#----------------------------------------------------------------------

tot_y_list <- big_list %>%
  purrr::map(~ .x$tot_y)

tot_u_list <- big_list %>%
  purrr::map(~ .x$tot_u)

# tot_y <- tot_y_list[["baptiste_lake"]]
# tot_u <- tot_u_list[["baptiste_lake"]]

# get all the yield data
my_y_data <- tibble()
for (i in names(tot_y_list)) {
  lake_data <- tot_y_list[[i]]
  y_df <- lake_data %>%
    as.data.frame.table(., responseName = "x_yield", dnn = c("cslope", "bmin")) %>%
    rename(
      "cslope" = "Var1",
      "bmin" = "Var2",
    ) %>%
    mutate(
      cslope = as.numeric(as.character(cslope)),
      bmin = as.numeric(as.character(bmin)),
      lake = gsub("_", " ", i)
    )
  if (names(tot_y_list)[1] == i) {
    my_y_data <- y_df
  } else {
    my_y_data <- rbind(my_y_data, y_df)
  }
}

# get all the utility data
my_u_data <- tibble()
for (i in names(tot_u_list)) {
  lake_data <- tot_u_list[[i]]
  u_df <- lake_data %>%
    as.data.frame.table(., responseName = "y_utility", dnn = c("cslope", "bmin")) %>%
    rename(
      "cslope" = "Var1",
      "bmin" = "Var2",
    ) %>%
    mutate(
      cslope = as.numeric(as.character(cslope)),
      bmin = as.numeric(as.character(bmin)),
      lake = gsub("_", " ", i)
    )
  if (names(tot_u_list)[1] == i) {
    my_u_data <- u_df
  } else {
    my_u_data <- rbind(my_u_data, u_df)
  }
}

frontier <- left_join(my_y_data, my_u_data)
frontier$bmin <- round(frontier$bmin, 1)

plot_list <- vector("list", length(unique(frontier$lake)))
for (i in unique(frontier$lake)) {
  yint <- frontier %>%
    filter(lake == i) %>%
    summarize(value = max(y_utility) * 0.8)
  xint <- frontier %>% filter(lake == i) %>%
    summarize(value = max(x_yield) * 0.8)
  
  n_bins <- frontier %>%
    filter(
      lake == i,
      x_yield > 0,
      bmin <= 20,
    ) %>%
    summarize(n_bins = length(unique(bmin)))
  my_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(n_bins$n_bins)
  my_colors[1] <- "black"
  p1 <- frontier %>%
    filter(
      lake == i,
      x_yield > 0,
      bmin <= 20
    ) %>%
    ggplot(aes(x = x_yield, y = y_utility, z = bmin, colour = as.factor(bmin))) +
    ggsidekick::theme_sleek() +
    geom_rect(xmin=xint$value, xmax=Inf, ymin=yint$value, ymax=Inf,
              colour = "black", fill = "grey", alpha=0.05, linetype = "dashed", inherit.aes=FALSE) + 
    geom_point(size=0.75) +
    scale_color_manual(values = my_colors) +
    xlab("Yield") +
    ylab("Utility") +
    ggtitle(i) +
    labs(color = "blim") +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5), 
      legend.title.align=0.5
    ) +
    guides(colour=guide_legend(override.aes = list(size=1.5)))  
  plot_list[[which(unique(my_data$lake) == i)]] <- p1
}

bigp <- gridExtra::grid.arrange(grobs = plot_list, ncol = 3, nrow=2)

ggsave("plots/frontier.pdf", bigp,
       width = 16,
       height = 8.0
)

# ml <- marrangeGrob(plot_list, nrow=1, ncol=1, newpage = T, top=NULL)
# 
# ggsave("plots/frontier.pdf", ml,
#   width = 11,
#   height = 8.0
# )

#---------------------------------------------------------------------
# take the frontiers, make a plot of all HCRs achieving purdy good utility, yield
vB <- 0:50 
my_df <- NA
for (l in unique(frontier$lake)) {
  if(l == "baptiste lake"){wut_percent <- 0.98}
  if(l == "calling lake"){wut_percent <- 0.98}
  if(l == "lac ste. anne"){wut_percent <- 0.994}
  if(l == "lake newell"){wut_percent <- 0.991}
  if(l == "moose lake"){wut_percent <- 0.993}
  if(l == "pigeon lake"){wut_percent <- 0.94}
  # find maximum obtainable yield, utility:
  yint <- frontier %>%
    filter(lake == l) %>%
    summarize(value = max(y_utility) )
  xint <- frontier %>% filter(lake == l) %>%
    summarize(value = max(x_yield) )
  
  d <- frontier %>%
    filter(lake == l) %>%
    filter(x_yield >= wut_percent * xint$value & 
             y_utility >= wut_percent * yint$value
    )
  
  TAC <- matrix(NA, nrow = length(vB), ncol = nrow(d))
  
  for (i in 1:nrow(d)) {
    TAC[, i] <- d$cslope[i] * (vB - d$bmin[i])
  }
  TAC[which(TAC < 0)] <- NA
  TAC <- setNames(melt(TAC), c("vB", "group", "TAC"))
  TAC$vB <- rep(0:50, nrow(d))
  TAC$lake <- l
  if (frontier$lake[1] == l) {
    my_df <- TAC
  } else {
    my_df <- rbind(my_df, TAC)
  }
}

p <- my_df %>% ggplot(aes(x = vB, y = TAC, group = group)) + 
  xlab("Biomass Vulnerable to Fishing (kg/ha)") +
  ylab("TAC (kg/ha)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 50), 
                     breaks = c(0,5,10,15,20,30,40,50)) + 
  geom_line(alpha = 1, size = 1.0) + 
  ggsidekick::theme_sleek() +
  theme(axis.text=element_text(size=8), 
        panel.spacing = unit(1.1, "lines"), 
        plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  facet_wrap(~lake) + 
  ggtitle("Best linear HCRs that achieve both HARA utility and yield") 
p 

ggsave("plots/purdy_good_hcrs_best.pdf",
       width = 8,
       height = 5
)

#----------------------------------------------------------------------
# let's see if we can visualize hcrs across simulation runs :(
# read in the hcr sim files
sim_files <- list.files("sims/", full.names = TRUE)
sim_files <- sim_files[grep("bh_cr_6_hcr", sim_files)]
sim_files <- sim_files[grep("sd_0.05_d_mort_0.3", sim_files)]
sim_files

names <- str_extract(
  string = sim_files,
  pattern = "(?<=sims/).*(?=_bh|_ricker)"
)

ass_ints <- str_extract(
  string = sim_files,
  pattern = "(?<=hcr_).*(?=_sd)"
)

new_names <- paste(names, ass_ints, sep="_")

big_list <-
  sim_files %>%
  purrr::set_names(new_names) %>%
  purrr::map(readRDS)

#extract best yield cslope and bmin for each 
names(big_list)

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

res_table <- data.frame(bmin = rep(NA, length(unique(my_data$lake))), 
                        cslope = rep(NA, length(unique(my_data$lake))), 
                        lake = rep(NA, length(unique(my_data$lake))),
                        value = rep(NA, length(unique(my_data$lake))))

for (i in unique(my_data$lake)) {
  highlight <- my_data %>%
    filter(lake == i) %>%
    filter(value == max(value))
  counter <- which(unique(my_data$lake)==i) 
  res_table$bmin[counter] <- highlight$bmin
  res_table$cslope[counter] <- highlight$cslope
  res_table$lake[counter] <- highlight$lake
  res_table$value[counter] <- highlight$value
}
res_table$ass_ints <- ass_ints
res_table$lake <- names

ass_ints <- str_extract(
  string = res_table$ass_ints,
  pattern = "(?<=ass_int_).*"
)

res_table$ass_ints <- as.numeric(ass_ints)

#re-order based on ass_ints
res_table <- 
  res_table %>%
  group_by(lake) %>%
  arrange(ass_ints, .by_group=T)

res_table %>%
  ggplot(aes(x=ass_ints, y=value, colour=lake)) + 
  geom_point() + 
  geom_line() + 
  scale_x_continuous(breaks=c(1,3,5,10)) + 
  ggtitle("Yield vs. assessment interval, sd = 0.4, dmort = 0.3") + 
  ylab("Yield") + 
  xlab("Assessment interval (yrs)")

res_table$y_int <- - res_table$cslope*res_table$bmin

res_table$lake <- gsub("_", " ", res_table$lake)

my_colors <- colorRampPalette(brewer.pal(8, "Greys"))(length(unique(res_table$ass_ints)))
#my_colors <- c('#7b3294','#c2a5cf','#a6dba0','#008837')
p <- 
  res_table %>%
  ggplot() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 50), 
                     breaks = c(0,5,10,15,20,25,30,35,40,45, 50)) + 
  geom_abline(aes(intercept = y_int, slope = cslope,
                  colour=as.factor(ass_ints)), size=1.15) +
  ggsidekick::theme_sleek() + 
  xlab("Biomass Vulnerable to Fishing (kg/ha)") +
  ylab("TAC (kg/ha)") +
  theme(axis.text=element_text(size=8), 
        panel.spacing = unit(1.1, "lines"), 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        legend.title.align=0.5
  ) +
  facet_wrap(~lake) +
  labs(color = "Assessment \ninterval") +
  scale_colour_grey(start = 0) + 
  #scale_color_manual(values = my_colors) +
  ggtitle("Best linear HCRs for yield objective given discard mortality = 0.3 and survey sd = 0.05")
p

ggsave("plots/yields_sensitivity_sd_0.05_dmort_0.3.pdf", p,
       width = 8,
       height = 5
)

res_table$objective <- "yield"

#----------------------------------------------------------------------
# let's see if we can visualize hcrs across simulation runs for HARA
# read in the hcr sim files
sim_files <- list.files("sims/", full.names = TRUE)
sim_files <- sim_files[grep("bh_cr_6_hcr", sim_files)]
sim_files <- sim_files[grep("sd_0.05_d_mort_0.3", sim_files)]
sim_files

names <- str_extract(
  string = sim_files,
  pattern = "(?<=sims/).*(?=_bh|_ricker)"
)

ass_ints <- str_extract(
  string = sim_files,
  pattern = "(?<=hcr_).*(?=_sd)"
)

new_names <- paste(names, ass_ints, sep="_")

big_list <-
  sim_files %>%
  purrr::set_names(new_names) %>%
  purrr::map(readRDS)

#extract best yield cslope and bmin for each 
names(big_list)

my_data <- NA
yield_list <- big_list %>%
  purrr::map(~ .x$tot_u) #NOTE I HAVE CHANGED TOT_Y TO TOT_U

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

res_table_hara <- data.frame(bmin = rep(NA, length(unique(my_data$lake))), 
                             cslope = rep(NA, length(unique(my_data$lake))), 
                             lake = rep(NA, length(unique(my_data$lake))),
                             value = rep(NA, length(unique(my_data$lake))))

for (i in unique(my_data$lake)) {
  highlight <- my_data %>%
    filter(lake == i) %>%
    filter(value == max(value))
  counter <- which(unique(my_data$lake)==i) 
  res_table_hara$bmin[counter] <- highlight$bmin
  res_table_hara$cslope[counter] <- highlight$cslope
  res_table_hara$lake[counter] <- highlight$lake
  res_table_hara$value[counter] <- highlight$value
}
res_table_hara$ass_ints <- ass_ints
res_table_hara$lake <- names

ass_ints <- str_extract(
  string = res_table_hara$ass_ints,
  pattern = "(?<=ass_int_).*"
)

res_table_hara$ass_ints <- as.numeric(ass_ints)

#re-order based on ass_ints
res_table_hara <- 
  res_table_hara %>%
  group_by(lake) %>%
  arrange(ass_ints, .by_group=T)

res_table_hara %>%
  ggplot(aes(x=ass_ints, y=value, colour=lake)) + 
  geom_point() + 
  geom_line() + 
  scale_x_continuous(breaks=c(1,3,5,10)) + 
  ggtitle("Yield vs. assessment interval, sd = 0.4, dmort = 0.3") + 
  ylab("Yield") + 
  xlab("Assessment interval (yrs)")

res_table_hara$y_int <- - res_table_hara$cslope*res_table_hara$bmin

res_table_hara$lake <- gsub("_", " ", res_table_hara$lake)

res_table_hara$objective <- "hara"
res_table_hara

# put yield and hara dfs together
res_table  <- rbind(res_table, res_table_hara)
res_table$linetype <- ifelse(res_table$objective=="yield", "solid", "dashed")
res_table$objective <- ifelse(res_table$objective == "yield", "Yield", "HARA")
#my_colors <- colorRampPalette(brewer.pal(8, "Blues"))(length(unique(res_table$ass_ints)))
# my_colors <- c('#7b3294','#c2a5cf','#a6dba0','#008837')

#next bit necessary for stupid legend
GeomAbline$draw_key <- function(data, params, size) 
{
  segmentsGrob(0, 0.5, 1, 0.5, gp = gpar(col = alpha(data$colour, 
                                                     data$alpha), lwd = data$size * .pt, lty = data$linetype, 
                                         lineend = "butt"))
}  

p <- 
  res_table %>%
  ggplot() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 30)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 30), 
                     breaks = c(0,5,10,15,20,25,30)) + 
  geom_abline(aes(intercept = y_int, slope = cslope,
                  colour=as.factor(ass_ints), 
                  linetype= objective), 
              size=0.5
  ) +
  ggsidekick::theme_sleek() + 
  scale_colour_grey(start = 0.2) +
  xlab("Biomass Vulnerable to Fishing (kg/ha)") +
  ylab("Total Allowable Catch (kg/ha)") +
  theme(axis.text=element_text(size=8), 
        panel.spacing = unit(1.1, "lines"), 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        legend.title.align=0.5 
  ) +
  facet_wrap(~lake) +
  labs(color = "Assessment \ninterval", 
       linetype = "Management \n objective") +
  #scale_color_manual(values = my_colors) +
  scale_linetype_manual(values=c("dotted", "solid")) + 
  ggtitle("Best linear HCRs for Yield and HARA management objectives \ngiven discard mortality = 0.3 and survey sd = 0.05")
p

ggsave("plots/Yield_HARA_dmort_0.3_sd_0.05.pdf", p,
       width = 8,
       height = 5
)

# p1 <- my_data %>%
#   fil



# p1 <- my_data %>%
#   filter(lake == i) %>%
#   ggplot(aes(bmin, cslope, z = value)) +
#   geom_contour_filled(bins = 15) +
#   ggsidekick::theme_sleek() +
#   labs(fill = "Yield") +
#   geom_point(data = highlight, aes(x = bmin, y = cslope), size = 1.75, show.legend = F) +
#   scale_color_manual(values = c(NA, "black")) +
#   xlab("Limit reference biomass (kg)") +
#   ggtitle(i) +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5)
#   )
# # geom_vline(xintercept=10)
# plot_list[[which(unique(my_data$lake) == i)]] <- p1
#----------------------------------------------------------------------
# Extra plotting code below:
#----------------------------------------------------------------------
# some checks
# lsa <- frontier %>%
#   filter(lake == "pigeon lake",
#          bmin == 2.7)
# lsa %>%
# ggplot(aes(x = cslope, y = y_utility)) +
#   geom_line() +
#   ggsidekick::theme_sleek() +
#   geom_point() 
# lsa %>%
#   ggplot(aes(x = cslope, y = x_yield)) +
#   geom_line() +
#   ggsidekick::theme_sleek() +
#   geom_point()  
#   #geom_rug()
# 
# plot(lsa$x_yield)

# yields_list <- big_list %>%
#   purrr::map(~ .x$MSY_yields)
# 
# # names(yields_list) <- gsub("_", " ", names(yields_list))
# 
# my_data <- NA
# 
# for (i in names(yields_list)) {
#   lake_data <- yields_list[[i]]
#   long_data <- tibble(
#     MSY_yield = lake_data,
#     year = retro_initial_yr:(retro_initial_yr + 8 * length(retro_initial_yr:retro_terminal_yr) - 1),
#     lake = gsub("_", " ", i)
#   )
#   if (names(hara_list)[1] == i) {
#     my_data <- long_data
#   } else {
#     my_data <- rbind(my_data, long_data)
#   }
# }
# 
# p <-
#   my_data %>%
#   ggplot(aes(x = year, y = MSY_yield)) +
#   geom_line(colour = "#80b1d3", size = 0.5, alpha = 1) +
#   ylab("Maximum yields from MSY policy") +
#   xlab("Year") +
#   facet_wrap(~lake, scales = "free") +
#   ggsidekick::theme_sleek() +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "none",
#     plot.title = element_text(face = "bold", hjust = 0.5)
#   )
# p
# 
# ggsave("plots/yields.pdf",
#   width = 8,
#   height = 5
# )
# 
# # yield plot
# tot_y2 <-
#   as.data.frame.table(tot_y, responseName = "value", dnn = c("cslope", "bmin")) %>%
#   rename(
#     "cslope" = "Var1",
#     "bmin" = "Var2"
#   ) %>%
#   mutate(
#     cslope = as.numeric(as.character(cslope)),
#     bmin = as.numeric(as.character(bmin))
#   )
# 
# tot_u2 <-
#   as.data.frame.table(tot_u, responseName = "value", dnn = c("cslope", "bmin")) %>%
#   rename(
#     "cslope" = "Var1",
#     "bmin" = "Var2"
#   ) %>%
#   mutate(
#     cslope = as.numeric(as.character(cslope)),
#     bmin = as.numeric(as.character(bmin))
#   )
# highlight <- tot_u2 %>%
#   filter(value == max(value))
# 
# # utility plot
# p2 <- tot_u2 %>%
#   ggplot(aes(bmin, cslope, z = value)) +
#   geom_contour_filled(bins = 15) +
#   ggsidekick::theme_sleek() +
#   labs(fill = "Utility") +
#   geom_point(data = highlight, aes(x = bmin, y = cslope), size = 1.75, show.legend = F) +
#   scale_color_manual(values = c(NA, "black")) +
#   xlab("Limit reference biomass (kg)")
# p2
# 
# prop_below2 <-
#   as.data.frame.table(prop_below, responseName = "value", dnn = c("cslope", "bmin")) %>%
#   rename(
#     "cslope" = "Var1",
#     "bmin" = "Var2"
#   ) %>%
#   mutate(
#     cslope = as.numeric(as.character(cslope)),
#     bmin = as.numeric(as.character(bmin))
#   ) %>%
#   mutate(color = max(value) == value)
# 
# # proportion failing plot
# p3 <- prop_below2 %>%
#   ggplot(aes(bmin, cslope, z = value)) +
#   geom_contour_filled(bins = 15) +
#   ggsidekick::theme_sleek() +
#   labs(fill = "Proportion of years \nbelow 10% of average \nunfished SSB") +
#   scale_fill_viridis_d(direction = -1) +
#   xlab("Limit reference biomass (kg)")
# 
# p3
# 
# TAC_zero2 <-
#   as.data.frame.table(TAC_zero, responseName = "value", dnn = c("cslope", "bmin")) %>%
#   rename(
#     "cslope" = "Var1",
#     "bmin" = "Var2"
#   ) %>%
#   mutate(
#     cslope = as.numeric(as.character(cslope)),
#     bmin = as.numeric(as.character(bmin))
#   ) %>%
#   mutate(color = max(value) == value)
# 
# # zero catch
# p4 <- TAC_zero2 %>%
#   ggplot(aes(bmin, cslope, z = value)) +
#   geom_contour_filled(bins = 15) +
#   ggsidekick::theme_sleek() +
#   labs(fill = "Proportion of years \nwith no harvest") +
#   scale_fill_viridis_d(direction = -1) +
#   xlab("Limit reference biomass (kg)")
# 
# p4
# 
# # make the comparison plot for policies
# msys <- data.frame(
#   "yield" = MSY_yields,
#   "Policy" = "MSY policy",
#   "year" = retro_initial_yr:(retro_initial_yr + n_sim_yrs - 1)
# ) %>%
#   filter(year <= retro_terminal_yr)
# 
# haras <- data.frame(
#   "yield" = HARA_yields,
#   "Policy" = "HARA policy",
#   "year" = retro_initial_yr:(retro_initial_yr + n_sim_yrs - 1)
# ) %>%
#   filter(year <= retro_terminal_yr)
# 
# # obs_yields <- yields %>%
# #   select(year, med) %>%
# #   filter(year <= terminal_yr) %>%
# #   filter(year >= initialization_yr) %>%
# #   mutate(
# #     "yield" = med,
# #     "Policy" = "historical yield",
# #     "year" = year
# #   ) %>%
# #   select("yield", "Policy", "year")
# 
# all_yields <- rbind(msys, haras) # , obs_yields
# p5 <- all_yields %>%
#   ggplot(aes(x = year, y = yield, linetype = Policy, color = Policy)) +
#   geom_line(size = 1.5) +
#   scale_linetype_manual(values = c("dotted", "solid")) + # , "solid"
#   scale_color_manual(values = c("black", "grey")) + # , "black"
#   xlab("Year") +
#   ylab("Yield (kg)") +
#   ggsidekick::theme_sleek() +
#   guides(fill = guide_legend(title = ""))
# p5
# 
# 
# 
# # plot title for area
# which_x <- min(all_yields$year)
# hjust <- 0
# size <- 3
# 
# p5 <- p5 +
#   annotate("text", which_x, Inf,
#     vjust = 3, hjust = hjust,
#     label = which_lake, size = size
#   )
# 
# my_plot <- cowplot::plot_grid(p1, p2, p3, p4,
#   nrow = 2
# )
# 
# filename <- paste0("plots/", which_lake, "_cr6_hcr_plot.pdf")
# filename <- gsub(" ", "_", filename)
# my_tableau <- cowplot::plot_grid(p5, my_plot, nrow = 2, rel_heights = c(0.4, 0.6))
# ggsave(
#   filename = filename,
#   width = 10, height = 11, units = "in"
# )

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
# get sbos for carl 

# d <- all_posts %>%
#   group_by(lake) %>%
#   mutate( sbro = unique(sbro), 
#           Ro = Ro, 
#           sbo = sbro*Ro)
# sbos <- d %>% group_by(lake) %>%
#   summarize(sbos = unique(sbo)) 
# 
# sbos <- as.data.frame(sbos)
# sbo_mat <- matrix(NA, nrow= sum(sbos$lake=="pigeon lake"), ncol = length(unique(sbos$lake)))
# colnames(sbo_mat) <- unique(sbos$lake)
# 
# for(i in unique(sbos$lake)){
#   sub_dat <- subset(sbos, lake ==i)
#   sbo_mat[,i] = sub_dat$sbos
# }
# sbo_avg <- colMeans(sbo_mat)

# write.csv(sbo_mat, "data/sbo_all_draws.csv")
# write.csv(sbo_avg, "data/sbo_avg.csv")
#----------------------------------------------------------------------------------------------------
# compare rectilinear, simple line 

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

#-------------------------------------------------------------------------
# This section compares the precautionary and linear harvest control rules
#-------------------------------------------------------------------------
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

max_yields <- 
  my_data %>%
  group_by(lake) %>%
  filter(value == max(value)) %>%
  distinct()

# now, get the rectilinear rules
# read in the hcr sim files
sim_files <- list.files("sims/", full.names = TRUE)
sim_files <- sim_files[grep("bh_cr_6_hcr_ass_int_1_sd_0.4_d_mort_0.3_rule_precautionary", sim_files)]

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

yield_list <- big_list %>%
  purrr::map_dfr(~ .x$tot_y)
max_yields$precautionary_yields <- NA

# add precautionary rules to max_yields df
for (i in names(yield_list)) {
  lake_data <- yield_list[[i]]
  idx <- grep(gsub("_", " ", i), max_yields$lake)
  max_yields$precautionary_yields[idx] <- lake_data[1]
}

max_yields %>%
  ggplot(aes(x=precautionary_yields, y = value))+ 
  geom_point() + geom_abline() + 
  xlab("Precautionary HCR maximum yield") + 
  ylab("Simple HCR maximum yield") + 
  ggtitle("maximum yield objective")


#-------------------------------------------------------------------------
# This section compares the precautionary and linear harvest control rules for HARA
#-------------------------------------------------------------------------
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

my_data <- NA
yield_list <- big_list %>%
  purrr::map(~ .x$tot_u)

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

max_yields <- 
  my_data %>%
  group_by(lake) %>%
  filter(value == max(value)) %>%
  distinct()

# now, get the rectilinear rules
# read in the hcr sim files
sim_files <- list.files("sims/", full.names = TRUE)
sim_files <- sim_files[grep("bh_cr_6_hcr_ass_int_1_sd_0.4_d_mort_0.3_rule_precautionary", sim_files)]

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

yield_list <- big_list %>%
  purrr::map(~ .x$tot_u)
max_yields$precautionary_yields <- NA

# add precautionary rules to max_yields df
for (i in names(yield_list)) {
  lake_data <- yield_list[[i]]
  idx <- grep(gsub("_", " ", i), max_yields$lake)
  max_yields$precautionary_yields[idx] <- lake_data[1]
}

max_yields %>%
  ggplot(aes(x=precautionary_yields, y = value))+ 
  geom_point() + geom_abline() + 
  xlab("Precautionary HCR HARA yield") + 
  ylab("Simple HCR HARA yield") + 
  ggtitle("HARA objective") + 
  scale_y_continuous(limits=c(0,3)) + 
  scale_x_continuous(limits=c(0,3))


#######################################################################################################
#######################################################################################################
#######################################################################################################

#Deliverable plots

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
# plot the recruitment anomaly w(t) sequences estimated by BERTA and used
# for hcr simulation

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

# ggsave("plots/hcr_wts.pdf",
#   width = 8,
#   height = 5
# )

my_data <- NA
wt_list <- big_list %>%
  purrr::map(~ .x$wt_seqs)

for (i in names(wt_list)) {
  lake_data <- wt_list[[i]]
  long_data <- setNames(melt(lake_data), c("sim_yr", "draw", "wt"))
  long_data$lake <- gsub("_", " ", i)
  if (names(wt_list)[1] == i) {
    my_data <- long_data
  } else {
    my_data <- rbind(my_data, long_data)
  }
}
my_data$sim_yr <- my_data$sim_yr + 1989

trace_data <- my_data %>%
  filter(draw == 1) %>%
  filter(sim_yr <= retro_terminal_yr)

trace_data2 <- my_data %>%
  filter(draw == 1)

p <-
  my_data %>%
  ggplot(aes(x = sim_yr, y = wt, colour = lake, group = draw)) +
  geom_line(color = "#80b1d3", size = 0.05, alpha = 0.35) +
  ylab(expression(ln(w[t]))) +
  xlab("Year") +
  facet_wrap(~lake) +
  ggsidekick::theme_sleek() +
  ggtitle("Recruitment anomaly sequences used for harvest control rule development") +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  geom_line(
    data = trace_data2, aes(x = sim_yr, y = wt), color = "black", size = 0.3,
    alpha = 0.4
  ) +
  geom_line(data = trace_data, aes(x = sim_yr, y = wt), color = "black", size = 0.3)
p

ggsave("plots/hcr_wts.jpeg",
       width = 8,
       height = 5
)

trace_data2 <- 
  trace_data2 %>%
  filter(lake == "baptiste lake")

trace_data <- 
  trace_data %>%
  filter(lake == "baptiste lake")
p <-
  my_data %>%
  filter(lake == "baptiste lake") %>%
  ggplot(aes(x = sim_yr, y = wt, colour = lake, group = draw)) +
  geom_line(color = "#80b1d3", size = 0.05, alpha = 0.35) +
  ylab(expression(ln(w[t]))) +
  xlab("Year") +
  ggsidekick::theme_sleek() +
  ggtitle("Recruitment anomaly sequences used for Baptiste Lake") + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  geom_line(
    data = trace_data2, aes(x = sim_yr, y = wt), color = "black", size = 0.75,
    alpha = 0.4
  ) +
  geom_line(data = trace_data, aes(x = sim_yr, y = wt), color = "black", size = 0.75)
p

ggsave("plots/hcr_wts_baptiste.jpeg",
       width = 8,
       height = 5
)


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
ggsave("plots/yield_isos_ass_int_1_dmort_0.3_sd_0.4.pdf", bigp,
       width = 8,
       height = 5
)

ggsave("plots/yield_isos_ass_int_1_dmort_0.3_sd_0.4.jpeg", bigp,
       width = 7.3,
       height = 4
)

#----------------------------------------------------------------------
# make the utility isopleth plots

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
# rm(my_data)

# make the utility isopleth plots
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
  plot_list[[which(unique(my_data$lake) == i)]] <- p1
}

bigp <- gridExtra::grid.arrange(grobs = plot_list, ncol = 3)
# bigp <- gridExtra::grid.arrange(grobs = plot_list, ncol=1)

# ggsave("plots/hara_isos.pdf", bigp,
#   width = 8,
#   height = 5
# )

ggsave("plots/hara_isos_ass_int_1_dmort_0.3_sd_0.4.jpeg", bigp,
       width = 7.3,
       height = 4
)

#----------------------------------------------------------------------------------------
# let's see if we can visualize hcrs across simulation runs :(
# read in the hcr sim files
#-----------------------------------------------------------------------------------------
sim_files <- list.files("sims/", full.names = TRUE)
sim_files <- sim_files[grep("bh_cr_6_hcr", sim_files)]
sim_files <- sim_files[grep("sd_0.4_d_mort_0.3", sim_files)] 
sim_files <- sim_files[-grep("precautionary", sim_files)]
# sim_files <- sim_files[grep("bh_cr_6_hcr_ass_int_1_sd_0.4_d_mort_0.3", sim_files)]

sim_files

names <- str_extract(
  string = sim_files,
  pattern = "(?<=sims/).*(?=_bh|_ricker)"
)

ass_ints <- str_extract(
  string = sim_files,
  pattern = "(?<=hcr_).*(?=_sd)"
)

new_names <- paste(names, ass_ints, sep="_")

big_list <-
  sim_files %>%
  purrr::set_names(new_names) %>%
  purrr::map(readRDS)

#extract best yield cslope and bmin for each 
names(big_list)

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

res_table <- data.frame(bmin = rep(NA, length(unique(my_data$lake))), 
                        cslope = rep(NA, length(unique(my_data$lake))), 
                        lake = rep(NA, length(unique(my_data$lake))),
                        value = rep(NA, length(unique(my_data$lake))))

for (i in unique(my_data$lake)) {
  highlight <- my_data %>%
    filter(lake == i) %>%
    filter(value == max(value))
  counter <- which(unique(my_data$lake)==i) 
  res_table$bmin[counter] <- highlight$bmin
  res_table$cslope[counter] <- highlight$cslope
  res_table$lake[counter] <- highlight$lake
  res_table$value[counter] <- highlight$value
}
res_table$ass_ints <- ass_ints
res_table$lake <- names

ass_ints <- str_extract(
  string = res_table$ass_ints,
  pattern = "(?<=ass_int_).*"
)

res_table$ass_ints <- as.numeric(ass_ints)

#re-order based on ass_ints
res_table <- 
  res_table %>%
  group_by(lake) %>%
  arrange(ass_ints, .by_group=T)

p <- res_table %>%
  ggplot(aes(x=ass_ints, y=value, colour=lake)) + 
  geom_point() + 
  geom_line() + 
  scale_x_continuous(breaks=c(1,3,5,10)) + 
  ggtitle("Performance (MAY) vs. assessment interval, sd = 0.4, dmort = 0.3") + 
  ylab("Yield") + 
  xlab("Assessment interval (yrs)") + 
  ggsidekick::theme_sleek() 

ggsave("plots/performance_v_assint_sd_0.4_dmort_0.3.jpeg", p,
       width = 5,
       height = 3
)


res_table$y_int <- - res_table$cslope*res_table$bmin

res_table$lake <- gsub("_", " ", res_table$lake)

my_colors <- colorRampPalette(brewer.pal(8, "Greys"))(length(unique(res_table$ass_ints)))
#my_colors <- c('#7b3294','#c2a5cf','#a6dba0','#008837')
p <- 
  res_table %>%
  ggplot() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 50), 
                     breaks = c(0,5,10,15,20,25,30,35,40,45, 50)) + 
  geom_abline(aes(intercept = y_int, slope = cslope,
                  colour=as.factor(ass_ints)), size=1.15) +
  ggsidekick::theme_sleek() + 
  xlab("Biomass Vulnerable to Fishing (kg/ha)") +
  ylab("TAC (kg/ha)") +
  theme(axis.text=element_text(size=8), 
        panel.spacing = unit(1.1, "lines"), 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        legend.title.align=0.5
  ) +
  facet_wrap(~lake) +
  labs(color = "Assessment \ninterval") +
  scale_colour_grey(start = 0) + 
  #scale_color_manual(values = my_colors) +
  ggtitle("Best linear HCRs for MAY objective given discard mortality = 0.3 and survey sd = 0.4")
p

ggsave("plots/yield_hcr_sensitivity_sd_0.4_dmort_0.3.jpeg", p,
       width = 7.3,
       height = 4
)


#-----------------------------------------------------------------------------------------
#now do it for HARA

sim_files <- list.files("sims/", full.names = TRUE)
sim_files <- sim_files[grep("bh_cr_6_hcr", sim_files)]
sim_files <- sim_files[grep("sd_0.4_d_mort_0.3", sim_files)] 
sim_files <- sim_files[-grep("precautionary", sim_files)]
# sim_files <- sim_files[grep("bh_cr_6_hcr_ass_int_1_sd_0.4_d_mort_0.3", sim_files)]

sim_files

names <- str_extract(
  string = sim_files,
  pattern = "(?<=sims/).*(?=_bh|_ricker)"
)

ass_ints <- str_extract(
  string = sim_files,
  pattern = "(?<=hcr_).*(?=_sd)"
)

new_names <- paste(names, ass_ints, sep="_")

big_list <-
  sim_files %>%
  purrr::set_names(new_names) %>%
  purrr::map(readRDS)

#extract best yield cslope and bmin for each 
names(big_list)

my_data <- NA
yield_list <- big_list %>%
  purrr::map(~ .x$tot_u) #note not tot_y 

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

res_table <- data.frame(bmin = rep(NA, length(unique(my_data$lake))), 
                        cslope = rep(NA, length(unique(my_data$lake))), 
                        lake = rep(NA, length(unique(my_data$lake))),
                        value = rep(NA, length(unique(my_data$lake))))

for (i in unique(my_data$lake)) {
  highlight <- my_data %>%
    filter(lake == i) %>%
    filter(value == max(value))
  counter <- which(unique(my_data$lake)==i) 
  res_table$bmin[counter] <- highlight$bmin
  res_table$cslope[counter] <- highlight$cslope
  res_table$lake[counter] <- highlight$lake
  res_table$value[counter] <- highlight$value
}
res_table$ass_ints <- ass_ints
res_table$lake <- names

ass_ints <- str_extract(
  string = res_table$ass_ints,
  pattern = "(?<=ass_int_).*"
)

res_table$ass_ints <- as.numeric(ass_ints)

#re-order based on ass_ints
res_table <- 
  res_table %>%
  group_by(lake) %>%
  arrange(ass_ints, .by_group=T)

p <- res_table %>%
  ggplot(aes(x=ass_ints, y=value, colour=lake)) + 
  geom_point() + 
  geom_line() + 
  scale_x_continuous(breaks=c(1,3,5,10)) + 
  ggtitle("Performance (HARA) vs. assessment interval, sd = 0.4, dmort = 0.3") + 
  ylab("Utility") + 
  xlab("Assessment interval (yrs)") + 
  ggsidekick::theme_sleek() 

ggsave("plots/performance_hara_v_assint_sd_0.4_dmort_0.3.jpeg", p,
       width = 5,
       height = 3
)
res_table$y_int <- - res_table$cslope*res_table$bmin

res_table$lake <- gsub("_", " ", res_table$lake)

my_colors <- colorRampPalette(brewer.pal(8, "Greys"))(length(unique(res_table$ass_ints)))
#my_colors <- c('#7b3294','#c2a5cf','#a6dba0','#008837')
p <- 
  res_table %>%
  ggplot() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 50), 
                     breaks = c(0,5,10,15,20,25,30,35,40,45, 50)) + 
  geom_abline(aes(intercept = y_int, slope = cslope,
                  colour=as.factor(ass_ints)), size=1.15) +
  ggsidekick::theme_sleek() + 
  xlab("Biomass Vulnerable to Fishing (kg/ha)") +
  ylab("TAC (kg/ha)") +
  theme(axis.text=element_text(size=8), 
        panel.spacing = unit(1.1, "lines"), 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        legend.title.align=0.5
  ) +
  facet_wrap(~lake) +
  labs(color = "Assessment \ninterval") +
  scale_colour_grey(start = 0) + 
  #scale_color_manual(values = my_colors) +
  ggtitle("Best linear HCRs for HARA objective given discard mortality = 0.3 and survey sd = 0.4")
p

ggsave("plots/hara_hcr_sensitivity_sd_0.4_dmort_0.3.jpeg", p,
       width = 7.3,
       height = 4
)

