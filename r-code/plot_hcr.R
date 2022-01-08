#----------------------------------------------------------------------
# plots
#----------------------------------------------------------------------
tot_y <- run$tot_y
tot_u <- run$tot_u

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

