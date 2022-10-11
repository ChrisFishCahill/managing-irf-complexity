# -----------------------------------------------------------
# Omniscient manager control aka open-loop optimization
# Carl Walters and Chris Cahill 9 Oct 2022
# -----------------------------------------------------------
library(devtools)
library(TMB)
devtools::install_github("ChrisFishCahill/gg-qfc")
library(ggqfc)
library(tidyverse)
library(ggtext)
library(cowplot)
library(ggpmisc)
library(future)
library(furrr)
# -----------------------------------------------------------
# function to get recmult sequences:

get_recmult <- function(pbig, Rbig, sdr) {
  urand <- runif(n_year, 0, 1)
  Nrand <- rnorm(n_year, 0, 1)
  recmult <- rep(1, n_year)
  rlow <- (1 - pbig * Rbig) / (1 - pbig)
  if (rlow < 0) rlow <- 0
  for (t in 1:n_year) {
    recmult[t] <- rlow
    if (urand[t] < pbig) {
      recmult[t] <- Rbig
    }
    recmult[t] <- recmult[t] * exp(sdr * Nrand[t])
  }
  out <- tibble(
    year = 1:n_year,
    urand, Nrand, recmult
  )
  list(dat = out)
}

# testing:
# years <- 1:200
# n_year <- length(years)
# pbig <- 0.05
# Rbig <- 10
# sdr <- 0.4
# set.seed(1)
# sim <- get_recmult(n_year, pbig = 0.02, Rbig, sdr)
# sim$dat %>%
#   ggplot(aes(x = year, y = recmult)) +
#   geom_line() +
#   theme_qfc()

# -----------------------------------------------------------

run_om <- function(pbig, Rbig, sdr, ahv, iter = NA) { # recruitment parameters
  # this if is used for parallel computations:
  if (!"om" %in% names(getLoadedDLLs())) {
    dyn.load(dynlib("om-sims/src/om"))
  }
  sim <- get_recmult(pbig, Rbig, sdr) # draw recruitment sequence
  tmb_data <- list(
    n_year = length(years),
    n_age = length(ages),
    vbk = vbk,
    s = s,
    cr = cr,
    rinit = rinit,
    ro = ro,
    uo = uo,
    asl = asl,
    ahv = ahv,
    ahm = ahm,
    upow = upow,
    ages = ages,
    recmult = sim$dat$recmult,
    obj_ctl = 0 # 0 = MAY, 1 = utility
  )
  tmb_pars <- list(Ut = rep(1, length(1:n_year))) # MAY 
  obj <- MakeADFun(tmb_data, tmb_pars, silent = T, DLL = "om")
  # run om simulation
  opt <- nlminb(obj$par, obj$fn, obj$gr,
    lower = rep(0, length(years)),
    upper = rep(1, length(years))
  )
  # re-run the optimization until convergence achieved
  while (opt$convergence == 1) {
    tmb_pars <- list(Ut = opt$par)
    obj <- MakeADFun(tmb_data, tmb_pars, silent = T, DLL = "om")
    opt <- nlminb(obj$par, obj$fn, obj$gr,
      lower = rep(0, length(years)),
      upper = rep(1, length(years))
    )
  }
  sim <- sim$dat %>% add_column(
    ssb = obj$report(opt$par)$`ssb`,
    abar = obj$report(opt$par)$`abar`,
    ut = opt$par,
    vulb = obj$report(opt$par)$`vulb`,
    objective = ifelse(tmb_data$obj_ctl == 0, "MAY", "utility"),
    convergence = opt$convergence,
    pbig, Rbig, sdr, iter = iter
  )
  yield <- sim
  # now do it for utility
  tmb_data$obj_ctl <- 1
  tmb_pars <- list(Ut = rep(0.1, length(1:n_year))) # utility 
  obj <- MakeADFun(tmb_data, tmb_pars, silent = T, DLL = "om")
  # run om simulation
  opt <- nlminb(obj$par, obj$fn, obj$gr,
    lower = rep(0, length(years)),
    upper = rep(1, length(years))
  )
  # re-run the optimization until convergence achieved
  while (opt$convergence == 1) {
    tmb_pars <- list(Ut = opt$par)
    obj <- MakeADFun(tmb_data, tmb_pars, silent = T, DLL = "om")
    opt <- nlminb(obj$par, obj$fn, obj$gr,
      lower = rep(0, length(years)),
      upper = rep(1, length(years))
    )
  }
  sim$ssb <- obj$report(opt$par)$`ssb`
  sim$abar <- obj$report(opt$par)$`abar`
  sim$ut <- opt$par
  sim$vulb <- obj$report(opt$par)$`vulb`
  sim$objective <- ifelse(tmb_data$obj_ctl == 0, "MAY", "utility")
  sim$convergence <- opt$convergence
  utility <- sim
  out <- bind_rows(yield, utility)
  out
}

# Set starting values:
# leading parameters/values for simulation
years <- 1:200
n_year <- length(years)
pbig <- 0.05
Rbig <- 10
sdr <- 0.4

ages <- 1:20 # slot 1 = recruits
cr <- 6
vbk <- .23
s <- .86
rinit <- 0.6
ro <- 1
uo <- 0.13
asl <- 0.5
ahv <- 5
ahm <- 6
upow <- 0.6

# compile and load the cpp
cppfile <- "om-sims/src/om.cpp"
compile(cppfile)
dyn.load(dynlib("om-sims/src/om"))

# simulate the om across these quantities
pbig <- seq(0.05)
Rbig <- c(20)
sdr <- c(0.4)
ahv <- c(5)
iter <- 1:12

to_sim <- expand.grid(pbig = pbig, Rbig = Rbig, sdr = sdr, ahv = ahv, iter = iter)
to_sim <- to_sim %>% distinct()
glimpse(to_sim)

# set.seed(1)
# system.time(
#  out <- furrr::pmap(to_sim, run_om) # testing
# )

future::plan(multisession)
system.time({
  out <- future_pmap(to_sim, run_om,
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )
})

data <-
  out %>%
  map_dfr(~.x)

# -----------------------------------------------------------
# now visualize solutions from the omniscient manager
# -----------------------------------------------------------

plot_dat_yield <- data %>% 
  filter(objective == "MAY") %>%
  select(year, ssb, abar, ut, iter) %>%
  pivot_longer(-c(year, iter))

#ot_dat_yield <- plot_dat_yield %>% pivot_longer(-year)

yield <- 
  plot_dat_yield %>%
  ggplot(aes(year, value, color = name)) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 1, lwd = 0.75, lty = 2) +
  scale_color_manual(
    name = NULL,
    values = c(ssb = "#0072B2", abar = "#009E73", ut = "#D55E00"),
    labels = c(
      ssb = "<i style='color:#0072B2'>SSB</i>",
      abar = "<i style='color:#009E73'>Abar</i>",
      ut = "<i style='color:#D55E00'>Ut</i>"
    ), 
  ) +
  labs(
    title = "Omniscient Manager
    <span style='font-size:11pt'>yield solutions for
    <span style='color:#0072B2;'>SSB</span>,
    <span style='color:#009E73;'>Abar</span>,
    <span style='font-size:11pt'> and
    <span style='color:#D55E00;'>Ut</span>
    </span>"
  ) +
  ylab("Value") +
  xlab("Year") +
  theme_qfc() +
  theme(
    plot.title = element_markdown(lineheight = 1.1, hjust = 0.5),
    legend.position = "non"
    # legend.text = element_markdown(size = 11)
  ) + facet_wrap(~as.factor(iter))
yield 
ggsave("plots/om-sims-yield.pdf", width = 11, height = 8)

plot_dat_utility <- data %>% 
  filter(objective == "utility") %>%
  select(year, ssb, abar, ut, iter) %>%
  pivot_longer(-c(year, iter))

utility <- 
  plot_dat_utility %>%
  ggplot(aes(year, value, color = name)) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 1, lwd = 0.75, lty = 2) +
  scale_color_manual(
    name = NULL,
    values = c(ssb = "#0072B2", abar = "#009E73", ut = "#D55E00"),
    labels = c(
      ssb = "<i style='color:#0072B2'>SSB</i>",
      abar = "<i style='color:#009E73'>Abar</i>",
      ut = "<i style='color:#D55E00'>Ut</i>"
    ), 
  ) +
  labs(
    title = "Omniscient Manager
    <span style='font-size:11pt'>yield solutions for
    <span style='color:#0072B2;'>SSB</span>,
    <span style='color:#009E73;'>Abar</span>,
    <span style='font-size:11pt'> and
    <span style='color:#D55E00;'>Ut</span>
    </span>"
  ) +
  ylab("Value") +
  xlab("Year") +
  theme_qfc() +
  theme(
    plot.title = element_markdown(lineheight = 1.1, hjust = 0.5),
    legend.position = "non"
    # legend.text = element_markdown(size = 11)
  ) + facet_wrap(~as.factor(iter))
utility 
ggsave("plots/om-sims-utility.pdf", width = 11, height = 8)

data %>%
  filter(year <= 175 & year >= 25) %>%
  ggplot(aes(x = vulb, y = ut, color = objective))+
  scale_color_manual(
    name = NULL,
    values = c(MAY = "#0072B2", utility = "#D55E00")
  ) +
  geom_point() + 
  facet_wrap(~iter) + 
  ylab(expression(Omniscient~manager~U[t])) + 
  xlab(expression(Vulnerable~biomass)) + 
  theme_qfc()
ggsave("plots/om-ut-vb.pdf", width = 11, height = 8)

data %>%
  ggplot(aes(x = year, y = recmult, group = iter))+
  geom_line() + 
  facet_wrap(~iter) + 
  ylab("Simulated Recruitment") + 
  xlab("Year") + 
  theme_qfc()
ggsave("plots/om-rec.pdf", width = 11, height = 8)



################################################################################
# Extra plotting code

# now do it for utility
plot_dat_utility <- data %>% 
  filter(objective == "utility") %>%
  select(year, ssb, abar, ut) %>%
  pivot_longer(-year)

hara <- ggplot(plot_dat_utility, aes(year, value, color = as.factor(name))) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 1, lwd = 0.75, lty = 2) +
  scale_color_manual(
    name = NULL,
    values = c(ssb = "#0072B2", abar = "#009E73", ut = "#D55E00"),
    labels = c(
      ssb = "<i style='color:#0072B2'>SSB</i>",
      abar = "<i style='color:#009E73'>Abar</i>",
      ut = "<i style='color:#D55E00'>Ut</i>"
    )
  ) +
  labs(
    title = "Omniscient Manager
    <span style='font-size:11pt'>utility solutions for
    <span style='color:#0072B2;'>SSB</span>,
    <span style='color:#009E73;'>Abar</span>,
    <span style='font-size:11pt'> and
    <span style='color:#D55E00;'>Ut</span>
    </span>"
  ) +
  ylab("Value") +
  xlab("Year") +
  theme_qfc() +
  theme(
    plot.title = element_markdown(lineheight = 1.1, hjust = 0.5),
    legend.position = "none"
    # legend.text = element_markdown(size = 11)
  )
hara

#-------------------------------------------------------------
# arrange plots into one big plot
#-------------------------------------------------------------
p1 <- plot_grid(yield, hara, ncol = 1)
p1

ggsave("plots/om-sims-tmb.pdf", width = 6, height = 5)

#-------------------------------------------------------------
# more yield plots
vulb <- obj$report(opt_yield$par)$`vulb`
plot_dat_yield <- plot_dat_yield %>% pivot_wider()
plot_dat_yield$vulb <- vulb

# remove last 25 rows to deal with risky omniscient manager:
n <- nrow(plot_dat_yield)
plot_dat_yield <- plot_dat_yield[1:(n - 25), ]
plot_dat_yield %>%
  ggplot(aes(x = abar, y = ut)) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
    after_stat(rr.label),
    sep = "*\", \"*"
  ))) +
  geom_point()

# plot_dat_yield[which(plot_dat_yield$ut > 0),] %>%
plot_dat_yield %>%
  ggplot(aes(x = vulb, y = ut)) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
    after_stat(rr.label),
    sep = "*\", \"*"
  ))) +
  geom_point()

# same plots for hara
#-------------------------------------------------------------
vulb <- obj$report(opt_hara$par)$`vulb`
plot_dat_hara <- plot_dat_hara %>% pivot_wider()
plot_dat_hara$vulb <- vulb

# remove last 25 rows to deal with risky omniscient manager:
n <- nrow(plot_dat_hara)
plot_dat_hara <- plot_dat_hara[1:(n - 25), ]
plot_dat_hara %>%
  ggplot(aes(x = abar, y = ut)) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
    after_stat(rr.label),
    sep = "*\", \"*"
  ))) +
  geom_point()

plot_dat_hara[which(plot_dat_hara$ut > 0), ] %>%
  ggplot(aes(x = vulb, y = ut)) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
    after_stat(rr.label),
    sep = "*\", \"*"
  ))) +
  geom_point()


#####################################################################
# vbk <- .23
# s <- .86
# age <- c(1:20)
# n <- array(0., 20) # numbers at age
# n[1] <- 2.
# for (a in 2:20) {
#   n[a] <- n[a - 1] * s
# } # initial numbers at age
# n[20] <- n[20] / (1 - s) # set n in age 20 as plus group
# w <- array(0, 20) # weight at age
# vul <- array(0, 20)
# w <- (1 - exp(-vbk * age))^3
# mwt <- array(0, 20)
# vul <- 1 / (1 + exp(-.5 * (age - 5)))
# mwt <- w / (1 + exp(-.5 * (age - 7))) # maturity x weight for SSB calculation
# spro <- sum(n * mwt)
# cr <- 6
# reca <- cr / spro
# ro <- 1
# recb <- (cr - 1) / (ro * spro)
# recmult <- array(1, 200)
# by <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190)
# for (t in by) {
#   recmult[t] <- 20
# } # big recruitment every 10 yrs
# ut <- array(.1, 200)
# yield <- array(0, 200)
# abar <- array(0, 200)
# for (t in 1:200) {
#   yield[t] <- ut[t] * sum(vul * n * w)
#   ssb <- sum(n * mwt)
#   abar[t] <- sum(age * n) / sum(n)
#   n <- s * n * (1 - vul * ut[t])
#
#   n[20] <- n[20] + n[19] # update plus group n
#   for (a in 19:2) {
#     n[a] <- n[a - 1]
#   } # move fish up one age
#   n[1] <- reca * ssb / (1 + recb * ssb) * recmult[t] # put in new recruits
#   # if (t/10-int(t/10)=0){n[1]=n[1]*20} #big cohort every 10 years
# }
# totalyield <- sum(yield)
# totalutility <- sum(yield^.6)
# plot(yield, type = "l")
# plot(abar, type = "l")
