# -----------------------------------------------------------
# Omniscient manager control aka open-loop optimization
# Carl Walters and Chris Cahill 9 Oct 2022
# -----------------------------------------------------------

vbk <- .23
s <- .86
age <- c(1:20)
n <- array(0., 20) # numbers at age
n[1] <- 2.
for (a in 2:20) {
  n[a] <- n[a - 1] * s
} # initial numbers at age
n[20] <- n[20] / (1 - s) # set n in age 20 as plus group
w <- array(0, 20) # weight at age
vul <- array(0, 20)
w <- (1 - exp(-vbk * age))^3
mwt <- array(0, 20)
vul <- 1 / (1 + exp(-.5 * (age - 5)))
mwt <- w / (1 + exp(-.5 * (age - 7))) # maturity x weight for SSB calculation
spro <- sum(n * mwt)
cr <- 6
reca <- cr / spro
ro <- 1
recb <- (cr - 1) / (ro * spro)
recmult <- array(1, 200)
by <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190)
for (t in by) {
  recmult[t] <- 20
} # big recruitment every 10 yrs
ut <- array(.1, 200)
yield <- array(0, 200)
abar <- array(0, 200)
for (t in 1:200) {
  yield[t] <- ut[t] * sum(vul * n * w)
  ssb <- sum(n * mwt)
  abar[t] <- sum(age * n) / sum(n)
  n <- s * n * (1 - vul * ut[t])

  n[20] <- n[20] + n[19] # update plus group n
  for (a in 19:2) {
    n[a] <- n[a - 1]
  } # move fish up one age
  n[1] <- reca * ssb / (1 + recb * ssb) * recmult[t] # put in new recruits
  # if (t/10-int(t/10)=0){n[1]=n[1]*20} #big cohort every 10 years
}
totalyield <- sum(yield)
totalutility <- sum(yield^.6)
plot(yield, type = "l")
plot(abar, type = "l")

# -----------------------------------------------------------
# Libraries -- do the same analysis in TMB
# -----------------------------------------------------------
library(devtools)
library(TMB)
devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
# -----------------------------------------------------------

# Set starting values:
# leading parameters/values for simulation
years <- 1:200
ages <- 1:20 # slot 1 = recruits
cr <- 6
vbk <- .23
s <- .86
recmult <- rep(1.0, length(years))
by <- seq(from = 10., to = max(years) - 10, by = 10)
recmult[by] <- 20 # set by years to recmult
obj_ctl <- 0 # 1 = utility, 0 = MAY

cppfile <- "om-sims/src/om.cpp"
compile(cppfile)
dyn.load(dynlib("om-sims/src/om"))

tmb_data <- list(
  n_year = length(years),
  n_age = length(ages),
  vbk = vbk,
  s = s,
  cr = cr,
  ages = ages,
  recmult = recmult,
  obj_ctl = obj_ctl
)

tmb_pars <- list(
  Ut = rep(0.1, length(years))
)

obj <- MakeADFun(tmb_data,
  tmb_pars,
  DLL = "om"
)

obj$report()$`yield`[2]
# obj$fn(obj$par)
# obj$gr(obj$par)

opt_yield <- nlminb(obj$par, obj$fn, obj$gr,
  control = list(eval.max = 1000, iter.max = 500),
  lower = rep(0, length(years)), upper = rep(1, length(years))
)

tmb_data$obj_ctl <- 1
obj <- MakeADFun(tmb_data,
  tmb_pars,
  DLL = "om"
)

opt_hara <- nlminb(obj$par, obj$fn, obj$gr,
  control = list(eval.max = 1000, iter.max = 500),
  lower = rep(0, length(years)), upper = rep(1, length(years))
)

# -----------------------------------------------------------
# now visualize solutions from the omniscient manager
# -----------------------------------------------------------
library(ggqfc)
library(tidyverse)
library(tidybayes)
library(ggtext)
library(cowplot)

ssb <- obj$report(opt_yield$par)$`ssb`
abar <- obj$report(opt_yield$par)$`abar`
ut <- opt_yield$par
plot_dat <- data.frame(ssb, abar, ut, year = years)

plot_dat <- plot_dat %>% pivot_longer(-year)

yield <- ggplot(plot_dat, aes(year, value, color = as.factor(name))) +
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
    <span style='font-size:11pt'>yield solutions for
    <span style='color:#0072B2;'>SSB</span>,
    <span style='color:#009E73;'>Abar</span>,
    <span style='font-size:11pt'> and
    <span style='color:#D55E00;'>Ut</span>
    </span>"
  ) +
  ylab("Value") + xlab("Year") + 
  theme_qfc() +
  theme(
    plot.title = element_markdown(lineheight = 1.1, hjust = 0.5),
    legend.text = element_markdown(size = 11)
  )
yield

# now do it for utility
ssb <- obj$report(opt_hara$par)$`ssb`
abar <- obj$report(opt_hara$par)$`abar`
ut <- opt_hara$par
plot_dat <- data.frame(ssb, abar, ut, year = years)

plot_dat <- plot_dat %>% pivot_longer(-year)

hara <- ggplot(plot_dat, aes(year, value, color = as.factor(name))) +
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
  ylab("Value") + xlab("Year") + 
  theme_qfc() +
  theme(
    plot.title = element_markdown(lineheight = 1.1, hjust = 0.5),
    legend.text = element_markdown(size = 11)
  )
hara

#-------------------------------------------------------------
# arrange plots into one big plot
#-------------------------------------------------------------
p1 <- plot_grid(yield, hara, ncol = 1)
p1

ggsave("plots/om-sims-tmb.pdf", width = 5, height = 5)
