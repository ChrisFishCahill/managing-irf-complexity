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

# testing recmult
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
# Set starting values
# leading parameters/values for simulation
years <- 1:2000
n_year <- length(years)
pbig <- 0.05
Rbig <- 10
sdr <- 0.4

ages <- 1:20 # slot 1 = recruits
cr <- 6
vbk <- .23
s <- .86
rinit <- 0.001
ro <- 1
uo <- 0.13
asl <- 0.5
ahm <- 6
upow <- 0.6
xinc <- 0.1

# simulate the om across these quantities
pbig <- 0.1
Rbig <- 0.99
sdr <- 0.3
ahv <- 5

# draw recruitment sequence
set.seed(1)
sim <- get_recmult(pbig, Rbig, sdr) 

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
  obj_ctl = 0, # 0 = MAY, 1 = utility
  xinc = 0.1, 
  hcr = 0 # linear = 0; moxnes = 1
)

#par=c(0.01,0.01,0.01,0.03,0.05,0.07,0.08,0.1,0.11,0.12,0.12)
tmb_pars <- list(par = par)
tmb_pars <- list(par = c(0.5, 0.5))


# compile and load the cpp
cppfile <- "om-sims/src/om_interp.cpp"
compile(cppfile)

dyn.load(dynlib("om-sims/src/om_interp"))
obj <- MakeADFun(tmb_data, tmb_pars,  silent = F, DLL = "om_interp")
obj$fn()
obj$gr()
#obj$hessian <- TRUE
#obj$he()

# run om simulation
opt <- nlminb(obj$par, obj$fn, obj$gr,  eval.max = 1000, iter.max = 1000,
              lower = rep(0, length(tmb_pars$par)),
              upper = rep(1, length(tmb_pars$par)))

#opt2$par
opt$par
obj$fn(opt$par)

to i*Xinc for i=0 to 10 and Y=par(i).
Y = X = rep(NA, length(opt$par))
for(i in 0:10){
 Y[i] = opt$par[i]
 X[i] = i*tmb_data$xinc
}
plot(Y~X, type = "b")

tmb_pars = list(par = opt2$par)
obj <- MakeADFun(tmb_data, tmb_pars,  silent = F, DLL = "om_interp")

prof <- tmbprofile(obj, lincomb = c(1,1,1,1,1,1,-1,1, 1, 1, 1))
plot(prof)

tmb_pars <- list(par = jitter(c(0.5, 0.5), 0.1))

tmb_pars <- function() {
  list(
    par = jitter(c(0.5, 0.5), 0.25)
  )
}
tmb_pars()
obj <- MakeADFun(tmb_data, tmb_pars,  silent = F, DLL = "om_interp")

options(mc.cores = parallel::detectCores())
library(tmbstan)
fit <- tmbstan(obj, lower = rep(0, length(tmb_pars$par)),
                    upper = rep(1, length(tmb_pars$par)), 
        silent = F, chains = 7, iter = 10000)


post <- extract(fit)
max(post$lp__)

fit %>% 
  spread_draws(par[ctl]) %>%
  ggplot(aes(x = par, colour = ctl))+
  geom_histogram()

