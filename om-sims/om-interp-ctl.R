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
rinit <- 0.6
ro <- 1
uo <- 0.13
asl <- 0.5
ahm <- 6
upow <- 0.6
xinc <- 0.5

# simulate the om across these quantities
pbig <- 0.01
Rbig <- 1.2
sdr <- 0.3
ahv <- 6

# draw recruitment sequence
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
  xinc = 0.2
)

tmb_pars <- list(par = rep(0.5, 15))

# compile and load the cpp
cppfile <- "om-sims/src/om_interp.cpp"
compile(cppfile)
dyn.load(dynlib("om-sims/src/om_interp"))
obj <- MakeADFun(tmb_data, tmb_pars,  silent = F, DLL = "om_interp")
obj$fn()
obj$gr()
obj$hessian <- TRUE
obj$he()
# run om simulation
opt <- nlminb(obj$par, obj$fn, obj$gr, 
              lower = rep(0, length(years)),
              upper = rep(1, length(years)))

opt$par


obj$report(opt$par)$`ut`

fit2 <- TMBhelper::fit_tmb(obj = obj, getsd = T, newtonsteps = 0, 
                   lower = rep(0, length(tmb_pars$par)),
                   upper = rep(1, length(tmb_pars$par)))




#############################
# code the fucking scheme in R
xinc <- 2
vulb <- 7
idx <- as.integer(vulb/xinc)
idx_plus1 <- idx + 1
diff <- vulb/xinc -idx


# re-run the optimization until convergence achieved
while (opt$convergence == 1) {
  tmb_pars <- list(Ut = opt$par)
  obj <- MakeADFun(tmb_data, tmb_pars, silent = T, DLL = "om")
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                lower = rep(0, length(years)),
                upper = rep(1, length(years))
  )
}










to_sim <- expand.grid(pbig = pbig, Rbig = Rbig, sdr = sdr, ahv = ahv, iter = iter)
to_sim <- to_sim %>% distinct()
glimpse(to_sim)

# set.seed(1)
out <- purrr::pmap(to_sim, run_om) # testing


future::plan(multisession)
system.time({
  out <- future_pmap(to_sim, run_om,
                     .options = furrr_options(seed = TRUE),
                     .progress = TRUE
  )
})
