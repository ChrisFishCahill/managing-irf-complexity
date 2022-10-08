vbk=.23
s=.86
age=c(1:20)
n=array(0.,20)  #numbers at age
n[1]=2.
for (a in 2:20){n[a]=n[a-1]*s}  #initial numbers at age
n[20]=n[20]/(1-s)  #set n in age 20 as plus group
w=array(0,20) #weight at age
vul=array(0,20)
w=(1-exp(-vbk*age))^3
mwt=array(0,20)
vul=1/(1+exp(-.5*(age-5)))
mwt=w/(1+exp(-.5*(age-7)))  #maturity x weight for SSB calculation
spro=sum(n*mwt)
cr=6
reca=cr/spro
ro=1
recb=(cr-1)/(ro*spro)
recmult=array(1,200)
by=c(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190)
for (t in by){recmult[t]=20}  #big recruitment every 10 yrs
ut=array(.1,200)
yield=array(0,200)
abar=array(0,200)
for (t in 1:200){
  yield[t]=ut[t]*sum(vul*n*w)
  ssb=sum(n*mwt)
  abar[t]=sum(age*n)/sum(n)
  n=s*n*(1-vul*ut[t])

  n[20]=n[20]+n[19]   #update plus group n
  for (a in 19:2){n[a]=n[a-1]}  #move fish up one age
  n[1]=reca*ssb/(1+recb*ssb)*recmult[t]   #put in new recruits
 # if (t/10-int(t/10)=0){n[1]=n[1]*20} #big cohort every 10 years
}
totalyield=sum(yield)
totalutility=sum(yield^.6)
plot(yield,type="l")
plot(abar,type="l")

# -----------------------------------------------------------
# Libraries and other stuff you'll need
# -----------------------------------------------------------
library(devtools)
library(TMB)
devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

# Set starting values:
#leading parameters/values for simulation
years = 1:200
ages = 1:20 # slot 1 = recruits 
cr = 6
vbk=.23
s=.86
recmult = rep(1.0, length(years))
by = seq(from = 10., to = max(years) - 10, by = 10)
recmult[by] = 20 # set by years to recmult
obj_ctl = 1 

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

tmb_data$obj_ctl <- 0
obj <- MakeADFun(tmb_data,
                 tmb_pars,
                 DLL = "om"
)

obj$report()$`yield`
yield

obj$fn(obj$par)
obj$gr(obj$par)

opt <- nlminb(obj$par, obj$fn, obj$gr,
       control = list(eval.max = 1000, iter.max = 500),
       lower = rep(0, length(years)), upper = rep(1, length(years))
)

plot(opt$par, type = "l")
opt_cccr <- TMBhelper::fit_tmb(obj = obj, getsd = T, newtonsteps = 1)
opt_cccr




###############################################################
###############################################################
###############################################################
# declare the tagged data list for stan
library(rstan)

# compile the stan model 
m = rstan::stan_model("om-sims/src/om-age.stan", verbose = F)

#leading parameters/values for simulation
years = 1:200
ages = 1:20 # slot 1 = recruits 
cr = 6

recmult = rep(1, length(years))
by = seq(from = 10, to = max(years) - 10, by = 10)
recmult[by] = 20 # set by years to recmult
obj_ctl = 1 

stan_data <- list(
  n_year = length(years),
  n_age = length(ages),
  vbk = vbk,
  s = s,
  cr = cr, 
  recmult = recmult, 
  obj_ctl = obj_ctl
)

inits <- function() {
  list(
    Ut = rep(0.1, length(years))
  )
}

stan_data$obj_ctl <- 1
#run the model
fit <-
  rstan::optimizing(
    m,
    data = stan_data,
    init = inits
  )
Ut_om = fit$par[grep("Ut", names(fit$par))]

inits <- function() {
  list(
    Ut = Ut_om
  )
}

fit <-
  rstan::optimizing(
    m,
    data = stan_data,
    init = inits
  )

# useful for debugging:
# m = rstan::stan_model("om-sims/src/om-age.stan", verbose = F)
# 
# #run the model
fit <-
  rstan::sampling(
    m,
    data = stan_data,
    init = inits,
    pars =
      c(
        "Ut", "ssb", "abar", "yield"
      ),
     iter = 1000,
    warmup = 500,
    chains = 1
  )

library(ggqfc)
library(tidyverse)
library(tidybayes)

fit %>%
  spread_draws(Ut[year]) %>%
  mutate(
    value = Ut,
    year = year
  ) %>%
  summarise(
    lwr = quantile(Ut, 0.1),
    med = quantile(Ut, 0.5),
    upr = quantile(Ut, 0.9),
    lwr2 = quantile(Ut, 0.25),
    upr2 = quantile(Ut, 0.75),
  ) %>%
  ggplot(aes(x = year, y = med)) + 
  geom_point() + 
  geom_line() + 
  theme_qfc()

fit %>%
  spread_draws(abar[year]) %>%
  mutate(
    value = abar,
    year = year
  ) %>%
  summarise(
    lwr = quantile(abar, 0.1),
    med = quantile(abar, 0.5),
    upr = quantile(abar, 0.9),
    lwr2 = quantile(abar, 0.25),
    upr2 = quantile(abar, 0.75),
  ) %>%
  ggplot(aes(x = year, y = med)) + 
  geom_point() + 
  geom_line() + 
  theme_qfc()

fit %>%
  spread_draws(yield[year]) %>%
  mutate(
    value = yield,
    year = year
  ) %>%
  summarise(
    lwr = quantile(yield, 0.1),
    med = quantile(yield, 0.5),
    upr = quantile(yield, 0.9),
    lwr2 = quantile(yield, 0.25),
    upr2 = quantile(yield, 0.75),
  ) %>%
  ggplot(aes(x = year, y = med)) + 
  geom_point() + 
  geom_line() + 
  theme_qfc()

fit %>%
  spread_draws(ssb[year]) 



Ut_om = fit$par[grep("Ut", names(fit$par))]
plot(Ut_om ~ years, type = "l", ylim=c(0,1))
B <- fit$par[grep("B", names(fit$par))]
lines(B, col = "blue")
lines(R, col = "orange")

abar_om = fit$par[grep("abar", names(fit$par))]
plot(abar_om ~ years, type = "l")

ssb_om = fit$par[grep("ssb", names(fit$par))]
plot(ssb_om ~ years, type = "l")

yield_om = fit$par[grep("yield", names(fit$par))]
plot(yield_om ~ years, type = "l")
