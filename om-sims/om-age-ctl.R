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

cppfile <- "om-sims/src/om_age.cpp"
compile(cppfile)
dyn.load(dynlib("om-sims/src/om_age"))

# read in carl's recmult
library(readxl)
Rmult = read_excel("C:/Users/Chris/Documents/manuscripts/alta harvest control rules/spasmodic age model optimization.xlsx","MaxU", range = "AG9:AG209")
years = 1:200
ages = 1:20
Ro = 1
cr = 6
su = 0.86
ahv = 5
ahm = 6
asl = 0.5
vbk = 0.5
uo = 0.13
Rinit = 0.6
obj_ctl = 0 # 0 = yield, 1 = hara

tmb_data <- list(
  n_year = length(years),
  n_age = length(ages),
  Ro = Ro, 
  cr = cr, 
  su = su, 
  ahv = ahv, 
  ahm = ahm, 
  asl = asl, 
  vbk = vbk,
  uo = uo, 
  Rinit = Rinit, 
  ages = ages,
  recmult = Rmult$Rmult,
  obj_ctl = obj_ctl
)

tmb_pars <- list(
  Ut = rep(0.1, length(years))
)

obj <- MakeADFun(tmb_data,
                 tmb_pars,
                 DLL = "om_age"
)
obj$report()$`sbro`

# sum(obj$report()$`yield` - yield)
# sum(obj$report()$`abar` - abar)
# obj$report()$`yield`[1] 

# obj$fn(obj$par)
# obj$gr(obj$par)