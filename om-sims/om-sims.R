#-----------------------------------------------------------------------
#          omniscient manager simulations  for spasmodic stocks
#                          from excel to R
#                   Walters & Cahill 7 Oct 2022
#-----------------------------------------------------------------------
vbk <- 0.23
pbig <- 0.05 # Pr(spasmodic cohort)
surv <- 0.84
cost <- 0.05 # small cost to get solver to stop fishing after peak years
seq <- 0 # --not used in our code, used to toggle between random number sequences in Excel

years <- 1:198
# create a vector indicating whether a big cohort happen
# carl does in excel using uniform deviate and algebra
# but I use rbinom()

set.seed(1)
R <- rbinom(n = length(years), size = 1, pbig)
Ut <- rep(0, length(years))

obj <- function(Ut) { # want optimizer to solve for Ut
  # bring in pars
  vbk <- leading_pars[1]
  pbig <- leading_pars[2]
  surv <- leading_pars[3]
  cost <- leading_pars[4]
  seq <- leading_pars[5]

  # set up age of dominant cohort vector
  Adom <- Wdom <- Ndom <- Nresid <- B <- rep(NA, length(years)) # dom = dominant cohort

  # initialize
  Adom[1] <- Ndom[1] <- 1
  Nresid[1] <- 0
  Wdom[1] <- (1 - exp(-vbk * Adom[1]))^3
  B[1] <- Ndom[1] * Wdom[1] + Nresid[1]

  for (t in 2:length(years)) {
    if (R[t - 1] == 0) {
      Adom[t] <- Adom[t - 1] + 1 # increment age by one year
      Ndom[t] <- Ndom[t - 1] * surv * (1 - Ut[t - 1])
      Nresid[t] <- Nresid[t - 1] * surv * (1 - Ut[t - 1]) # carry over residual (if present)
    } else {
      Adom[t] <- 1 # new cohort, start a counter over
      Ndom[t] <- R[t - 1]
      Nresid[t] <- Ndom[t - 1] * surv * (1 - Ut[t - 1])
    }
    Wdom[t] <- (1 - exp(-vbk * Adom[t]))^3
    B[t] <- Ndom[t] * Wdom[t] + Nresid[t]
  }
  if (which_obj == "catch") {
    obj <- sum(Ut * B) - cost * sum(Ut)
  }
  if (which_obj == "utility") {
    obj <- sum((Ut * B)^0.6) - cost * sum(Ut)
  }
  # my_data <- data.frame(years, Ut, Ndom, Nresid, Adom, Wdom, R, B)
  -obj
}
which_obj <- "utility"
leading_pars <- c(vbk, pbig, surv, cost, seq)

obj(Ut)

# I had to run a few times for utility, catch was easy
opt <- nlminb(Ut, obj,
              control = list(eval.max = 1000, iter.max = 500),
              lower = rep(0, length(Ut)), upper = rep(1, length(Ut)),
)
opt$objective

opt2 <- nlminb(opt$par, obj,
              control = list(eval.max = 1000, iter.max = 500),
              lower = rep(0, length(Ut)), upper = rep(1, length(Ut)),
)
opt2$objective

while(opt2$objective < opt$objective){
  opt$objective <- opt2$objective
  opt2 <- nlminb(opt2$par, obj,
                 control = list(eval.max = 1000, iter.max = 500),
                 lower = rep(0, length(Ut)), upper = rep(1, length(Ut)),
  )
  message(paste0("nlminb call completed"))
}

Ut <- opt2$par

plot(my_data$B ~ my_data$years, type = "l", col = "blue", ylim = c(0, 1), 
     lwd = 1.5, main="total biomass, Rt for Utility objective", ylab="", 
     xlab = "Year")
lines(my_data$R ~ my_data$years, type = "l", col = "black", lwd = 1.5)
lines(opt$par ~ my_data$years, type = "l", col = "orange", lwd = 1.5) 
legend(x = "topright",          # Position
       legend = c("biomass", "Recruits", "Ut"),  # Legend texts
       col = c("blue", "black", "orange"),       # Line colors
       lwd = 1.5)               

pdf("plots/om-utility.pdf",
 width = 12, height = 8
)
plot(my_data$Adom, opt$par, type="l", main = "Optimum U vs. Adom", 
     ylab="Optimum Ut", xlab = "Age of Dominant Cohort")

text(my_data$Adom, opt$par, labels = my_data$years)
dev.off()

