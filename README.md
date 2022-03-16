## Developing harvest control rules for Alberta Walleye fisheries

This directory contains Stan and R code to fit age structured population dynamics models to Alberta FWIN data as per Cahill et al. 2021.  It also contains .R simulation code that develops harvest control rules for Special Harvest License (SHL) Walleye fisheries given these model fits. There have been three changes to this .stan model from what was published in Cahill et al. 2021:

* We now assume logistic vulnerability to fishing, which was not used in the paper.
* We have assumed a higher Fearly mortality rate prior of 0.5 (vs. 0.3 used in the paper). 
* We have corrected the estimates of MSY and Fmsy so that they now account for the lognormal bias. 

We made these changes to ensure the development of our harvest control rules was robust and more conservative.  None of the major findings from Cahill et al. 2021 changed when we make these additional changes. 

The directory is straightforward, with a "data" folder, an "r-code" folder, a "plots" folder, and an "src" (i.e., source) folder for .stan files. Tutorials are in the "Rmd" folder, and harvest control rule simulations are stored in the "sims" folder. The plots are mostly a work in progress, so be careful trying to interpret them (i.e., don't do this).  It was simply a location for me to share plots with Carl and others.  The scripts in the "r-code" folder that are potentially useful to folks are the `run.R` (calls and runs the stan file, saves fits), `hcr.R` (runs the harvest control rule simulator, saves simulations), and plotting files (calls either the Bayesian model fits or the simulation results). 

The project plan is as follows:

* Create a single lake version of the BERTA models from Cahill et al. 2021 so that one can easily run models for a single lake.
* Develop the necessary programs to evaluate harvest control rules for a collection of lakes given model fits from part 1. 
* Build capacity within AEP so that they can implement assessments and harvest control rule development on their own (via tutorials).

Additional notes:

* There are two .stan files in the /src folder.  The first is BERTA.stan and the second is BERTA_single_lake.stan.  These correspond to the original and single-lake versions of BERTA used in Cahill et al. 2021, respectively. The BERTA.stan file is not used, but is here in case folks want to compare with the single lake version. 

TODO for Chris: 
* work on markdown lessons 
* weep, drink more espresso 