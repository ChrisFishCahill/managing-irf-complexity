## Developing harvest control rules for Alberta Walleye fisheries
Chris Cahill & Carl Walters, March 2022

### Background
In Alberta there are Walleye fisheries managed using a so-called "Special Harvest Licence" or SHL program.  This program has regional managers allocate some number of harvest tags, which anglers can then purchase through a lottery system identical to the system currently used to manage hunting tags for critters like Whitetail Deer, Moose, or tasty Bighorn Sheep.  In these systems, if one wants to retain a Walleye one must first apply for a tag, be selected in the harvest lottery, and then purchase that  harvest tag.  Otherwise, any fish an angler catches must be released (note this means total effort is NOT limited in these fisheries).

This management system arose from the career-long work efforts of Michael Sullivan (2003), Stephen Spencer, and many other Alberta biologists.  In our opinion, this program represents an important step forward for managing open access inland fisheries in North America. The SHL program has been running since the early 2000s, and started in a few (< 5) high effort lakes.  The program has since expanded to more ~ 20 systems throughout the province. 

A recent paper used age-structured population dynamics to assess Walleye populations (not just SHL managed fisheries) throughout the province (see Cahill et al. 2021 CJFAS).  A major finding of this work was the documentation of foregone harvest opportunities in many lakes during 2000-2018, with the caveat that assessments estimated low Fmsy values for all fisheries (i.e., the instantaneous fishing mortality that would achieve maximum sustainable yield).  These low Fmsy values mean that while there were foregone harvest opportunities during 2000-2018, Alberta Walleye still remain vulnerable to overharvest.  Additionally, their paper showed that Alberta Walleye recruitment was highly variable and strongly pulsed in some lakes.  Given these findings, we argued that management strategy evaluation would be a powerful way forward for informing defensible management strategies (MSE) in Alberta fisheries. 

MSE is a powerful technique with which to design effective fisheries management system (see Punt et al. 2016).  This is particularly true for Walleye fisheries managed via a total allowable catch in terms of the number of tags managers allocate in a given year. The primary goal of this project, then, was to develop defensible harvest control rules for Alberta managers that were robust to this (as of yet) unpredictable recruitment variability.  To do this, we use the assessment results from Cahill et al. 2021 (i.e., the Bayesian posteriors from their fitted models) to parameterize a harvest control rule simulation.  The simulation tests a ton of different simple management policies to find policies that maximize simple objectives like maximum average yield or HARA utility (defined in tutorial on harvest control rules).

### Directory organization
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