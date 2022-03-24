## Developing harvest control rules for Alberta Walleye fisheries with highly variable recruitment dynamics
Chris Cahill and Carl Walters, March 2022

### Directory organization
This directory contains data, Stan, and R code to fit age structured population dynamics models to Alberta FWIN data as per Cahill et al. 2021.  It also contains .R simulation code that develops harvest control rules for Special Harvest License (SHL) Walleye fisheries given these model fits. 

The directory is straightforward, with a "data" folder, an "r-code" folder, a "plots" folder, and an "src" (i.e., source) folder for .stan files. Tutorials are in the "Rmd" folder, and harvest control rule simulations are stored in the "sims" folder. The plots are mostly a work in progress, so be careful trying to interpret them (i.e., don't do this).  It was simply a location for me to share plots with Carl and others.  The scripts in the "r-code" folder that are potentially useful to folks are the `run.R` (calls and runs the stan file, saves fits), `hcr.R` (runs the harvest control rule simulator, saves simulations), and plotting files (calls either the Bayesian model fits or the simulation results). 

Note this directory has stan fit objects on it, which make it pretty big and thus cumbersome to download.  It is probably best to download the repostiory as a zip and play with it that way. Also, it is best to work from the .pdf tutorials in the `Rmd/` folder rather than the `.rmd` files directly. FYI.

### Background
In Alberta there are Walleye fisheries managed using a "Special Harvest Licence" or SHL harvest tag program.  This program has regional managers allocate some number of harvest tags, which anglers can then purchase through a lottery system identical to the system currently used to manage hunting tags for non-finned critters like Whitetail Deer, Moose, or tasty Bighorn Sheep.  In these systems, if one wants to retain a Walleye one must first apply for a tag, be selected in the harvest lottery, and then purchase that  harvest tag.  Otherwise, any fish an angler catches must be released (note this means total effort is NOT limited in these fisheries).

This management system arose from the career-long efforts of Michael Sullivan (2003), Stephen Spencer, and many other Alberta biologists.  These stringent restrictions were implemented in response to recognition that Alberta Walleye apparently underwent an "Invisible Collapse" (Post et al. 2001; Sullivan 2003). The SHL program has been in place since the early 2000s, and started in a few (< 5) high effort lakes.  The program has since expanded to ~ 20 systems in the province. 

A recent paper used age-structured population dynamics models to assess Walleye populations throughout Alberta (see Cahill et al. 2021 CJFAS).  A major finding of our work was the documentation of foregone harvest opportunities in many lakes during 2000-2018, with the caveat that we estimated low Fmsy values for all fisheries (i.e., the instantaneous fishing mortality that would achieve maximum sustainable yield).  These low Fmsy values mean that while there were foregone harvest opportunities during 2000-2018, Alberta Walleye still remain vulnerable to overharvest.  Additionally, we showed that Alberta Walleye recruitment was highly variable and strongly pulsed in some lakes.  Given these findings, we argued that management strategy evaluation (MSE) would be useful for informing managers in Alberta. 

MSE is a powerful technique with which to design effective fisheries management systems (see Punt et al. 2016 for details of the approach).  This is particularly true for Walleye fisheries managed via a total allowable catch in terms of the number of tags managers allocate in a given year. The primary goal of this project, then, was to develop defensible harvest control rules for Alberta managers that were robust to this (as of yet) unpredictable recruitment variability.  To do this, we used the assessment results from Cahill et al. 2021 (i.e., the Bayesian posteriors from their fitted models) to parameterize a harvest control rule simulation. 

The simulation tests simple linear management policies or harvest control rules to find policies that maximize simple objectives like maximum average yield or HARA utility (defined explicitly in tutorial on harvest control rules).  We chose conflicting but simple objectives to demonstrate how these fisheries could be managed rather than show folks how they should be managed.

### Changes 

There have been three changes to the .stan model from what was published in Cahill et al. 2021:

* We now assume logistic vulnerability to fishing, which was not used in the paper.
* We have assumed a higher Fearly mortality rate prior of 0.5 (vs. 0.3 used in the paper). 
* We have corrected the estimates of MSY and Fmsy so that they now account for the lognormal bias. 

We made these changes to ensure the development of our harvest control rules was robust and more conservative from a fishing mortality perspective.  None of the major findings from Cahill et al. 2021 changed when we made these additional changes. 

The project plan is as follows:

* Create a single lake version of the BERTA models from Cahill et al. 2021 so that one can easily run models for a single lake.
* Develop the necessary programs to evaluate harvest control rules for a collection of lakes given model fits from part 1. 
* Build capacity within AEP so that they can implement assessments and harvest control rule development on their own (via tutorials).
* Teach Andrew Paul and Laura MacPherson how to do this via Zoom (taught during Jan - March 2022).

Additional notes:

* There are two .stan files in the /src folder.  The first is BERTA.stan and the second is BERTA_single_lake.stan.  These correspond to the original and single-lake versions of BERTA used in Cahill et al. 2021, respectively. The BERTA.stan file is not used, but is here in case folks want to compare with the single lake version. 

TODO for Chris: 
* drink more espresso 