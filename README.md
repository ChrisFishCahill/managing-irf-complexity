## Developing harvest control rules for Alberta Walleye fisheries

This directory contains Stan and R code to fit age structured population dynamics models to Alberta FWIN data as per Cahill et al. 2021.  It also contains .R simulation code that develops harvest control rules for Special Harvest License (SHL) Walleye fisheries given these model fits. The directory is straightforward, with a "data" folder, an "r-files" folder, a "plots" folder, and an "src" (i.e., source) folder for .stan files. 

The project plan is as follows:

* Create a single lake version of the BERTA models from Cahill et al. 2021 so that one can easily run models for a single lake 
* Develop the necessary programs to evaluate harvest control rules for a collection of lakes given model fits from part 1 
* Build capacity within AEP so that they can implement assessments and harvest control rule development on their own