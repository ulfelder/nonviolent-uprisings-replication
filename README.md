This repository contains replication materials for Erica Chenoweth and Jay Ulfelder, "Can Structural Conditions Explain the Onset of Nonviolent Uprisings?" *Journal of Conflict Resolution* Vol. 61, No 2 (February 2017), pp. 298-324 ([here](https://journals.sagepub.com/doi/full/10.1177/0022002715576574)). Please direct questions about the data or code to ulfelder at gmail dot com.

The replication materials comprise three R scripts and two data sets.

Two of the R scripts ("nvc.crossvalidation.R" and "nvc.prediction.R") execute different parts of the analysis; a third ("nvc.models.R") creates formula objects for the various models and is called by the other two. The order in which the two executable scripts are run does not matter.

The executing scripts call two versions of the project's data set. For cross-validation, we use "valdat.csv", a subset of our original analysis file that a) includes variables marking samples for 10 iterations of five-fold cross-validation and b) has already gone through listwise deletion of incomplete cases. For quasi-prediction, we use "nvc.transformed.csv", the original country-year data set from which "valdat.csv" is derived.

This work was done in R Version 3.1.1 for 64-bit Windows OS. Package dependencies are indicated at the top of the executable scripts.

**NOTE:** In early 2020, I (Jay) refactored the "nvc.crossvalidation.R" script, mostly as a practice project to address some inefficiencies in the original approach. The refactored version is stored here as "R/nvc.crossvalidation-2020-refactor.R". Its output is virtually identical to the original's (see comments in the script for why some results differ at the 0.001s), but it gets there with one-third as many lines of code. If your objective is to precisely reproduce the estimates reported in the paper, please use the original version. If you're looking to replicate or extend the process instead, please use the refactored version.
