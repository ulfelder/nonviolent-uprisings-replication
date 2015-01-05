This repository contains replication materials for Erica Chenoweth and Jay Ulfelder, "Can Structural Conditions Explain the Onset of Nonviolent Uprisings?" Please direct questions or suggestions to Jay at ulfelder at gmail dot com.

The replication materials comprise three R scripts and two data sets.

Two of the R scripts ("nvc.crossvalidation.R" and "nvc.prediction.R") execute different parts of the analysis; a third ("nvc.models.R") creates formula objects for the various models and is called by the other two. The order in which the two executable scripts are run does not matter.

The executing scripts call two versions of the project's data set. For cross-validation, we use "valdat.csv", a subset of our original analysis file that a) includes variables marking samples for 10 iterations of five-fold cross-validation and b) has already gone through listwise deletion of incomplete cases. For quasi-prediction, we use "nvc.transformed.csv", the original country-year data set from which "valdat.csv" is derived.

This work was done in R Version 3.1.1 for 64-bit Windows OS. Package dependencies are indicated at the top of the two executable scripts.
