# CHENOWETH PROJECT: MODEL COMPARISON IN CROSS-VALIDATION
# Jay Ulfelder
# 2015-01-05 --> 2020 refactor

options(stringsAsFactors = FALSE)

# Load required packages
library(tidyverse)
library(verification)
library(psych)

# Load full country-year data set with transformations but no listwise deletion
mash <- read.csv("data/nvc.transformed.csv")

# Load replication data set with folds for cross-validation and listwise deletion already done
valdat <- read.csv("data/valdat.csv")

####################################
# Plots of DV
####################################

# Generate annual counts
nvcbyyr <- tapply(mash$nvc.start, mash$year, sum, na.rm = TRUE)

# Plot for full data set
png(file = "figs/nvc.onsets.by.year.1955.2013.png",
   width=14, height=8, units="cm", bg="white", res=150)
par(mar=c(3,2,1,1))
par(cex.axis = 0.8)
plot(seq(1955,2013,1), nvcbyyr, type="h", col="burlywood4", lwd=3.5, ylim=c(0,max(nvcbyyr) + 1),
  frame.plot=FALSE, axes = FALSE, xlab="", ylab="")
axis(2, at = seq(0, max(nvcbyyr) + 1, by = 2), tick = FALSE, las = 2, line = -1.5)
axis(1, pos = 0.5, at = seq(1960, 2010, by = 10), las = 2, tick = FALSE)
dev.off()

# Plot for 1980-
png(file = "figs/nvc.onsets.by.year.1980.2013.png",
     width=10, height=8, units="cm", bg="white", res=150)
par(mar=c(3,2,1,1))
par(cex.axis = 0.8)
plot(seq(1980,2013,1), nvcbyyr[which(as.numeric(names(nvcbyyr)) >= 1980)],
  type="h", col="burlywood4", lwd=4, ylim=c(0, max(nvcbyyr) + 1),
  frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
axis(2, at = seq(0, max(nvcbyyr) + 1, by = 2), tick = FALSE, las = 2, line = -1)
axis(1, pos = 0.5, at = c(1980,1985,1990,1995,2000,2005,2010), las = 2, tick = FALSE)
dev.off()

################################################
# Cross-Validation
################################################

# This block of code implements 5 iterations of 10-fold cross-validation for each of the specified models, then
# summarizes the results overall and for some subsets of interest. It depends on a sampling process that was
# previously implemented and stored in the data set, to ensure replicability across different software tools.
# Each row's membership in a particular fold (k) for a particular iteration (i) of k-fold CV is indicated by 
# the values of the columns 'ik1' through 'ik10'. The suffixes in those column names identify the folds, while
# the values in those columns identify the iteration. So, for example ik2 == 3 indicates that that row belongs
# to the 2nd fold in the 3rd iteration.

# Use cross and pmap from purrr to create a table with all the desired combinations of i, k,
# and model and then iterate over them a function that generates oos predictions for each of those combos.
# The result is tidy long, which is handy.
cv_predit <- function(my_i = 5, my_k = 10, model_list) {

    require(dplyr)
    require(purrr)

    predit <- function(i, k, model_name) {

        var <- sprintf("ik%s", k)

        train <- valdat[ -which(valdat[var] == i),]
        test <- valdat[ which(valdat[var] == i),]

        mod <- glm(model_list[[model_name]], family = "binomial", data = train, na.action = na.exclude)

        test$pred <- predict(mod, newdata = test, type = "response")

        test$iteration <- i
        test$k <- k
        test$model <- model_name

        output <- subset(test, select = c(country, year, nvc.start.1, iteration, k, model, pred))

        return(output)

    }

    crossArg <- cross_df(list(i = seq(my_i), k = seq(my_k), model = names(model_list)))
    results <- pmap_dfr(with(crossArg, list(i, k, model)), predit)

    return(results)

}

# helper functions for measuring accuracy
f_brier <- function(y, y_hat) { mean((y - y_hat)^2) }
f_log <- function(y, y_hat) { mean((y * log(y_hat)) + ((1 - y) * log(1 - y_hat))) }
f_auc <- function(y, y_hat) { verification::roc.area(y, y_hat)$A }

# Load model formulae
source("r/nvc.models.r")

# make a named list with model formulae over which we can iterate
list_formulae <- list(base = base.f, gripes = gripes.f, modnzn = modnzn.f, rscmob = rscmob.f, polopp = polopp.f)

# apply the cv functions to it
cv_results <- cv_predit(model_list = list_formulae)

# Summarize accuracy of each model by iteration
acc.by.iter <- cv_results %>%
    group_by(iteration, model) %>%
    summarise(brier = f_brier(nvc.start.1, pred), log = f_log(nvc.start.1, pred), auc = f_auc(nvc.start.1, pred))

# Summarize overall accuracy of each model
# NOTE: Original version did this with results pooled by model. The correct approach is to average the stats
# for each model across all iterations of k-fold cv, which is what's done here. There results are virtually
# identical, but I'm a stickler, so...
acc.overall <- acc.by.iter %>% 
    ungroup() %>%
    group_by(model) %>%
    summarise_at(vars(brier, log, auc), mean) 

# Summarize accuracy by model by year
acc.by.year <- cv_results %>%
    group_by(model, year) %>%
    summarise(brier = f_brier(nvc.start.1, pred), log = f_log(nvc.start.1, pred))

# Plot accuracy by year
png(file = "figs/nvc.validation.auc.by.model.by.year.png",
    width=12, height=10, units='cm', bg='white', res=150)
par(mar=c(3,2,1,1))
plot(acc.by.year$log[acc.by.year$model=="base"], type = "n", axes = FALSE, ylim = c(-0.3,0), ylab = "Brier score", xlab = "year")
points(acc.by.year$log[acc.by.year$model=="base"], col = "black", pch = 19, cex = 0.5)
lines(predict(loess(acc.by.year$log[acc.by.year$model=="base"] ~ acc.by.year$year[acc.by.year$model=="base"])), col = "black", lwd = 2)
points(acc.by.year$log[acc.by.year$model=="modnzn"], col = "dodgerblue", pch = 19, cex = 0.5)
lines(predict(loess(acc.by.year$log[acc.by.year$model=="modnzn"] ~ acc.by.year$year[acc.by.year$model=="modnzn"])), col = "dodgerblue", lwd = 2)
points(acc.by.year$log[acc.by.year$model=="gripes"], col = "red", pch = 19, cex = 0.5)
lines(predict(loess(acc.by.year$log[acc.by.year$model=="gripes"] ~ acc.by.year$year[acc.by.year$model=="gripes"])), col = "red", lwd = 2)
points(acc.by.year$log[acc.by.year$model=="rscmob"], col = "forestgreen", pch = 19, cex = 0.5)
lines(predict(loess(acc.by.year$log[acc.by.year$model=="rscmob"] ~ acc.by.year$year[acc.by.year$model=="rscmob"])), col = "forestgreen", lwd = 2)
points(acc.by.year$log[acc.by.year$model=="polopp"], col = "goldenrod", pch = 19, cex = 0.5)
lines(predict(loess(acc.by.year$log[acc.by.year$model=="polopp"] ~ acc.by.year$year[acc.by.year$model=="polopp"])), col = "goldenrod", lwd = 2)
axis(1, at = c(4,9,14,19,24,29), labels = c("1985", "1990", "1995", "2000", "2005", "2010"), tick = FALSE)
axis(2, at = c(-0.3,-0.2,-0.1,0), labels = TRUE, tick = FALSE, las = 2, line = -1)
text(31, -0.24, "population only", col = "black", cex = 0.8, pos = 2)
text(31, -0.255, "grievances", col = "red", cex = 0.8, pos = 2)
text(31, -0.27, "modernization", col = "blue", cex = 0.8, pos = 2)
text(31, -0.285, "resource mobilization", col = "forestgreen", cex = 0.8, pos = 2)
text(31, -0.3, "political opportunity", col = "goldenrod", cex = 0.8, pos = 2)
dev.off()

#######################################
# LOO exercise
#######################################

# The code in this section iteratively reruns the accuracy-stats exercise, leaving out a different
# one of the model variables at each step, and then selects and combines the relevant stats from
# each iteration into a single table with labels. It produces Table 6 in the paper.
#
# 2020 refactor: for each model, 1) make a named list of formulae with all the permutations (and names that
# are informative about them, e.g., 'gripes - im', to use as labels in the results); 2) use the
# same cross/pmap/predit combo used in the previous section to get lists or dfs with the results; then
# 3) create and combine tables summarizing accuracy.
#
# NOTE: In the refactoring, I changed the summarization of the cross-validation from stats for results
# pooled across all iterations to averages of stats for each iteration of 10-fold CV. The results are
# virtually indistinguishable, but the two-step approach is the correct one.

#### GRIEVANCES ####

loo_mods_gripes <- list(

    `gripes - im` = formula(nvc.start.1 ~ log(wdi.pop) + wdi.gdpchg.s + sqrt(wdi.cpi) + log1p(bnn.yroff) + elceleth.c + dispota4.c + cir.physint + I(cir.physint^2) ),

    `gripes - gdpchg` = formula(nvc.start.1 ~ log(wdi.pop) + log(xxxcimr) + sqrt(wdi.cpi) + log1p(bnn.yroff) + elceleth.c + dispota4.c + cir.physint + I(cir.physint^2) ),

    `gripes - inflation` = formula(nvc.start.1 ~ log(wdi.pop) + log(xxxcimr) + wdi.gdpchg.s + log1p(bnn.yroff) + elceleth.c + dispota4.c + cir.physint + I(cir.physint^2) ),

    `gripes - tenure` =  formula(nvc.start.1 ~ log(wdi.pop) + log(xxxcimr) + wdi.gdpchg.s + sqrt(wdi.cpi) + elceleth.c + dispota4.c + cir.physint + I(cir.physint^2) ),

    `gripes - ethinicity` = formula(nvc.start.1 ~ log(wdi.pop) + log(xxxcimr) + wdi.gdpchg.s + sqrt(wdi.cpi) + log1p(bnn.yroff) + dispota4.c + cir.physint + I(cir.physint^2) ),

    `gripes - discrimination` = formula(nvc.start.1 ~ log(wdi.pop) + log(xxxcimr) + wdi.gdpchg.s + sqrt(wdi.cpi) + log1p(bnn.yroff) + elceleth.c + cir.physint + I(cir.physint^2) ),

    `gripes - repression` = formula(nvc.start.1 ~ log(wdi.pop) + log(xxxcimr) + wdi.gdpchg.s + sqrt(wdi.cpi) + log1p(bnn.yroff) + elceleth.c + dispota4.c )

)

loo_acc_gripes <- cv_predit(model_list = loo_mods_gripes) %>%
    group_by(model, iteration) %>%
    summarise(brier = f_brier(nvc.start.1, pred), log = f_log(nvc.start.1, pred), auc = f_auc(nvc.start.1, pred)) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise_at(vars(brier, log, auc), mean)

#### MODERNIZATION ####

loo_mods_modnzn <- list(

    `modnzn - urbanization` = formula(nvc.start.1 ~ log(wdi.pop) + I(wdi.manuf.mi + wdi.servs.mi ) + wdi.sch2.mi + log1p(wdi.mobp100) + ios.gattwto ),

    `modnzn - manufacturing` = formula(nvc.start.1 ~ log(wdi.pop) + wdi.popurb.mi + wdi.sch2.mi + log1p(wdi.mobp100) + ios.gattwto ),

    `modnzn - education` = formula(nvc.start.1 ~ log(wdi.pop) + wdi.popurb.mi + I(wdi.manuf.mi + wdi.servs.mi ) + log1p(wdi.mobp100) + ios.gattwto ),

    `modnzn - cellphones` = formula(nvc.start.1 ~ log(wdi.pop) + wdi.popurb.mi + I(wdi.manuf.mi + wdi.servs.mi ) + wdi.sch2.mi + ios.gattwto ),

    `modnzn - gattwto` = formula(nvc.start.1 ~ log(wdi.pop) + wdi.popurb.mi + I(wdi.manuf.mi + wdi.servs.mi ) + wdi.sch2.mi + log1p(wdi.mobp100) )

)

loo_acc_modnzn <- cv_predit(model_list = loo_mods_modnzn) %>%
    group_by(model, iteration) %>%
    summarise(brier = f_brier(nvc.start.1, pred), log = f_log(nvc.start.1, pred), auc = f_auc(nvc.start.1, pred)) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise_at(vars(brier, log, auc), mean)

#### RESOURCE MOBILIZATION ####

loo_mods_rscmob <- list(

    `rscmob - urbanization` = formula(nvc.start.1 ~ log(wdi.pop) + ythbul4 + log1p(bnk.unrest) + log1p(bnk.strikes) + log1p(nvc.dosregt) + nvc.ongoing + civilwar ),

    `rscmob - youthbulge` = formula(nvc.start.1 ~ log(wdi.pop) + wdi.popurb.mi + log1p(bnk.unrest) + log1p(bnk.strikes) + log1p(nvc.dosregt) + nvc.ongoing + civilwar ),

    `rscmob - socialunrest` = formula(nvc.start.1 ~ log(wdi.pop) + wdi.popurb.mi + ythbul4 + log1p(bnk.strikes) + log1p(nvc.dosregt) + nvc.ongoing + civilwar ),

    `rscmob - strikes` = formula(nvc.start.1 ~ log(wdi.pop) + wdi.popurb.mi + ythbul4 + log1p(bnk.unrest) + log1p(nvc.dosregt) + nvc.ongoing + civilwar ),

    `rscmob - contagion` = formula(nvc.start.1 ~ log(wdi.pop) + wdi.popurb.mi + ythbul4 + log1p(bnk.unrest) + log1p(bnk.strikes) + nvc.ongoing + civilwar ),

    `rscmob - ongoingnvc` = formula(nvc.start.1 ~ log(wdi.pop) + wdi.popurb.mi + ythbul4 + log1p(bnk.unrest) + log1p(bnk.strikes) + log1p(nvc.dosregt) + civilwar ),

    `rscmob - civilwar` = formula(nvc.start.1 ~ log(wdi.pop) + wdi.popurb.mi + ythbul4 + log1p(bnk.unrest) + log1p(bnk.strikes) + log1p(nvc.dosregt) + nvc.ongoing )

)

loo_acc_rscmob <- cv_predit(model_list = loo_mods_rscmob) %>%
    group_by(model, iteration) %>%
    summarise(brier = f_brier(nvc.start.1, pred), log = f_log(nvc.start.1, pred), auc = f_auc(nvc.start.1, pred)) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise_at(vars(brier, log, auc), mean)

#### POLITICAL OPPORTUNITY ####

loo_mods_polopp <- list(

    `polopp - age` = formula(nvc.start.1 ~ log(wdi.pop) + postcoldwar + ios.iccpr1 + nld.any.1 + pitfdem + I(pitfdem * nld.any.1) + I(postcoldwar * nld.any.1) + as.factor(fiw.cl) + log1p(pol.durable) + log1p(cou.tries5) ),

    `polopp - postcoldwar` = formula(nvc.start.1 ~ log(wdi.pop) + log1p(age) + ios.iccpr1 + nld.any.1 + pitfdem + I(pitfdem * nld.any.1) + as.factor(fiw.cl) + log1p(pol.durable) + log1p(cou.tries5) ),

    `polopp - iccpr` = formula(nvc.start.1 ~ log(wdi.pop) + log1p(age) + postcoldwar + nld.any.1 + pitfdem + I(pitfdem * nld.any.1) + I(postcoldwar * nld.any.1) + as.factor(fiw.cl) + log1p(pol.durable) + log1p(cou.tries5) ),

    `polopp - election` = formula(nvc.start.1 ~ log(wdi.pop) + log1p(age) + postcoldwar + ios.iccpr1 + pitfdem + as.factor(fiw.cl) + log1p(pol.durable) + log1p(cou.tries5) ),

    `polopp - democracy` = formula(nvc.start.1 ~ log(wdi.pop) + log1p(age) + postcoldwar + ios.iccpr1 + nld.any.1 + I(postcoldwar * nld.any.1) + as.factor(fiw.cl) + log1p(pol.durable) + log1p(cou.tries5) ),

    `polopp - civlibs` = formula(nvc.start.1 ~ log(wdi.pop) + log1p(age) + postcoldwar + ios.iccpr1 + nld.any.1 + pitfdem + I(pitfdem * nld.any.1) + I(postcoldwar * nld.any.1) + log1p(pol.durable) + log1p(cou.tries5) ),

    `polopp - regimedur` = formula(nvc.start.1 ~ log(wdi.pop) + log1p(age) + postcoldwar + ios.iccpr1 + nld.any.1 + pitfdem + I(pitfdem * nld.any.1) + I(postcoldwar * nld.any.1) + as.factor(fiw.cl) + log1p(cou.tries5) ),

    `polopp - coups` = formula(nvc.start.1 ~ log(wdi.pop) + log1p(age) + postcoldwar + ios.iccpr1 + nld.any.1 + pitfdem + I(pitfdem * nld.any.1) + I(postcoldwar * nld.any.1) + as.factor(fiw.cl) + log1p(pol.durable) )

)

loo_acc_polopp <- cv_predit(model_list = loo_mods_polopp) %>%
    group_by(model, iteration) %>%
    summarise(brier = f_brier(nvc.start.1, pred), log = f_log(nvc.start.1, pred), auc = f_auc(nvc.start.1, pred)) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise_at(vars(brier, log, auc), mean)

### TABLES ###

# Create table with accuracy stats for full models and models leaving one var out
loo <- bind_rows(list(

    rbind(acc.overall[acc.overall$model == "gripes",], loo_acc_gripes),

    rbind(acc.overall[acc.overall$model == "modnzn",], loo_acc_modnzn),
    
    rbind(acc.overall[acc.overall$model == "rscmob",], loo_acc_rscmob),

    rbind(acc.overall[acc.overall$model == "polopp",], loo_acc_polopp)

))

# Create vector of labels for those stats; order matters!
loo.names <- c(

  # grievances
  "Full model",
  "Less discrimination",
  "Less elite ethnicity",
  "Less GDP growth",
  "Less infant mortality",
  "Less consumer price inflation",
  "Less physical integrity",
  "Less leader's tenure",
  
  # modernization
  "Full model",
  "Less mobile phone subscribers",
  "Less education",
  "Less GATT/WTO membership",
  "Less manufacturing and services",
  "Less urbanization",

  # resource mobilization
  "Full model",
  "Less ongoing civil war",
  "Less contagion",
  "Less ongoing nonviolent campaign",
  "Less social unrest",
  "Less general strikes",
  "Less urbanization",
  "Less youth bulge",

  # political opportunities
  "Full model",
  "Less country age",
  "Less civil liberties",
  "Less coup activity",
  "Less democracy",
  "Less election year",
  "Less ICCPR 1st Opt Protocol",
  "Less post-cold war indicator",
  "Less regime durability"

)

# Create table with labels next to stats
loo <- cbind(variation = loo.names, loo)

# Add cleaner indicator for model to make next step easier
loo$model[grepl("gripes", loo$model)] <- "Grievances"
loo$model[grepl("modnzn", loo$model)] <- "Modernization"
loo$model[grepl("rscmob", loo$model)] <- "Resource mobilization"
loo$model[grepl("polopp", loo$model)] <- "Political opportunities"

# Calculate percent change in AUC
loo$auc.pctchange <- NA
loo$auc.pctchange[loo$model == "Grievances"] <- 100 * (loo$auc[loo$model == "Grievances"] - loo$auc[loo$model == "Grievances" & loo$variation == "Full model"])
loo$auc.pctchange[loo$model == "Modernization"] <- 100 * (loo$auc[loo$model == "Modernization"] - loo$auc[loo$model == "Modernization" & loo$variation == "Full model"])
loo$auc.pctchange[loo$model == "Resource mobilization"] <- 100 * (loo$auc[loo$model == "Resource mobilization"] - loo$auc[loo$model == "Resource mobilization" & loo$variation == "Full model"])
loo$auc.pctchange[loo$model == "Political opportunities"] <- 100 * (loo$auc[loo$model == "Political opportunities"] - loo$auc[loo$model == "Political opportunities" & loo$variation == "Full model"])

# tidy up the results for presentation
loo <- loo %>%
    dplyr::select(model, variation, everything()) %>%
    mutate(auc.pctchange = sprintf("%s %%", round(auc.pctchange, 1))) %>%
    mutate(auc.pctchange = ifelse(variation == "Full model", "", auc.pctchange)) %>%
    mutate_if(is.numeric, ~ round(., 3))

# Inspect the results [Table 6 in the paper]
print(loo)

##########################################################
# Parameter estimates & descriptive stats for full panel
##########################################################

# Re-estimate the models
full_results <- purrr::map(list_formulae, ~ glm(., family = binomial, data = valdat, na.action = na.exclude))

# Inspect the results (Tables A2-A5 in the paper's appendix)
purrr::map(full_results, summary)

# Generate descriptive statistics for model variables from original data set (final table in appendix)
descstats <- mash %>%
    filter(year >= 1981) %>%
    mutate(age = year - yrborn) %>%  # this country age measure used in gripes is not precomputed in 'mash' for some reason
    dplyr::select(one_of(c("nvc.start.1", unique(c(all.vars(gripes.f[[3]]), all.vars(modnzn.f[[3]]), all.vars(rscmob.f[[3]]), all.vars(polopp.f[[3]])))))) %>%
    psych::describe(.)
