# CHENOWETH PROJECT: MODEL COMPARISON IN CROSS-VALIDATION
# Jay Ulfelder
# 2015-01-05

# Housekeeping: clear workspace and set working directory
rm(list=ls(all=TRUE))
setwd("c:/users/jay/documents/nonviolent uprisings/replication/")

# Load required packages
library(caret)
library(reshape2)
library(plyr)
library(verification)
library(ROCR)
library(scoring)
library(DataCombine)

# Load full country-year data set with transformations but no listwise deletion
mash <- read.table("data/nvc.tranformed.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Load replication data set with folds for cross-validation and listwise deletion already done
valdat <- read.table("data/valdat.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

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
axis(1, pos = 0.5, at = c(1960,1970,1980,1990,2000,2010), las = 2, tick = FALSE)
dev.off()

# Plot for 1980-
png(file = "figs/nvc.onsets.by.year.1980.2013.png",
     width=10, height=8, units="cm", bg="white", res=150)
par(mar=c(3,2,1,1))
par(cex.axis = 0.8)
plot(seq(1980,2013,1), nvcbyyr[26:59], type="h", col="burlywood4", lwd=4, ylim=c(0, max(nvcbyyr) + 1),
  frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
axis(2, at = seq(0, max(nvcbyyr) + 1, by = 2), tick = FALSE, las = 2, line = -1)
axis(1, pos = 0.5, at = c(1980,1985,1990,1995,2000,2005,2010), las = 2, tick = FALSE)
dev.off()

################################################
# Cross-Validation
################################################

# Load model formulae
source("r/nvc.models.r")

# Create vector of variable names used in all models plus IDs
varlist <- c(c("country", "sftgcode", "year", "nvc.start.1"), unique(c(all.vars(gripes.f[[3]]), all.vars(modnzn.f[[3]]),
  all.vars(rscmob.f[[3]]), all.vars(polopp.f[[3]]))))

# Create function to generate predicted probabilities by iteration and fold of k-fold CV
predit <- function(i, x) {
  var <- paste0("ik", as.character(i))
  train <- valdat[ -which(valdat[var] == x),]
  test <- valdat[ which(valdat[var] == x),]
  test$base.p <- predict(glm(base.f, family = binomial, data = train, na.action = na.exclude), newdata = test, type = "response")
  test$gripes.p <- predict(glm(gripes.f, family = binomial, data = train, na.action = na.exclude), newdata = test, type = "response")
  test$modnzn.p <- predict(glm(modnzn.f, family = binomial, data = train, na.action = na.exclude), newdata = test, type = "response")
  test$rscmob.p <- predict(glm(rscmob.f, family = binomial, data = train, na.action = na.exclude), newdata = test, type = "response")
  test$polopp.p <- predict(glm(polopp.f, family = binomial, data = train, na.action = na.exclude), newdata = test, type = "response")
  test$iteration <- i
  test$k <- x
  out <- subset(test, select = c(country, year, nvc.start.1, base.p, gripes.p, modnzn.p, rscmob.p, polopp.p, iteration, k ))
  return(out)
}

# Apply that function iteratively to all 50 folds (10 iterations of 5-fold splitting)
out <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))

# Melt results to make data tidy
out.melt <- melt(out, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)

# Summarize accuracy by iteration
acc.by.iter <- ddply(out.melt, .(iteration, variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

# Summarize overall accuracy for each model
acc.overall2 <- ddply(out.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

# Summarize accuracy by model by year
acc.by.year <- ddply(out.melt, .(variable, year), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))))

# Plot accuracy by year
png(file = "figs/nvc.validation.auc.by.model.by.year.png",
    width=12, height=10, units='cm', bg='white', res=150)
par(mar=c(3,2,1,1))
plot(acc.by.year$log[acc.by.year$variable=="base.p"], type = "n", axes = FALSE, ylim = c(-0.3,0), ylab = "Brier score", xlab = "year")
points(acc.by.year$log[acc.by.year$variable=="base.p"], col = "black", pch = 19, cex = 0.5)
lines(predict(loess(acc.by.year$log[acc.by.year$variable=="base.p"] ~ acc.by.year$year[acc.by.year$variable=="base.p"])), col = "black", lwd = 2)
points(acc.by.year$log[acc.by.year$variable=="modnzn.p"], col = "dodgerblue", pch = 19, cex = 0.5)
lines(predict(loess(acc.by.year$log[acc.by.year$variable=="modnzn.p"] ~ acc.by.year$year[acc.by.year$variable=="modnzn.p"])), col = "dodgerblue", lwd = 2)
points(acc.by.year$log[acc.by.year$variable=="gripes.p"], col = "red", pch = 19, cex = 0.5)
lines(predict(loess(acc.by.year$log[acc.by.year$variable=="gripes.p"] ~ acc.by.year$year[acc.by.year$variable=="gripes.p"])), col = "red", lwd = 2)
points(acc.by.year$log[acc.by.year$variable=="rscmob.p"], col = "forestgreen", pch = 19, cex = 0.5)
lines(predict(loess(acc.by.year$log[acc.by.year$variable=="rscmob.p"] ~ acc.by.year$year[acc.by.year$variable=="rscmob.p"])), col = "forestgreen", lwd = 2)
points(acc.by.year$log[acc.by.year$variable=="polopp.p"], col = "goldenrod", pch = 19, cex = 0.5)
lines(predict(loess(acc.by.year$log[acc.by.year$variable=="polopp.p"] ~ acc.by.year$year[acc.by.year$variable=="polopp.p"])), col = "goldenrod", lwd = 2)
axis(1, at = c(4,9,14,19,24,29), labels = c("1985", "1990", "1995", "2000", "2005", "2010"), tick = FALSE)
axis(2, at = c(-0.3,-0.2,-0.1,0), labels = TRUE, tick = FALSE, las = 2, line = -1)
text(31, -0.24, "population only", col = "black", cex = 0.8, pos = 2)
text(31, -0.255, "grievances", col = "red", cex = 0.8, pos = 2)
text(31, -0.27, "modernization", col = "blue", cex = 0.8, pos = 2)
text(31, -0.285, "resource mobilization", col = "forestgreen", cex = 0.8, pos = 2)
text(31, -0.3, "political opportunity", col = "goldenrod", cex = 0.8, pos = 2)
dev.off()

#######################################
# Brute force version of LOO exercise
#######################################

# The code in this section iteratively reruns the accuracy-stats exercise, leaving out a different
# one of the model variables at each step, and then selects and combines the relevant stats from
# each iteration into a single table with labels. There has to be some way to do this more efficiently,
# but my R skills aren't up to the task.

#### GRIEVANCES ####

# gripes.f <- formula(nvc.start.1 ~ log(wdi.pop) +
#   log(xxxcimr) + wdi.gdpchg.s + sqrt(wdi.cpi) + log1p(bnn.yroff) + elceleth.c + dispota4.c + cir.physint + I(cir.physint^2) )

## im rate ##
gripes.f <- formula(nvc.start.1 ~ log(wdi.pop) + 
  wdi.gdpchg.s + sqrt(wdi.cpi) + log1p(bnn.yroff) + elceleth.c + dispota4.c + cir.physint + I(cir.physint^2) )
gripes.im <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
gripes.im.melt <- melt(gripes.im, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.gripes.im <- ddply(gripes.im.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

## gdp growth ##
gripes.f <- formula(nvc.start.1 ~ log(wdi.pop) + 
  log(xxxcimr) + sqrt(wdi.cpi) + log1p(bnn.yroff) + elceleth.c + dispota4.c + cir.physint + I(cir.physint^2) )
gripes.grow <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
gripes.grow.melt <- melt(gripes.grow, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.gripes.grow <- ddply(gripes.grow.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

## inflation ##
gripes.f <- formula(nvc.start.1 ~ log(wdi.pop) + 
  log(xxxcimr) + wdi.gdpchg.s + log1p(bnn.yroff) + elceleth.c + dispota4.c + cir.physint + I(cir.physint^2) )
gripes.cpi <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
gripes.cpi.melt <- melt(gripes.cpi, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.gripes.cpi <- ddply(gripes.cpi.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

## leader's tenure ##
gripes.f <- formula(nvc.start.1 ~ log(wdi.pop) + 
  log(xxxcimr) + wdi.gdpchg.s + sqrt(wdi.cpi) + elceleth.c + dispota4.c + cir.physint + I(cir.physint^2) )
gripes.yroff <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
gripes.yroff.melt <- melt(gripes.yroff, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.gripes.yroff <- ddply(gripes.yroff.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

## elite ethnicity ##
gripes.f <- formula(nvc.start.1 ~ log(wdi.pop) + 
  log(xxxcimr) + wdi.gdpchg.s + sqrt(wdi.cpi) + log1p(bnn.yroff) + dispota4.c + cir.physint + I(cir.physint^2) )
gripes.eleth <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
gripes.eleth.melt <- melt(gripes.eleth, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.gripes.eleth <- ddply(gripes.eleth.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

## discrimination ##
gripes.f <- formula(nvc.start.1 ~ log(wdi.pop) + 
  log(xxxcimr) + wdi.gdpchg.s + sqrt(wdi.cpi) + log1p(bnn.yroff) + elceleth.c + cir.physint + I(cir.physint^2) )
gripes.dis <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
gripes.dis.melt <- melt(gripes.dis, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.gripes.dis <- ddply(gripes.dis.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

## repression ##
gripes.f <- formula(nvc.start.1 ~ log(wdi.pop) + 
  log(xxxcimr) + wdi.gdpchg.s + sqrt(wdi.cpi) + log1p(bnn.yroff) + elceleth.c + dispota4.c )
gripes.physint <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
gripes.physint.melt <- melt(gripes.physint, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.gripes.physint <- ddply(gripes.physint.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

#### MODERNIZATION ####

# Reset model formulae
source("r/nvc.models.r")

# modnzn.f <- formula(nvc.start.1 ~ log(wdi.pop) +
#   wdi.popurb.mi + I(wdi.manuf.mi + wdi.servs.mi ) + wdi.sch2.mi + log1p(wdi.mobp100) + ios.gattwto ) 

## urbanization ##
modnzn.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  I(wdi.manuf.mi + wdi.servs.mi ) + wdi.sch2.mi + log1p(wdi.mobp100) + ios.gattwto )
modnzn.urb <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
modnzn.urb.melt <- melt(modnzn.urb, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.modnzn.urb <- ddply(modnzn.urb.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

## manufacturing and services ##
modnzn.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  wdi.popurb.mi + wdi.sch2.mi + log1p(wdi.mobp100) + ios.gattwto )
modnzn.econ <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
modnzn.econ.melt <- melt(modnzn.econ, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.modnzn.econ <- ddply(modnzn.econ.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

## education ##
modnzn.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  wdi.popurb.mi + I(wdi.manuf.mi + wdi.servs.mi ) + log1p(wdi.mobp100) + ios.gattwto )
modnzn.educ <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
modnzn.educ.melt <- melt(modnzn.educ, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.modnzn.educ <- ddply(modnzn.educ.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

## mobiles ##
modnzn.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  wdi.popurb.mi + I(wdi.manuf.mi + wdi.servs.mi ) + wdi.sch2.mi + ios.gattwto )
modnzn.mobl <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
modnzn.mobl.melt <- melt(modnzn.mobl, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.modnzn.mobl <- ddply(modnzn.mobl.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

## GATT/WTO ##
modnzn.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  wdi.popurb.mi + I(wdi.manuf.mi + wdi.servs.mi ) + wdi.sch2.mi + log1p(wdi.mobp100) )
modnzn.wto <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
modnzn.wto.melt <- melt(modnzn.wto, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.modnzn.wto <- ddply(modnzn.wto.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

#### RESOURCE MOBILIZATION ####

# Reset model formulae
source("r/nvc.models.r")

# rscmob.f <- formula(nvc.start.1 ~ log(wdi.pop) +
#   wdi.popurb.mi + ythbul4 + log1p(bnk.unrest) + log1p(bnk.strikes) + log1p(nvc.dosregt) + nvc.ongoing + civilwar )

# Urbanization #
rscmob.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  ythbul4 + log1p(bnk.unrest) + log1p(bnk.strikes) + log1p(nvc.dosregt) + nvc.ongoing + civilwar )
rscmob.urb <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
rscmob.urb.melt <- melt(rscmob.urb, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.rscmob.urb <- ddply(rscmob.urb.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

# Youth bulge #
rscmob.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  wdi.popurb.mi + log1p(bnk.unrest) + log1p(bnk.strikes) + log1p(nvc.dosregt) + nvc.ongoing + civilwar )
rscmob.ybul <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
rscmob.ybul.melt <- melt(rscmob.ybul, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.rscmob.ybul <- ddply(rscmob.ybul.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

# social unrest #
rscmob.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  wdi.popurb.mi + ythbul4 + log1p(bnk.strikes) + log1p(nvc.dosregt) + nvc.ongoing + civilwar )
rscmob.unrest <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
rscmob.unrest.melt <- melt(rscmob.unrest, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.rscmob.unrest <- ddply(rscmob.unrest.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

# general strikes #
rscmob.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  wdi.popurb.mi + ythbul4 + log1p(bnk.unrest) + log1p(nvc.dosregt) + nvc.ongoing + civilwar )
rscmob.strikes <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
rscmob.strikes.melt <- melt(rscmob.strikes, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.rscmob.strikes <- ddply(rscmob.strikes.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

# contagion #
rscmob.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  wdi.popurb.mi + ythbul4 + log1p(bnk.unrest) + log1p(bnk.strikes) + nvc.ongoing + civilwar )
rscmob.contagion <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
rscmob.contagion.melt <- melt(rscmob.contagion, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.rscmob.contagion <- ddply(rscmob.contagion.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

# ongoing #
rscmob.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  wdi.popurb.mi + ythbul4 + log1p(bnk.unrest) + log1p(bnk.strikes) + log1p(nvc.dosregt) + civilwar )
rscmob.ongoing <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
rscmob.ongoing.melt <- melt(rscmob.ongoing, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.rscmob.ongoing <- ddply(rscmob.ongoing.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

# civil war #
rscmob.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  wdi.popurb.mi + ythbul4 + log1p(bnk.unrest) + log1p(bnk.strikes) + log1p(nvc.dosregt) + nvc.ongoing )
rscmob.cwar <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
rscmob.cwar.melt <- melt(rscmob.cwar, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.rscmob.cwar <- ddply(rscmob.cwar.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

#### POLITICAL OPPORTUNITY ####

# Reset model formulae
source("r/nvc.models.r")

# polopp.f <- formula(nvc.start.1 ~ log(wdi.pop) +
#   log1p(age) + postcoldwar + ios.iccpr1 + nld.any.1 + pitfdem + I(pitfdem * nld.any.1) + I(postcoldwar * nld.any.1) +
#   as.factor(fiw.cl) + log1p(pol.durable) + log1p(cou.tries5) )

# age #
polopp.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  postcoldwar + ios.iccpr1 + nld.any.1 + pitfdem + I(pitfdem * nld.any.1) + I(postcoldwar * nld.any.1) +
  as.factor(fiw.cl) + log1p(pol.durable) + log1p(cou.tries5) )
polopp.age <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
polopp.age.melt <- melt(polopp.age, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.polopp.age <- ddply(polopp.age.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

# post-cold war #
polopp.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  log1p(age) + ios.iccpr1 + nld.any.1 + pitfdem + I(pitfdem * nld.any.1) +
  as.factor(fiw.cl) + log1p(pol.durable) + log1p(cou.tries5) )
polopp.pcw <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
polopp.pcw.melt <- melt(polopp.pcw, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.polopp.pcw <- ddply(polopp.pcw.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

# iccpr #
polopp.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  log1p(age) + postcoldwar + nld.any.1 + pitfdem + I(pitfdem * nld.any.1) + I(postcoldwar * nld.any.1) +
  as.factor(fiw.cl) + log1p(pol.durable) + log1p(cou.tries5) )
polopp.iccpr <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
polopp.iccpr.melt <- melt(polopp.iccpr, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.polopp.iccpr <- ddply(polopp.iccpr.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

# election year #
polopp.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  log1p(age) + postcoldwar + ios.iccpr1 + pitfdem +
  as.factor(fiw.cl) + log1p(pol.durable) + log1p(cou.tries5) )
polopp.elex <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
polopp.elex.melt <- melt(polopp.elex, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.polopp.elex <- ddply(polopp.elex.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

# democracy #
polopp.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  log1p(age) + postcoldwar + ios.iccpr1 + nld.any.1 + I(postcoldwar * nld.any.1) +
  as.factor(fiw.cl) + log1p(pol.durable) + log1p(cou.tries5) )
polopp.dem <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
polopp.dem.melt <- melt(polopp.dem, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.polopp.dem <- ddply(polopp.dem.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

# civil liberties #
polopp.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  log1p(age) + postcoldwar + ios.iccpr1 + nld.any.1 + pitfdem + I(pitfdem * nld.any.1) + I(postcoldwar * nld.any.1) +
  log1p(pol.durable) + log1p(cou.tries5) )
polopp.cl <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
polopp.cl.melt <- melt(polopp.cl, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.polopp.cl <- ddply(polopp.cl.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

# regime durability #
polopp.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  log1p(age) + postcoldwar + ios.iccpr1 + nld.any.1 + pitfdem + I(pitfdem * nld.any.1) + I(postcoldwar * nld.any.1) +
  as.factor(fiw.cl) + log1p(cou.tries5) )
polopp.regdur <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
polopp.regdur.melt <- melt(polopp.regdur, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.polopp.regdur <- ddply(polopp.regdur.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

# coup activity #
polopp.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  log1p(age) + postcoldwar + ios.iccpr1 + nld.any.1 + pitfdem + I(pitfdem * nld.any.1) + I(postcoldwar * nld.any.1) +
  as.factor(fiw.cl) + log1p(pol.durable) )
polopp.coup <- rbind(predit(1,1), predit(1,2), predit(1,3), predit(1,4), predit(1,5),
  predit(2,1), predit(2,2), predit(2,3), predit(2,4), predit(2,5),
  predit(3,1), predit(3,2), predit(3,3), predit(3,4), predit(3,5),
  predit(4,1), predit(4,2), predit(4,3), predit(4,4), predit(4,5),
  predit(5,1), predit(5,2), predit(5,3), predit(5,4), predit(5,5),
  predit(6,1), predit(6,2), predit(6,3), predit(6,4), predit(6,5),
  predit(7,1), predit(7,2), predit(7,3), predit(7,4), predit(7,5),
  predit(8,1), predit(8,2), predit(8,3), predit(8,4), predit(8,5),
  predit(9,1), predit(9,2), predit(9,3), predit(9,4), predit(9,5),
  predit(10,1), predit(10,2), predit(10,3), predit(10,4), predit(10,5))
polopp.coup.melt <- melt(polopp.coup, id = c("country", "year", "nvc.start.1", "k", "iteration"), na.rm = FALSE)
acc.polopp.coup <- ddply(polopp.coup.melt, .(variable), summarise,
  brier = mean((nvc.start.1 - value)^2),
  log = mean((nvc.start.1 * log(value)) + ((1 - nvc.start.1) * log(1 - value))),
  auc = roc.area(nvc.start.1, value)$A )

### TABLES ###

# Create table with accuracy stats for full models and models leaving one var out
loo <- rbind(
  acc.overall2[2,2:4],
  acc.gripes.im[2,2:4],
  acc.gripes.grow[2,2:4],
  acc.gripes.cpi[2,2:4],
  acc.gripes.yroff[2,2:4],
  acc.gripes.eleth[2,2:4],
  acc.gripes.dis[2,2:4],
  acc.gripes.physint[2,2:4],

  acc.overall2[3,2:4],
  acc.modnzn.urb[3,2:4],
  acc.modnzn.econ[3,2:4],
  acc.modnzn.educ[3,2:4],
  acc.modnzn.mobl[3,2:4],
  acc.modnzn.wto[3,2:4],

  acc.overall2[4,2:4],
  acc.rscmob.urb[4,2:4],
  acc.rscmob.ybul[4,2:4],
  acc.rscmob.unrest[4,2:4],
  acc.rscmob.strikes[4,2:4],
  acc.rscmob.contagion[4,2:4],
  acc.rscmob.ongoing[4,2:4],
  acc.rscmob.cwar[4,2:4],

  acc.overall2[5,2:4],
  acc.polopp.age[5,2:4],
  acc.polopp.pcw[5,2:4],
  acc.polopp.iccpr[5,2:4],
  acc.polopp.elex[5,2:4],
  acc.polopp.dem[5,2:4],
  acc.polopp.cl[5,2:4],
  acc.polopp.regdur[5,2:4],
  acc.polopp.coup[5,2:4] )

# Create vector of labels for those stats
loo.names <- c(
  "Grievances",
  "Less infant mortality",
  "Less GDP growth",
  "Less consumer price inflation",
  "Less leader's tenure",
  "Less elite ethnicity",
  "Less discrimination",
  "Less physical integrity",
  
  "Modernization",
  "Less urbanization",
  "Less manufacturing and services",
  "Less education",
  "Less mobile phone subscribers",
  "Less GATT/WTO membership",

  "Resource mobilization",
  "Less urbanization",
  "Less youth bulge",
  "Less social unrest",
  "Less general strikes",
  "Less contagion",
  "Less ongoing nonviolent campaign",
  "Less ongoing civil war",

  "Political opportunities",
  "Less country age",
  "Less post-cold war indicator",
  "Less ICCPR 1st Opt Protocol",
  "Less election year",
  "Less democracy",
  "Less civil liberties",
  "Less regime durability",
  "Less coup activity" )

# Create table with labels next to stats
loo <- cbind(loo.names, loo)

# Add categorical indicator for model to make next step easier
loo$model <- c(rep("Grievances", 8), rep("Modernization", 6), rep("Resource mobilization", 8), rep("Political opportunities", 9))

# Calculate percent change in AUC
loo$auc.pctchange <- NA
loo$auc.pctchange[loo$model == "Grievances"] <- 100 * (loo$auc[loo$model == "Grievances"] - acc.overall2$auc[2])
loo$auc.pctchange[loo$model == "Modernization"] <- 100 * (loo$auc[loo$model == "Modernization"] - acc.overall2$auc[3])
loo$auc.pctchange[loo$model == "Resource mobilization"] <- 100 * (loo$auc[loo$model == "Resource mobilization"] - acc.overall2$auc[4])
loo$auc.pctchange[loo$model == "Political opportunities"] <- 100 * (loo$auc[loo$model == "Political opportunities"] - acc.overall2$auc[5])

# Move string variables to the front to make next step easier
loo2 <- MoveFront(loo, c("loo.names", "model"))

# Round all numeric variables to nearest 1000th
loo2[,-1:-2] <- round(loo2[,-1:-2], 3)

# Inspect the results [Table 6 in the paper]
loo2

##########################################################
# Parameter estimates & descriptive stats for full panel
##########################################################

#### Models ####

# Reset model formulae
source("r/nvc.models.r")

# Re-estimate the models
gripes <- glm(gripes.f, family = binomial, data = valdat, na.action = na.exclude)
modnzn <- glm(modnzn.f, family = binomial, data = valdat, na.action = na.exclude)
rscmob <- glm(rscmob.f, family = binomial, data = valdat, na.action = na.exclude)
polopp <- glm(polopp.f, family = binomial, data = valdat, na.action = na.exclude)

# Inspect the results (Tables A2-A5 in the paper's appendix)
summary(gripes)
summary(modnzn)
summary(rscmob)
summary(polopp)

# Generate descriptive statistics for model variables from original data set (final table in appendix)
descstats <- summary(subset(mash, year >= 1981, select = varlist))

