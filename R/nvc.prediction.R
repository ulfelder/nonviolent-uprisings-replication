# CHENOWETH PROJECT: QUASI-PREDICTION
# Jay Ulfelder
# 2015-01-05

# Clear workspace.
rm(list=ls(all=TRUE))

# Load required packages
library(verification)
library(scoring)
library(lattice)
library(Hmisc)

# Load data
mash <- read.csv("data/nvc.transformed.csv")

# Load model formulae
source("r/nvc.models.r")

# Add columns with country age and prediction year (observation year + 1)
mash$age <- mash$year - mash$yrborn
mash$predyr <- mash$year + 1

# Create formula for "culled" model using vars selected on basis of cross-validation
culled.f <- formula(nvc.start.1 ~ log(wdi.pop) +
  log(xxxcimr) + log1p(bnn.yroff) + elceleth.c + 
  log1p(wdi.mobp100) +  
  wdi.popurb.mi + ythbul4 + log1p(bnk.unrest) + log1p(nvc.dosregt) +
  postcoldwar + nld.any.1 + I(postcoldwar * nld.any.1) + as.factor(fiw.cl) + ios.iccpr1)

#############################################
# Some easy updates
#############################################

# International organizations
# Updates source: http://treaties.un.org/Pages/ViewDetails.aspx?src=TREATY&mtdsg_no=IV-5&chapter=4&lang=en
for (i in 1:dim(mash)[1]) mash$ios.iccpr1[i] <- ifelse(mash$year[i] > 2010, mash$ios.iccpr1[i-1], mash$ios.iccpr1[i])
mash$ios.iccpr1[mash$sftgcode=="TUN" & mash$year>=2011] <- 1 # acceded in 2011
mash$ios.iccpr1[mash$sftgcode=="USA"] <- 0 
mash$ios.iccpr1[mash$sftgcode=="SSD"] <- 0 
mash$ios.iccpr1[mash$sftgcode=="MNG"] <- 0 
mash$ios.iccpr1[mash$sftgcode=="SWA"] <- 0 
mash$ios.iccpr1[mash$sftgcode=="GNB" & mash$year>=2013] <- 1  # acceded in 2013
# Source: http://www.wto.org/english/thewto_e/whatis_e/tif_e/org6_e.htm
for (i in 1:dim(mash)[1]) mash$ios.gattwto[i] <- ifelse(mash$year[i] > 2010, mash$ios.gattwto[i-1], mash$ios.gattwto[i])
mash$ios.gattwto[mash$sftgcode=="MNG" & mash$year>=2012] <- 1  # acceded in 2012
mash$ios.gattwto[mash$sftgcode=="RUS" & mash$year>=2012] <- 1  # acceded in 2012
mash$ios.gattwto[mash$sftgcode=="SRB"] <- 0  # fill in missing
mash$ios.gattwto[mash$sftgcode=="SSD"] <- 0  # fill in missing
mash$ios.gattwto[mash$sftgcode=="USA"] <- 1  # fill in missing
mash$ios.gattwto[mash$sftgcode=="LAO" & mash$year>=2013] <- 1  # acceded in 2013
mash$ios.gattwto[mash$sftgcode=="TAJ" & mash$year>=2013] <- 1  # acceded in 2013

# Election years
# Source: IFES (http://www.electionguide.org/elections/?inst=&cont=&yr=2012)
# Source: NDI (http://www.ndi.org/2012calendar)
# Source: Wikipedia (http://en.wikipedia.org/wiki/List_of_elections_in_2012)
mash$nld.any.1 <- ifelse(mash$predyr <= 2011, mash$nld.any.1,
    ifelse((mash$sftgcode=="EGY" & mash$predyr==2012) |
           (mash$sftgcode=="FIN" & mash$predyr==2012) |
           (mash$sftgcode=="KZK" & mash$predyr==2012) |
           (mash$sftgcode=="KUW" & mash$predyr==2012) |
           (mash$sftgcode=="TKM" & mash$predyr==2012) |
           (mash$sftgcode=="YEM" & mash$predyr==2012) |
           (mash$sftgcode=="SEN" & mash$predyr==2012) |
           (mash$sftgcode=="IRN" & mash$predyr==2012) |
           (mash$sftgcode=="RUS" & mash$predyr==2012) |
           (mash$sftgcode=="SLO" & mash$predyr==2012) |
           (mash$sftgcode=="SAL" & mash$predyr==2012) |
           (mash$sftgcode=="GNB" & mash$predyr==2012) |  # Delayed, but planned in 2011
           (mash$sftgcode=="MYA" & mash$predyr==2012) |
           (mash$sftgcode=="ETM" & mash$predyr==2012) |
           (mash$sftgcode=="ARM" & mash$predyr==2012) |
           (mash$sftgcode=="GAM" & mash$predyr==2012) |
           (mash$sftgcode=="ROK" & mash$predyr==2012) |
           (mash$sftgcode=="FRN" & mash$predyr==2012) |
           (mash$sftgcode=="GRC" & mash$predyr==2012) |
           (mash$sftgcode=="SRB" & mash$predyr==2012) |
           (mash$sftgcode=="BHM" & mash$predyr==2012) |
           (mash$sftgcode=="SYR" & mash$predyr==2012) |
           (mash$sftgcode=="ALG" & mash$predyr==2012) |
           (mash$sftgcode=="DOM" & mash$predyr==2012) |
           (mash$sftgcode=="LES" & mash$predyr==2012) |
           (mash$sftgcode=="PNG" & mash$predyr==2012) |
           (mash$sftgcode=="MON" & mash$predyr==2012) |
           (mash$sftgcode=="MEX" & mash$predyr==2012) |
           (mash$sftgcode=="LIB" & mash$predyr==2012) |
           (mash$sftgcode=="CON" & mash$predyr==2012) |
           (mash$sftgcode=="ANG" & mash$predyr==2012) |
           (mash$sftgcode=="SOM" & mash$predyr==2012) |
           (mash$sftgcode=="NTH" & mash$predyr==2012) |
           (mash$sftgcode=="BLR" & mash$predyr==2012) |
           (mash$sftgcode=="GRG" & mash$predyr==2012) |
           (mash$sftgcode=="CZR" & mash$predyr==2012) |
           (mash$sftgcode=="LIT" & mash$predyr==2012) |
           (mash$sftgcode=="MNE" & mash$predyr==2012) |
           (mash$sftgcode=="UKR" & mash$predyr==2012) |
           (mash$sftgcode=="USA" & mash$predyr==2012) |
           (mash$sftgcode=="VEN" & mash$predyr==2012) |
           (mash$sftgcode=="SLV" & mash$predyr==2012) |
           (mash$sftgcode=="SIE" & mash$predyr==2012) |
           (mash$sftgcode=="BFO" & mash$predyr==2012) |
           (mash$sftgcode=="GHA" & mash$predyr==2012) |
           (mash$sftgcode=="RUM" & mash$predyr==2012) |
           (mash$sftgcode=="JPN" & mash$predyr==2012) |
           (mash$sftgcode=="ARM" & mash$predyr==2013) |
           (mash$sftgcode=="BNG" & mash$predyr==2013) |
           (mash$sftgcode=="IRN" & mash$predyr==2013) |
           (mash$sftgcode=="ISR" & mash$predyr==2013) |
           (mash$sftgcode=="KUW" & mash$predyr==2013) |
           (mash$sftgcode=="MON" & mash$predyr==2013) |
           (mash$sftgcode=="NEP" & mash$predyr==2013) |
           (mash$sftgcode=="QAT" & mash$predyr==2013) |
           (mash$sftgcode=="JOR" & mash$predyr==2013) |
           (mash$sftgcode=="PAK" & mash$predyr==2013) |
           (mash$sftgcode=="CAM" & mash$predyr==2013) | 
           (mash$sftgcode=="AZE" & mash$predyr==2013) |
           (mash$sftgcode=="TAJ" & mash$predyr==2013) | 
           (mash$sftgcode=="TKM" & mash$predyr==2013) | 
           (mash$sftgcode=="BHU" & mash$predyr==2013) |
           (mash$sftgcode=="MAL" & mash$predyr==2013) |
           (mash$sftgcode=="SIN" & mash$predyr==2013) |
           (mash$sftgcode=="JPN" & mash$predyr==2013) |
           (mash$sftgcode=="PHI" & mash$predyr==2013) |
           (mash$sftgcode=="DJI" & mash$predyr==2013) |
           (mash$sftgcode=="KEN" & mash$predyr==2013) |
           (mash$sftgcode=="MAG" & mash$predyr==2013) | 
           (mash$sftgcode=="GUI" & mash$predyr==2013) |
           (mash$sftgcode=="GNB" & mash$predyr==2013) |
           (mash$sftgcode=="TUN" & mash$predyr==2013) |
           (mash$sftgcode=="CAO" & mash$predyr==2013) |
           (mash$sftgcode=="ZAI" & mash$predyr==2013) |
           (mash$sftgcode=="ZIM" & mash$predyr==2013) |
           (mash$sftgcode=="MLI" & mash$predyr==2013) |
           (mash$sftgcode=="TOG" & mash$predyr==2013) |
           (mash$sftgcode=="RWA" & mash$predyr==2013) |
           (mash$sftgcode=="SWA" & mash$predyr==2013) | 
           (mash$sftgcode=="EQG" & mash$predyr==2013) |
           (mash$sftgcode=="AUS" & mash$predyr==2013) |
           (mash$sftgcode=="CZE" & mash$predyr==2013) |
           (mash$sftgcode=="GER" & mash$predyr==2013) |
           (mash$sftgcode=="BUL" & mash$predyr==2013) | 
           (mash$sftgcode=="ITA" & mash$predyr==2013) |
           (mash$sftgcode=="NOR" & mash$predyr==2013) |
           (mash$sftgcode=="CYP" & mash$predyr==2013) |
           (mash$sftgcode=="MNE" & mash$predyr==2013) | 
           (mash$sftgcode=="ALB" & mash$predyr==2013) |
           (mash$sftgcode=="GRG" & mash$predyr==2013) | 
           (mash$sftgcode=="POR" & mash$predyr==2013) |
           (mash$sftgcode=="HON" & mash$predyr==2013) |
           (mash$sftgcode=="CUB" & mash$predyr==2013) |
           (mash$sftgcode=="ARG" & mash$predyr==2013) |
           (mash$sftgcode=="CHL" & mash$predyr==2013) |
           (mash$sftgcode=="PAR" & mash$predyr==2013) |
           (mash$sftgcode=="VEN" & mash$predyr==2013) |
           (mash$sftgcode=="AUL" & mash$predyr==2013),
           1, ifelse(mash$year==2014, NA, 0) ) )

#########################################
# Modeling
#########################################

# Set 2010 as the year to split into training and test sets
splityear <- 2010

# Split into training and test sets
train <- subset(mash, year < splityear)
test <- subset(mash, year >= splityear)

# Estimate models in training set and generate predicted probabilities in test set
gripes.mod <- glm(gripes.f, family = binomial, data = train, na.action = na.exclude)
test$gripes.p <- predict(gripes.mod, newdata = test, type = "response")
modnzn.mod <- glm(modnzn.f, family = binomial, data = train, na.action = na.exclude)
test$modnzn.p <- predict(modnzn.mod, newdata = test, type = "response")
rscmob.mod <- glm(rscmob.f, family = binomial, data = train, na.action = na.exclude)
test$rscmob.p <- predict(rscmob.mod, newdata = test, type = "response")
polopp.mod <- glm(polopp.f, family = binomial, data = train, na.action = na.exclude)
test$polopp.p <- predict(polopp.mod, newdata = test, type = "response")
culled.mod <- glm(culled.f, family = binomial, data = train, na.action = na.exclude)
test$culled.p <- predict(culled.mod, newdata = test, type = "response")
base.mod <- glm(nvc.start.1 ~ log(wdi.pop), family = binomial, data = train, na.action = na.exclude)
test$base.p <- predict(base.mod, newdata = test, type = "response")
test$mean3.p <- 1/3 * (test$gripes.p + test$rscmob.p + test$polopp.p)

# Tidy up and write out the results
testsub <- subset(test, predyr < 2014)
testsub <- testsub[order(testsub$year, -testsub$culled.p),]
write.csv(testsub, "data/nvc.predictions.csv", row.names = FALSE)

# AUC scores by model
test2 <- subset(test, year < 2013)  # Limit comparison to years with preds for all models
base <- roc.area(test2$nvc.start.1, test2$base.p)
gripes <- roc.area(test2$nvc.start.1, test2$gripes.p)
modnzn <- roc.area(test2$nvc.start.1, test2$modnzn.p)
rscmob <- roc.area(test2$nvc.start.1, test2$rscmob.p)
polopp <- roc.area(test2$nvc.start.1, test2$polopp.p)
culled <- roc.area(test2$nvc.start.1, test2$culled.p)
mean3 <- roc.area(test2$nvc.start.1, test2$mean3.p)
aucs <- c(base$A, gripes$A, modnzn$A, rscmob$A, polopp$A, culled$A, mean3$A )
names(aucs) <- c("base", "gripes", "modnzn", "rscmob", "polopp", "culled", "mean")
aucs

# Brier scores by model
bss <- c(mean((test2$nvc.start.1 - test2$base.p)^2 , na.rm = TRUE),
  mean((test2$nvc.start.1 - test2$gripes.p)^2 , na.rm = TRUE),
  mean((test2$nvc.start.1 - test2$modnzn.p)^2 , na.rm = TRUE),
  mean((test2$nvc.start.1 - test2$rscmob.p)^2 , na.rm = TRUE),
  mean((test2$nvc.start.1 - test2$polopp.p)^2 , na.rm = TRUE),
  mean((test2$nvc.start.1 - test2$culled.p)^2 , na.rm = TRUE),
  mean((test2$nvc.start.1 - test2$mean3.p)^2 , na.rm = TRUE))
names(bss) <- c("base", "gripes", "modnzn", "rscmob", "polopp", "culled", "mean")
bss

# Logarithmic scores by model
logs <- c(mean((test2$nvc.start.1 * log(test2$base.p)) + ((1 - test2$nvc.start.1) * log(1 - test2$base.p)), na.rm = TRUE),
  mean((test2$nvc.start.1 * log(test2$gripes.p)) + ((1 - test2$nvc.start.1) * log(1 - test2$gripes.p)), na.rm = TRUE),
  mean((test2$nvc.start.1 * log(test2$modnzn.p)) + ((1 - test2$nvc.start.1) * log(1 - test2$modnzn.p)), na.rm = TRUE),
  mean((test2$nvc.start.1 * log(test2$rscmob.p)) + ((1 - test2$nvc.start.1) * log(1 - test2$rscmob.p)), na.rm = TRUE),
  mean((test2$nvc.start.1 * log(test2$polopp.p)) + ((1 - test2$nvc.start.1) * log(1 - test2$polopp.p)), na.rm = TRUE),
  mean((test2$nvc.start.1 * log(test2$culled.p)) + ((1 - test2$nvc.start.1) * log(1 - test2$culled.p)), na.rm = TRUE),
  mean((test2$nvc.start.1 * log(test2$mean3.p)) + ((1 - test2$nvc.start.1) * log(1 - test2$mean3.p)), na.rm = TRUE))
names(logs) <- c("base", "gripes", "modnzn", "rscmob", "polopp", "culled", "mean")
logs

all <- rbind(aucs, bss, logs)
all

# Bar chart of AUC by model
png(file = "figs/nvc.prediction.barplot.png",
     width=6, height=5, units='in', bg='white', res=150)
par(mai = c(0.5,0.5,0.25,0.25))
barplot(c(base$A - 0.5, gripes$A - 0.5, modnzn$A - 0.5, rscmob$A - 0.5, polopp$A - 0.5, culled$A - 0.5), ylim = c(0,0.25),
        names.arg=c("pop only", "grievance", "modnzn", "rsc mob", "pol opp", "culled"), cex.names = 0.75,
        col=c("black", "red", "dodgerblue", "forestgreen", "goldenrod", "gray75"), border = NA,
        axes = FALSE, main = "", xlab = "", ylab = "Area under ROC Curve" )
axis(2, at=c(0,0.1,0.2,0.3,0.4,0.5), labels = c("0.5","0.6","0.7","0.8","0.9","1"), tick = FALSE, las = 2)
dev.off()

# Bar chart of AUC by model by year
test.11 <- subset(test, predyr == 2011)
cull.11 <- round(roc.area(test.11$nvc.start.1, test.11$culled.p)$A, 4)
base.11 <- round(roc.area(test.11$nvc.start.1, test.11$base.p)$A, 4)
gripes.11 <- round(roc.area(test.11$nvc.start.1, test.11$gripes.p)$A, 4)
modnzn.11 <- round(roc.area(test.11$nvc.start.1, test.11$modnzn.p)$A, 4)
rscmob.11 <- round(roc.area(test.11$nvc.start.1, test.11$rscmob.p)$A, 4)
polopp.11 <- round(roc.area(test.11$nvc.start.1, test.11$polopp.p)$A, 4)
test.12 <- subset(test, predyr == 2012)
cull.12 <- round(roc.area(test.12$nvc.start.1, test.12$culled.p)$A, 4)
base.12 <- round(roc.area(test.12$nvc.start.1, test.12$base.p)$A, 4)
gripes.12 <- round(roc.area(test.12$nvc.start.1, test.12$gripes.p)$A, 4)
modnzn.12 <- round(roc.area(test.12$nvc.start.1, test.12$modnzn.p)$A, 4)
rscmob.12 <- round(roc.area(test.12$nvc.start.1, test.12$rscmob.p)$A, 4)
polopp.12 <- round(roc.area(test.12$nvc.start.1, test.12$polopp.p)$A, 4)

png(file = "figs/nvc.prediction.by.year.barplot.png",
     width=7.5, height=5, units='in', bg='white', res=150)
par(mai = c(0.75,0.5,0.25,0.25))
par(cex.axis = 0.75)
barplot(c(base.11, base.12,
          gripes.11, gripes.12,
          modnzn.11, modnzn.12,
          rscmob.11, rscmob.12,
          polopp.11, polopp.12,
          cull.11, cull.12),
  col=c("black", "black", "red", "red", "dodgerblue", "dodgerblue", "forestgreen", "forestgreen",
    "goldenrod", "goldenrod", "gray75", "gray75"), border = NA,
  ylim = c(0,1),
  names.arg=c("2011", "2012", "2011", "2012", "2011", "2012", "2011", "2012", "2011", "2012", "2011", "2012"), cex.names = 0.75,
  axes = FALSE, main = "", xlab = "", ylab = "Area under ROC Curve" )
axis(2, at=c(0,0.25,0.5,0.75,1), labels = c("0","0.25","0.50","0.75","1"), tick = FALSE, las = 2, pos = 0.1)
mtext(side = 1, at = 1.25, line = 2.5, "pop only", cex = 0.8)
mtext(side = 1, at = 3.75, line = 2.5, "grievance", cex = 0.8)
mtext(side = 1, at = 6.1, line = 2.5, "modernization", cex = 0.8)
mtext(side = 1, at = 8.55, line = 2.5, "resource mob", cex = 0.8)
mtext(side = 1, at = 11, line = 2.5, "pol opp", cex = 0.8)
mtext(side = 1, at = 13.25, line = 2.5, "culled", cex = 0.8)
dev.off()

###############################################
# Dot plots by year
###############################################

## CULLED MODEL ##

# 2011
pred11 <- subset(test, predyr==2011, select=c(country, nvc.start.1, culled.p))
pred11 <- pred11[order(-pred11$culled.p),]
row.names(pred11) <- NULL
pred11$rank <- as.numeric(as.character(row.names(pred11)))
pred11$country <- as.character(pred11$country)
auc11 <- roc.area(pred11$nvc.start.1, pred11$culled.p)
condcol <- ifelse(pred11$nvc.start.1==1, "deeppink1", "gray")
png(file = "figs/nvc.prediction.dotplot.2011.culled.png",
    width=14, height=18, units='cm', bg='white', res=150)
dotchart2(pred11$culled.p[1:40], labels=pred11$country[1:40],
    lines=TRUE, lwd=0.05, lty=3, sort=FALSE, dotsize=1.25, pch=20, col=condcol, cex.labels=0.75, xlim=c(0,0.75))
title(main=list("Probability of Campaign Onset in 2011 (Top 40)", cex=1),
    sub = list(paste("AUC", as.character(substr(auc11$A,1,4)), sep = " = "), cex=0.8) )
dev.off()

# 2012
pred12 <- subset(test, predyr==2012, select=c(country, nvc.start.1, culled.p))
pred12 <- pred12[order(-pred12$culled.p),]
row.names(pred11) <- NULL
pred12$rank <- as.numeric(as.character(row.names(pred12)))
pred12$country <- as.character(pred12$country)
auc12 <- roc.area(pred12$nvc.start.1, pred12$culled.p)
condcol <- ifelse(pred12$nvc.start.1==1, "deeppink1", "gray")
png(file = "figs/nvc.prediction.dotplot.2012.culled.png",
     width=14, height=18, units='cm', bg='white', res=150)
dotchart2(pred12$culled.p[1:40], labels=pred12$country[1:40],
     lines=TRUE, lwd=0.05, lty=3, sort=FALSE, dotsize=1.25, pch=20, col=condcol, cex.labels=0.75, xlim=c(0,0.75))
title(main=list("Probability of Campaign Onset in 2012 (Top 40)", cex=1),
    sub = list(paste("AUC", as.character(substr(auc12$A,1,4)), sep = " = "), cex=0.8) )
dev.off()


# 2013
pred13 <- subset(test, predyr==2013, select=c(country, nvc.start.1, culled.p))
pred13 <- pred13[order(-pred13$culled.p),]
row.names(pred11) <- NULL
pred13$rank <- as.numeric(as.character(row.names(pred13)))
pred13$country <- as.character(pred13$country)
auc13 <- roc.area(pred13$nvc.start.1, pred13$culled.p)
condcol <- ifelse(pred13$nvc.start.1==1, "deeppink1", "gray")
png(file = "figs/nvc.prediction.dotplot.2013.culled.png",
     width=14, height=18, units='cm', bg='white', res=150)
dotchart2(pred13$culled.p[1:40], labels=pred13$country[1:40],
     lines=TRUE, lwd=0.05, lty=3, sort=FALSE, dotsize=1.25, pch=20, col=condcol, cex.labels=0.75, xlim=c(0,0.75))
title(main=list("Probability of Campaign Onset in 2013 (Top 40)", cex=1),
    sub = list(paste("AUC", as.character(substr(auc13$A,1,4)), sep = " = "), cex=0.8) )
dev.off()

## GRIEVANCES ##

# 2011
pred11 <- subset(test, predyr==2011, select=c(country, nvc.start.1, gripes.p))
pred11 <- pred11[order(-pred11$gripes.p),]
row.names(pred11) <- NULL
pred11$rank <- as.numeric(as.character(row.names(pred11)))
pred11$country <- as.character(pred11$country)
auc11 <- roc.area(pred11$nvc.start.1, pred11$gripes.p)
condcol <- ifelse(pred11$nvc.start.1==1, "deeppink1", "gray")
png(file = "figs/nvc.prediction.dotplot.2011.gripes.png",
    width=14, height=18, units='cm', bg='white', res=150)
dotchart2(pred11$gripes.p[1:40], labels=pred11$country[1:40],
    lines=TRUE, lwd=0.05, lty=3, sort=FALSE, dotsize=1.25, pch=20, col=condcol, cex.labels=0.75, xlim=c(0,0.75))
title(main=list("Probability of Campaign Onset in 2011 (Top 40)", cex=1),
    sub = list(paste("AUC", as.character(substr(auc11$A,1,4)), sep = " = "), cex=0.8) )
dev.off()

# 2012
pred12 <- subset(test, predyr==2012, select=c(country, nvc.start.1, gripes.p))
pred12 <- pred12[order(-pred12$gripes.p),]
row.names(pred11) <- NULL
pred12$rank <- as.numeric(as.character(row.names(pred12)))
pred12$country <- as.character(pred12$country)
auc12 <- roc.area(pred12$nvc.start.1, pred12$gripes.p)
condcol <- ifelse(pred12$nvc.start.1==1, "deeppink1", "gray")
png(file = "figs/nvc.prediction.dotplot.2012.gripes.png",
     width=14, height=18, units='cm', bg='white', res=150)
dotchart2(pred12$gripes.p[1:40], labels=pred12$country[1:40],
     lines=TRUE, lwd=0.05, lty=3, sort=FALSE, dotsize=1.25, pch=20, col=condcol, cex.labels=0.75, xlim=c(0,0.75))
title(main=list("Probability of Campaign Onset in 2012 (Top 40)", cex=1),
    sub = list(paste("AUC", as.character(substr(auc12$A,1,4)), sep = " = "), cex=0.8) )
dev.off()


# 2013
pred13 <- subset(test, predyr==2013, select=c(country, nvc.start.1, gripes.p))
pred13 <- pred13[order(-pred13$gripes.p),]
row.names(pred11) <- NULL
pred13$rank <- as.numeric(as.character(row.names(pred13)))
pred13$country <- as.character(pred13$country)
auc13 <- roc.area(pred13$nvc.start.1, pred13$gripes.p)
condcol <- ifelse(pred13$nvc.start.1==1, "deeppink1", "gray")
png(file = "figs/nvc.prediction.dotplot.2013.gripes.png",
     width=14, height=18, units='cm', bg='white', res=150)
dotchart2(pred13$gripes.p[1:40], labels=pred13$country[1:40],
     lines=TRUE, lwd=0.05, lty=3, sort=FALSE, dotsize=1.25, pch=20, col=condcol, cex.labels=0.75, xlim=c(0,0.75))
title(main=list("Probability of Campaign Onset in 2013 (Top 40)", cex=1),
    sub = list(paste("AUC", as.character(substr(auc13$A,1,4)), sep = " = "), cex=0.8) )
dev.off()

## POLITICAL OPPORTUNITY ##

# 2011
pred11 <- subset(test, predyr==2011, select=c(country, nvc.start.1, polopp.p))
pred11 <- pred11[order(-pred11$polopp.p),]
row.names(pred11) <- NULL
pred11$rank <- as.numeric(as.character(row.names(pred11)))
pred11$country <- as.character(pred11$country)
auc11 <- roc.area(pred11$nvc.start.1, pred11$polopp.p)
condcol <- ifelse(pred11$nvc.start.1==1, "deeppink1", "gray")
png(file = "figs/nvc.prediction.dotplot.2011.polopp.png",
    width=14, height=18, units='cm', bg='white', res=150)
dotchart2(pred11$polopp.p[1:40], labels=pred11$country[1:40],
    lines=TRUE, lwd=0.05, lty=3, sort=FALSE, dotsize=1.25, pch=20, col=condcol, cex.labels=0.75, xlim=c(0,0.75))
title(main=list("Probability of Campaign Onset in 2011 (Top 40)", cex=1),
    sub = list(paste("AUC", as.character(substr(auc11$A,1,4)), sep = " = "), cex=0.8) )
dev.off()

# 2012
pred12 <- subset(test, predyr==2012, select=c(country, nvc.start.1, polopp.p))
pred12 <- pred12[order(-pred12$polopp.p),]
row.names(pred11) <- NULL
pred12$rank <- as.numeric(as.character(row.names(pred12)))
pred12$country <- as.character(pred12$country)
auc12 <- roc.area(pred12$nvc.start.1, pred12$polopp.p)
condcol <- ifelse(pred12$nvc.start.1==1, "deeppink1", "gray")
png(file = "figs/nvc.prediction.dotplot.2012.polopp.png",
     width=14, height=18, units='cm', bg='white', res=150)
dotchart2(pred12$polopp.p[1:40], labels=pred12$country[1:40],
     lines=TRUE, lwd=0.05, lty=3, sort=FALSE, dotsize=1.25, pch=20, col=condcol, cex.labels=0.75, xlim=c(0,0.75))
title(main=list("Probability of Campaign Onset in 2012 (Top 40)", cex=1),
    sub = list(paste("AUC", as.character(substr(auc12$A,1,4)), sep = " = "), cex=0.8) )
dev.off()

# 2013
pred13 <- subset(test, predyr==2013, select=c(country, nvc.start.1, polopp.p))
pred13 <- pred13[order(-pred13$polopp.p),]
row.names(pred11) <- NULL
pred13$rank <- as.numeric(as.character(row.names(pred13)))
pred13$country <- as.character(pred13$country)
auc13 <- roc.area(pred13$nvc.start.1, pred13$polopp.p)
condcol <- ifelse(pred13$nvc.start.1==1, "deeppink1", "gray")
png(file = "figs/nvc.prediction.dotplot.2013.polopp.png",
     width=14, height=18, units='cm', bg='white', res=150)
dotchart2(pred13$polopp.p[1:40], labels=pred13$country[1:40],
     lines=TRUE, lwd=0.05, lty=3, sort=FALSE, dotsize=1.25, pch=20, col=condcol, cex.labels=0.75, xlim=c(0,0.75))
title(main=list("Probability of Campaign Onset in 2013 (Top 40)", cex=1),
    sub = list(paste("AUC", as.character(substr(auc13$A,1,4)), sep = " = "), cex=0.8) )
dev.off()
