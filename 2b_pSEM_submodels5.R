################################################################################### pSEM submodels #####
##### Construct the individual models that could make up the global SEM #####
##### Submodel including expansion only #####
##### Written by Regan Early #####
##### Written on: 19th August 2021 #####
##### Modified by: Regan Early #####
##### Modified on: 22nd August 2022 (checking results and ready for github) #####
##############################################################################

.libPaths("C:/SOFTWARE/R/R-4.1.2/library")
library(ggplot2)
library(jtools) ## effect_plot() - plots partial residuals in multiple regressions - and summ() - super neat model summary. Needs latest version of ggplot2
library(ggfortify) ## nice diagnostic plots
library(MuMIn) ## helpful if want to identify best model subsets that respect marginality of quadratic relationships
library(car) ## VIF
library(betareg)
        
setwd("E:\\NON_PROJECT\\TEACHING\\CORNWALL\\UG_PROJECTS\\2014-15\\JACK_HARNESS\\PATH_ANALYSIS\\2022")

dat <- read.csv("regan_dataset_analyse2_transformed.csv")

### Use 'woody' or not rather than all three growth forms. 
## There are only 9 shrubs, so little point coding them distinctly. 
## Having a binary variable makes it easier to examine interactions in SEM. 
## See powerpoint for previous analyses with shrub and tree.
## Results not meaningfully changed by replacing three growth forms with two.
dat$woody <- 0
dat$woody[dat$growth %in% c("Shrub","Tree")] <- 1

##### Submodel 5 #####
par(mfrow=c(2,3))
par(mar=c(4, 4, 1, 1)) # c(bottom, left, top, right)
plot(expansion_t ~ rsEURm2, data=dat)
plot(expansion_t ~ rsUSAm2, data=dat)
plot(expansion_t ~ yr1, data=dat)
plot(expansion_t ~ en_rel_t, data=dat)
plot(expansion_t ~ growth, data=dat)

### Linear model
summary(m5 <- lm(expansion_t ~ rsEURm2_t + rsUSAm2_t + yr1 + en_rel_t + woody, data=dat)) ## rsUSA and yr1 NOT significant
autoplot(m5, label.size = 1) ## Bit suspicious - residuals vs fitted line stillcurved

summary(m5.logit <- lm(logit(expansion_t) ~ rsEURm2_t + rsUSAm2_t + yr1 + en_rel_t + woody, data=dat)) ## yr1 NOT significant
autoplot(m5.logit, label.size = 1) ## Bit suspicious. Residuals vs. fitted looking much better.

summary(m5.binomial <- glm(expansion_t ~ rsEURm2_t + rsUSAm2_t + yr1 + en_rel_t + woody, data=dat, family="binomial")) ## Nothing except European range size significant - very surprising
autoplot(m5.binomial, label.size = 1) ## Residuals vs fitted line still curved. Beta regression might be better.

summary(m5.beta <- betareg(expansion_t ~ rsEURm2_t + rsUSAm2_t + yr1 + en_rel_t + woody, data=dat)) ## yr1 not significant
par(mfrow=c(2,3)); par(mar=c(3,3,3,2)); plot(m5.beta, which=c(1:2,4:6), cex.lab=0.5, cex.main=0.5)
vif(m5.beta) ## all fine

### logit transformed lm actually performs well
par(mfrow=c(2,2))
par(mar=c(3.8,3.8,3,1))
plot(m5$residuals ~ m5$fitted.values, main="lm")
plot(m5.logit$residuals ~ m5.logit$fitted.values, main="logit trans")
plot(m5.binomial$residuals ~ m5.binomial$fitted.values, main="binomial")
plot(m5.beta$residuals ~ m5.beta$fitted.values, main="betareg")

## Quadratic model with logit transformation
## Using q and I(q^2) is problematic because they will be correlated and correlated variables can cause problems. 
## The use of poly() avoids this by producing orthogonal polynomials
## However, in order to test for interactions between linear effects and to dredge, have to create new variables for first and second order polynomial
dat$rsEURm2_t1 <- poly(dat$rsEURm2_t, 2)[,1]
dat$rsEURm2_t2 <- poly(dat$rsEURm2_t, 2)[,2]
dat$rsUSAm2_t1 <- poly(dat$rsUSAm2_t, 2)[,1]
dat$rsUSAm2_t2 <- poly(dat$rsUSAm2_t, 2)[,2]
dat$yr1_1 <- poly(dat$yr1, 2)[,1]
dat$yr1_2 <- poly(dat$yr1, 2)[,2]
dat$en_rel_t1 <- poly(dat$en_rel_t, 2)[,1]
dat$en_rel_t2 <- poly(dat$en_rel_t, 2)[,2]

summary(m5q <- lm(logit(expansion_t) ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + en_rel_t1 + en_rel_t2 + woody, data=dat, na.action=na.fail)) ## Only range sizes significant, yr and enemy release not. Growth *was* sig when coded as three forms
par(mfrow=c(2,3)); par(mar=c(3,3,3,2)); plot(m5q, which=1:6) ## Seems ok
vif(m5q) ## ok

## Stepwise quadratic model
m5s <- dredge(m5q)
m5s[1] ## best model includes en_rel_t, woody, European (linear and quadratic) and USA range size.
m5sb <- get.models(m5s, subset=1)$'94'
summ(m5sb) ## adj R2=0.49
AIC(m5sb) ## 445.3154

### Test best model with betareg
summary(m5q.betareg <- betareg(expansion_t ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + en_rel_t1 + en_rel_t2 + woody, data=dat, na.action=na.fail)) ## Only Europe range sizes USA range size (linear) and growth significant, yr and enemy release not. 

## Stepwise quadratic model (need to use I(...^2) format in order to remove quadratic terms seperately from linear terms
m5s.betareg <- dredge(m5q.betareg)
m5s.betareg[1] ## best model includes en_rel_t, woody, European (linear and quadratic) and USA range size 
m5s.betaregb <- get.models(m5s.betareg, subset=1)$'94'
AIC(m5s.betaregb) ## -271.1076
## So classical betaregression supports the best model selected by logit transforming the response variable.

##### Interactions #####
summary(m5i <- lm(logit(expansion_t) ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + woody +
                       rsEURm2_t1*rsUSAm2_t1 + rsEURm2_t1*en_rel_t1 + rsUSAm2_t1*en_rel_t1 + 
                       woody*rsEURm2_t1 + woody*rsUSAm2_t1 + woody*en_rel_t1, 
                     data=dat, na.action=na.fail)) ## 
m5i.d <- dredge(m5i) 
m5i.d[1] ## AIC best model 443.6. Includes interaction woody*rsEUR.. en_rel_t retained in six of seven models in best model subset
m5i.db <- get.models(m5i.d, subset=1)$'544'
summ(m5i.db) ## adj R2=0.51
AIC(m5i.db) ## 442.3803

vif(m5i.db) ## good
par(mfrow=c(2,3)); par(mar=c(3,3,3,2)); plot(m5i.db, which=1:6) ## Seems ok

### Test interactions in the best model selected with logit, using betareg
summary(m5i.beta <- betareg(expansion_t ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + woody +
                              rsEURm2_t1*rsUSAm2_t1 + rsEURm2_t1*en_rel_t1 + rsUSAm2_t1*en_rel_t1 + 
                              woody*rsEURm2_t1 + woody*rsUSAm2_t1 + woody*en_rel_t1, 
                               data=dat, na.action=na.fail)) ## 
dredge(m5i.beta) ## best model is the same as the logit version AIC=-272.7