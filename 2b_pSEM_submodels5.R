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
        
setwd("E:\\NON_PROJECT\\TEACHING\\CORNWALL\\UG_PROJECTS\\2014-15\\JACK_HARNESS\\PATH_ANALYSIS\\2021")

dat <- read.csv("regan_dataset_analyse2_transformed.csv")

##### Submodel 5 #####
par(mfrow=c(2,3))
par(mar=c(4, 4, 1, 1)) # c(bottom, left, top, right)
plot(expansion_t ~ rsEURm2, data=dat)
plot(expansion_t ~ rsUSAm2, data=dat)
plot(expansion_t ~ yr1, data=dat)
plot(expansion_t ~ en_rel_t, data=dat)
plot(expansion_t ~ growth, data=dat)

### Linear model
summary(m5 <- lm(expansion_t ~ rsEURm2_t + rsUSAm2_t + yr1 + en_rel_t + growth, data=dat)) ## rsUSA and yr1 NOT significant
autoplot(m5, label.size = 1) ## Bit suspicious - residuals vs fitted line stillcurved

summary(m5.logit <- lm(logit(expansion_t) ~ rsEURm2_t + rsUSAm2_t + yr1 + en_rel_t + growth, data=dat)) ## yr1 NOT significant
autoplot(m5.logit, label.size = 1) ## Bit suspicious. Residuals vs. fitted looking much better.

summary(m5.binomial <- glm(expansion_t ~ rsEURm2_t + rsUSAm2_t + yr1 + en_rel_t + growth, data=dat, family="binomial")) ## Nothing except European range size significant - very surprising
autoplot(m5.binomial, label.size = 1) ## Residuals vs fitted line still curved. Beta regression might be better.

summary(m5.beta <- betareg(expansion_t ~ rsEURm2_t + rsUSAm2_t + yr1 + en_rel_t + growth, data=dat)) ## yr1 not significant
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
# summary(m5q <- lm(logit(expansion_t) ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail)) ## Only range sizes and growth significant, yr and enemy release not. 
summary(m5q <- lm(logit(expansion_t) ~ rsEURm2_t + I(rsEURm2_t^2) + rsUSAm2_t + I(rsUSAm2_t^2) + yr1 + I(yr1^2) + en_rel_t + I(en_rel_t^2) + growth, data=dat, na.action=na.fail)) ## Only range sizes and growth significant, yr and enemy release not. 
par(mfrow=c(2,3)); par(mar=c(3,3,3,2)); plot(m5q, which=1:6) ## Seems ok
vif(m5q) ## inflation

## Stepwise quadratic model (need to use I(...^2) format in order to remove quadratic terms seperately from linear terms
m5s <- dredge(m5q)
m5s[1] ## best model includes en_rel_t, growth, European (linear and quadratic) and USA range size AIC = 445.4
m5sb <- get.models(m5s, subset=1)$'62'
summ(m5sb) ## adj R2=0.50
AIC(m5sb) ## 444.21

### Test best model with poly
AIC(m5sb.poly <- lm(logit(expansion_t) ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 1) + en_rel_t + growth, data=dat, na.action=na.fail)) ## 444.2102. 

### Test best model with betareg
# summary(m5q.betareg <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail)) ## ## Only range sizes and shrub significant, yr and enemy release not.
summary(m5q.betareg <- betareg(expansion_t ~ rsEURm2_t + I(rsEURm2_t^2) + rsUSAm2_t + I(rsUSAm2_t^2) + yr1 + I(yr1^2) + en_rel_t + I(en_rel_t^2) + growth, data=dat, na.action=na.fail)) ## Only Europe range sizes and shrub significant, yr and enemy release not. 

## Stepwise quadratic model (need to use I(...^2) format in order to remove quadratic terms seperately from linear terms
m5s.betareg <- dredge(m5q.betareg)
m5s.betareg[1] ## best model includes en_rel_t, growth, European (linear and quadratic) and USA range size AIC = -270.6
m5s.betaregb <- get.models(m5s.betareg, subset=1)$'62'
AIC(m5s.betaregb) ## -271.8066

### So classical betaregression supports the best model selected by logit transforming the response variable.

##### Interactions #####
summary(m5i <- lm(logit(expansion_t) ~ rsEURm2_t + I(rsEURm2_t^2) + rsUSAm2_t + en_rel_t + growth +
                       rsEURm2_t*rsUSAm2_t + rsEURm2_t*en_rel_t + rsUSAm2_t*en_rel_t + 
                       growth*rsEURm2_t + growth*rsUSAm2_t + growth*en_rel_t, 
                     data=dat, na.action=na.fail)) ## 449.8271. 
dredge(m5i) ## AIC best model 444.3. Includes interaction growth*rsEUR.
vif(m5i) ## inflation a little high - because of quadratics?

## Best model with poly:
AIC(m5i.poly <- lm(logit(expansion_t) ~ poly(rsEURm2_t, 2) + rsUSAm2_t + en_rel_t + growth +
                    growth*poly(rsEURm2_t, 1), 
                  data=dat, na.action=na.fail)) ## AIC 442.5015

# no interactions between en_rel_t and anything other than growth retained in best model subset 
# growth interacts with something (usually USA or EUR RS) in all but one of the ten best models
# the specific interaction with growth has very little effect on the coefficients of the other terms
# en_rel_t retained in seven of ten

***# Suggest making SEM without interactions and check whether adding them in makes a difference

par(mfrow=c(2,3)); par(mar=c(3,3,3,2)); plot(m5i, which=1:6) ## Seems ok


### Test interactions in the best model selected with logit, using betareg
summary(m5i.beta <- betareg(expansion_t ~ rsEURm2_t + I(rsEURm2_t^2) + rsUSAm2_t + en_rel_t + growth +
                                 rsEURm2_t*rsUSAm2_t + rsEURm2_t*en_rel_t + rsUSAm2_t*en_rel_t + 
                                 growth*rsEURm2_t + growth*rsUSAm2_t + growth*en_rel_t, 
                               data=dat, na.action=na.fail)) ## -262.7203
dredge(m5i.beta) 

# best model is the same as the logit version and includes interaction growth*rsEUR. AIC=-271.0 
summary(m5i.betab <- betareg(expansion_t ~ rsEURm2_t + I(rsEURm2_t^2) + rsUSAm2_t + en_rel_t + growth +
                                  growth*rsEURm2_t, 
                                data=dat, na.action=na.fail)) 
AIC(m5i.betab.poly <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + rsUSAm2_t + en_rel_t + growth, data=dat, na.action=na.fail)) ## -271.8066
