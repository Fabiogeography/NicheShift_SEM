################################################################################### pSEM submodels #####
##### Investigate logit model against classical and zero-inflated betaregression #####
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
library(lmtest) ## lrtest
library(withr)
library(brms) ## zero-augmented models
library(rstan) ## Note that installing this was problematic. I had to use Sys.getenv("R_LIBS_USER"); file.edit("~/.Renviron") and set R_LIBS_USER = D pathway
# To avoid recompilation of unchanged Stan programs, we recommend calling:
# rstan_options(auto_write = TRUE)
# Do not specify '-march=native' in 'LOCAL_CPPFLAGS' or a Makevars file
library(SimDesign) ## rmvnorm

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
summary(m5 <- betareg(expansion_t ~ rsEURm2_t + rsUSAm2_t + yr1 + en_rel_t + growth, data=dat)) ## yr1 not significant
par(mfrow=c(2,3)); par(mar=c(3,3,3,2)); plot(m5, which=c(1:2,4:6), cex.lab=0.5, cex.main=0.5)
vif(m5.beta) ## all fine

## Quadratic model with betaregression
summary(m5q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail)) ## ## Only range sizes significant, yr and enemy release not. 
par(mfrow=c(2,3)); par(mar=c(3,3,3,2)); plot(m5q, which=1:6) ## Seems ok
vif(m5q) ## marginal - en_rel_t is 5.003

## Is the big difference in numerical values of some explanatory variables causing a problem?
summary(m5.1q <- betareg(expansion_t ~ poly(scale(rsEURm2_t), 2) + poly(scale(rsUSAm2_t), 2) + poly(scale(yr1), 2) + poly(scale(en_rel_t), 2) + growth, data=dat, na.action=na.fail)) ## Making numerical values of explanatory variables more similar does not change significance.
par(mfrow=c(2,3)); par(mar=c(3,3,3,2)); plot(m5.1q, which=1:6) ## Seems almost identical to results with unscaled variables

## Are two species with very high expansion affecting results?
dat2 <- dat[!dat$Host%in%c("Nelumbo nucifera", "Sagina subulata"),]
summary(m5.2q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, data=dat2, na.action=na.fail)) ## ## Only range sizes significant, yr and enemy release not.
par(mfrow=c(2,3)); par(mar=c(3,3,3,2)); plot(m5.2q, which=1:6) ## Seems similar to results with these species

## Stepwise model
m5s <- dredge(m5q)
m5s[1] ## best model (but linear and quad terms included / excluded together). AIC = -269.4

## Manual model selection to investigate similar models to best model
AIC(m5q) ## -266.9716
AIC(m5.3q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + yr1 + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail)) ## -268.5119
AIC(m5.4q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail)) ## -268.5471
AIC(m5.5q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + en_rel_t + growth, data=dat, na.action=na.fail)) ## -269.9637
AIC(m5.6q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + growth, data=dat, na.action=na.fail)) ## -270.5813 Removing en_rel_t altogether is best? Why not found in dredge?
AIC(m5.7q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2), data=dat, na.action=na.fail)) ## -268.7519
AIC(m5.8q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + rsUSAm2_t, data=dat, na.action=na.fail)) ## -267.7897
AIC(m5.9q <- betareg(expansion_t ~ poly(rsEURm2_t, 2), data=dat, na.action=na.fail)) ## -236.9898
AIC(m5.10q <- betareg(expansion_t ~ rsEURm2_t + poly(rsUSAm2_t, 2), data=dat, na.action=na.fail)) ## -233.4178
AIC(m5.11q <- betareg(expansion_t ~ poly(rsUSAm2_t, 2), data=dat, na.action=na.fail)) ## -183.93
AIC(m5.12q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + growth, data=dat, na.action=na.fail)) ## -268.966
AIC(m5.13q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + yr1 + growth, data=dat, na.action=na.fail)) ## -270.5499
AIC(m5.14q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + yr1 + en_rel_t + growth, data=dat, na.action=na.fail)) ## -269.6209

AIC(m5.15q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2), data=dat, na.action=na.fail)) ## -264.573
AIC(m5.16q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + yr1 + poly(en_rel_t, 2), data=dat, na.action=na.fail)) ## -265.9659
AIC(m5.17q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(en_rel_t, 2), data=dat, na.action=na.fail)) ## -265.7859
AIC(m5.18q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + en_rel_t, data=dat, na.action=na.fail)) ## -267.3905

## AIC of best model is -269.4, so dAIC2 threshold is -267.4. Thus the best model subset always includes linear and quadratic fx of both range sizes, growth in all but one model, and includes models with linear and quadratic effects of yr1 and en_rel_t, but not with their quadratic effects together. poly(yr1,2) gives slightly better AIC than poly(en_rel_t, 2) in each other's absence.

# ## Try without quadratic effect of USA rs
# AIC(m5.30q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + rsUSAm2_t + yr1 + growth, data=dat, na.action=na.fail)) ## -270.7992
# AIC(m5.31q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + rsUSAm2_t + growth, data=dat, na.action=na.fail)) ## -270.8262
# AIC(m5.32q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + rsUSAm2_t + poly(yr1, 2) + growth, data=dat, na.action=na.fail)) ## -269.0386
# AIC(m5.33q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + rsUSAm2_t + en_rel_t + growth, data=dat, na.action=na.fail)) ## -271.8066
# AIC(m5.34q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + rsUSAm2_t + yr1 + en_rel_t + growth, data=dat, na.action=na.fail)) ## -271.3893
# AIC(m5.35q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + rsUSAm2_t + yr1 + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail)) ## -270.0947
# AIC(m5.36q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + rsUSAm2_t + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail)) ## -270.2743
# AIC(m5.37q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + rsUSAm2_t, data=dat, na.action=na.fail)) ## -267.7897

##### Interactions #####
### betareg doesn't seem able to test for interactions when poly() is used, so test with I(^2)
summary(m5.12i <- betareg(expansion_t ~ rsEURm2_t + I(rsEURm2_t^2) + rsUSAm2_t + I(rsUSAm2_t^2) + 
                            yr1 + I(yr1^2) + growth +
                            rsEURm2_t*rsUSAm2_t + rsEURm2_t*yr1 + rsUSAm2_t*yr1,
                          data=dat, na.action=na.fail)) ## AIC and p-values indicate no sig effect of interactions

summary(m5.4i <- betareg(expansion_t ~ rsEURm2_t + I(rsEURm2_t^2) + rsUSAm2_t + 
                           I(rsUSAm2_t^2) + en_rel_t + I(en_rel_t^2) + growth +
                           rsEURm2_t*rsUSAm2_t + rsEURm2_t*en_rel_t + rsUSAm2_t*en_rel_t,
                         data=dat, na.action=na.fail)) ## AIC and p-values indicate no sig effect of interactions


##### Try different links with full quadratic model #####
summary(m5.19q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail, link="probit")) ## ## Only range sizes significant, enemy release^2 near sig. 
par(mfrow=c(2,3)); par(mar=c(3,3,3,2)); plot(m5.19q, which=1:6) ## Very similar to logit link
vif(m5.19q) ## just ok

summary(m5.20q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail, link="cloglog")) ## ## Only range sizes significant
par(mfrow=c(2,3)); par(mar=c(3,3,3,2)); plot(m5.20q, which=1:6) ## Very similar to logit link
vif(m5.20q) ## en_rel_t now exceeds 5

summary(m5.21q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail, link="cauchit")) ## ## Range sizes and growth sig. yr approaching sig. 
par(mfrow=c(2,3)); par(mar=c(3,3,3,2)); plot(m5.21q, which=1:6) ## Very similar to logit link
vif(m5.21q) ## both range sizes and en_rel_t now exceed 5

summary(m5.22q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail, link="log")) ## fails 

summary(m5.23q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail, link="loglog")) ## ## Range sizes, growth en_rel_t sig. 
par(mfrow=c(2,3)); par(mar=c(3,3,3,2)); plot(m5.23q, which=1:6) ## Very similar to logit link
vif(m5.23q) ## All <4.06

## Compare logit and loglog results:
par(mfrow=c(1,2))
plot(m5.23q, which=4, main="loglog")
plot(m5q, which=4, main="logit")

## Compare AICs between models with different link functions. AICs are comparable between models with different link functions as long as the models are fit with ML (true) https://stats.stackexchange.com/questions/386275/is-the-use-of-loglik-or-aic-to-compare-logit-probit-cloglog-models-valid
AIC(m5q) ## -266.9716
AIC(m5.19q) ## -263.5092
AIC(m5.20q) ## -273.7979 cloglog
AIC(m5.21q) ## -275.072 cauchit
AIC(m5.23q) ## -254.3733. Investigate removing year: 
AIC(betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail, link="loglog")) ## -257.0977. Beneficial.

## Compare logit and loglog results with two high expansion points removed:
summary(m5.24q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, data=dat2, na.action=na.fail, link="loglog")) ## Range sizes, growth en_rel_t sig, which is good - same as model with the two points.
par(mfrow=c(2,3)); par(mar=c(3,3,3,2)); plot(m5.24q, which=1:6) ## Very similar to logit link, points a bit stretched out on the right
vif(m5.24q) ## just ok

par(mfrow=c(1,2))
plot(m5.24q, which=4, main="loglog, 2 pts rm")
plot(m5.2q, which=4, main="logit, 2 pts rm")

par(mfrow=c(1,2)); par(mar=c(5,3,3,2))
hist(dat$expansion_t, main="all data")
hist(dat2$expansion_t, main="two big expanders removed")

### See pages 15 and 16 of (https://cran.r-project.org/web/packages/betareg/vignettes/betareg.pdf) for deciding if log-log link function better
par(mfrow=c(1,2))
plot(abs(residuals(m5.23q, type = "response")),
     + abs(residuals(m5q, type = "response")), main="all data")
abline(0, 1, lty = 2)

plot(abs(residuals(m5.24q, type = "response")),
     + abs(residuals(m5.2q, type = "response")), main="two big expanders removed")
abline(0, 1, lty = 2)

plot(m5q, which = 6, xlim=c(0,0.5), main="logit"); plot(m5.23q, which = 6, xlim=c(0,0.5), main="loglog")

## Add powers of linear predictors to the model. In well-specified models this should not lead to an improvement.
lrtest(m5q, . ~ . + I(predict(m5q, type = "link")^2)) ## Does improve the model - mispecified?
lrtest(m5.23q, . ~ . + I(predict(m5.23q, type = "link")^2)) ## Does improve the model - mispecified?

## Compare log-likelihoods of the different link-functions
sapply(c("logit", "probit", "cloglog", "cauchit", "loglog"), function(x) logLik(update(m5q, link = x))) ## cauchit and cloglog perform best, as they did with AIC (obviously, since each model has teh same number of parameters). loglog performs worst

####### Alternatives #####
### Allow phi to vary with predictors. Is this akin to a random effect?
summary(m5.25q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth | growth, data=dat, na.action=na.fail, link="loglog")) ## Shrub might have a lower phi than other growth forms - just non-sig effect
lrtest(m5.23q, m5.25q) ## no significant difference between models 

### Test bias correction method
# Biased estimators can lead to erroneous inference, and methods for bias-reduction and bias-correction 
# are an active research area in statistics and in beta regression in particular (Grün et al., 2012; 
# Kosmidis, 2014; Kosmidis & Firth, 2010; Ospina, Cribari-Neto, & Vasconcellos, 2006; Simas et al., 2010).
# In particular, the precision parameter in beta regression models is prone to overestimation bias 
# (Kosmidis & Firth, 2010) which leads directly to underestimation of the width of confidence intervals 
# for other model parameters. Two main types of solutions are available to reduce bias: bias correction 
# and bias reduction (Firth, 1993). Bias correction methods correct for bias in a separate step following
# maximum likelihood estimation, while bias reduction methods modify the maximum likelihood estimation 
# procedure such that the resulting estimator is less biased. See Appendix S3 accompanying Case study 1
# for a demonstration of bias correction and bias reduction and a bootstrap-based technique for assessing 
# the degree of bias for any given data-model combination (Kosmidis, 2014).
summary(m5.26q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail, type="BC")) ## ## Only range sizes significant, yr and enemy release not. 
summary(m5.27q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail, link="loglog", type="BC")) ## ## Range sizes, growth en_rel_t sig. 

summary(m5.28q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail, type="BR")) ## ## Only range sizes significant, yr and enemy release not. 
summary(m5.29q <- betareg(expansion_t ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail, link="loglog", type="BR")) ## ## Range sizes, growth en_rel_t sig. 

AIC(m5q, m5.23q, m5.26q, m5.27q, m5.28q, m5.29q) ## The ML models perform best - continue to use default

###### Quantifying bias - didn't finish. Probably no bias as the bias correction methods gave the same values. #####
# extract means and variance-covariance matrix from model object
means <- unlist(coef(m5q))
vc <- vcov(m5q)

# sample from distribution of parameters assuming multivariate normal (accounts for correlation of parameters)
N <- 1000 # global parameter for number of replications
rnd <- rmvnorm(N, mean = means, sigma = vc) # returns an N x p matrix of coefficients

dm <- model.matrix(m5q)

# make predictions at the scale of the linear predictor
per <- rnd[, 1:11] %*% t(dm)
precision_lp <- rnd[, 12] %*% t(dm)

# transform to mu and precision by using the appropriate inverse link functions
library(boot)
mu <- inv.logit(per)
prec <- exp(precision_lp)

## simulate from each new set of predictors and calculate ratios
# initialize output matrix
output <- array(NA, dim = c(N, 10))

# clone the original dataset
newdata <- dat

for (i in 1:N) {
  # prediction
  newdata$ALGAE.scaled <- rbeta2(16, mu[i, ], prec[i, ])
  # fit model
  bet1 <- try(betareg(formula = ALGAE.scaled ~ treat | treat, data = newdata))
  # store output
  output[i, 1:8] <- coef(bet1)
  output[i, 9] <- logLik(bet1)
  output[i, 10] <- AIC(bet1)
}

# remove failed optimisations, with AIC >0
output <- output[output[, 10] < 0, ]

# compute bootstrap mean for each parameter
boots_mean <- apply(output, 2, mean)

##### Zero-inflation? https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13234 #####
sort(dat$expansion) # 9 0-points and 1 1-point, so try zero-augmented beta regression?
## The added value of a zero-inflated/augmented models will be larger when both zeros and relatively high proportions are observed within replicates of the same treatment - evidence that two data-generating processes are operating.
## Since the application of regular beta regression to data with zeros (and/or ones) requires transformation of the data, formal model selection criteria such as AIC or Bayesian Information Criterion (BIC) cannot be applied to compare the fit of a beta regression model fitted to a transformed response to zero-and/or-one inflated beta regression fit to an untransformed response. Therefore, model selection needs to be based on other criteria such as visual inspection of residuals and comparison of model predictions with observations
## Because response variables differ, cannot directly compare classical and zero-inflated models using standard model selection criteria such as AIC or likelihood ratio tests.

## Can't include values of 1
m5.23q.zoi <- brm(brmsformula(expansion ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, family=zero_inflated_beta()), data=dat[dat$expansion<1,], iter=5000)

## check trace plots 
traceplot(m5.23q.zoi$fit) ## all seem fine. , pars=c("b_Intercept", "b_phi_Intercept", "b_TREATrem", "b_TREATt0.33"))

## check posterior mean predictions
summary(m5.23q.zoi)

## Compare to results from classical model - different numerically and in order
coefs.zoi <- summary(m5.23q.zoi)$fixed[,"Estimate"]
coefs.clas <- summary(m5.23q)$coefficients$mean[,"Estimate"]

## Quadratic model with logit transformation
summary(m5.logitq <- lm(logit(expansion_t) ~ poly(rsEURm2_t, 2) + poly(rsUSAm2_t, 2) + poly(yr1, 2) + poly(en_rel_t, 2) + growth, data=dat, na.action=na.fail)) ## Only range sizes and growth significant, yr and enemy release not. 

plot(coefs.zoi ~ coefs.clas,
     ylab = "Parameter zero-inflated",
     xlab = "Parameter classical model"
)

par(mar=c(9, 4, 4, 2) + 0.1)
plot(coefs.zoi ~ I(1:11),
     pch = 16, cex = 1.5,
     ylim=c(-12,6),
     axes=F,
     ylab = "Parameter estimate",
     xlab = "",
     main="Coefficients from classical and zero-inflated beta-regressions\nand logit transformed linear models"
)
abline(a=0, b=0, lty=3)
axis(side=1, at=c(1:11), labels=rownames(summary(m5.2q)$coefficients), las=2)
axis(side=2, at=c(-10,-5,0,5), labels=c(-10,-5,0,5))
box()

points(coefs.clas ~ I((1:11) + 0.1), pch = 16, col = "red", cex = 1.5)
points(summary(m5q)$coefficients$mean[,"Estimate"] ~ I((1:11) + 0.1), pch = 16, col = "green", cex = 1.5) ## Quadratic model with betaregression
points(summary(m5.logitq)$coefficients[,"Estimate"] ~ I((1:11) + 0.1), pch = 16, col = "blue", cex = 1.5) ## Quadratic model with lm and logit transformed response

legend("bottomright", legend = c("z-infl", "classical loglog betareg", "classical logit betareg", "logit transform lm"), pch = 16, pt.cex = 1.5, col = c("black","red","green","blue"), bty = "n")

### Parameter estimates very similar with classical and beta-inflated betaregression, though the linear effect of rsEUR estimated by logit transform lm is a bit lower than other estimates
