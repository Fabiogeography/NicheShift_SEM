################################################################################### 
##### Make pSEM from submodels #####
##### Written by Regan Early #####
##### Written on: 26th August 2021 #####
##### Modified by: Regan Early #####
##### Modified on: 22nd August 2022 #####
##############################################################################
##### https://cran.r-project.org/web/packages/piecewiseSEM/vignettes/piecewiseSEM.html#correlated-errors #####

.libPaths("C:/SOFTWARE/R/R-4.1.2/library")
library(piecewiseSEM)
library(jtools) ## effect_plot() - plots partial residuals in multiple regressions - and summ() - super neat model summary. Needs latest version of ggplot2
library(ggplot2)
library(ggfortify) ## nice diagnostic plots
library(betareg)
library(car) ## logit transformation

setwd("E:\\NON_PROJECT\\TEACHING\\CORNWALL\\UG_PROJECTS\\2014-15\\JACK_HARNESS\\PATH_ANALYSIS\\2022")

dat <- read.csv("regan_dataset_analyse2_transformed.csv")

##### Prepare a new data frame #####
## Quadratic effects need to be coded as a variable from the data frame, otherwise pSEM doesn't recognise them.
## Assign the values of poly() to the data frame
## using poly() doesn't change the model from the one specified in the submodel code using (I(...^2)), although linear effects become more significant.
dat$rsEURm2_t1 <- poly(dat$rsEURm2_t, 2)[,1]
dat$rsEURm2_t2 <- poly(dat$rsEURm2_t, 2)[,2]
dat$rsUSAm2_t1 <- poly(dat$rsUSAm2_t, 2)[,1]
dat$rsUSAm2_t2 <- poly(dat$rsUSAm2_t, 2)[,2]
dat$yr1_1 <- poly(dat$yr1, 2)[,1]
dat$yr1_2 <- poly(dat$yr1, 2)[,2]
dat$density_europe_t1 <- poly(dat$density_europe_t, 2)[,1]
dat$density_europe_t2 <- poly(dat$density_europe_t, 2)[,2]
dat$citations_t1 <- poly(dat$citations_t, 2)[,1]
dat$citations_t2 <- poly(dat$citations_t, 2)[,2]
dat$en_rel_t1 <- poly(dat$en_rel_t, 2)[,1]
dat$en_rel_t2 <- poly(dat$en_rel_t, 2)[,2]
dat$expansion_tl <- logit(dat$expansion_t)
dat$woody <- 0
dat$woody[dat$growth %in% c("Shrub","Tree")] <- 1

##### Calculate submodels. #####
## Note that responses and explanatory variables needs to be the same throughout, i.e. the linear terms of the polynomial calculation

### Submodel 1
m1 <- lm(yr1_1 ~ rsEURm2_t1, data=dat) ## Needs to be the same explanatory variable as used (and calculated) below

### Submodel 2. 
m2q <- lm(rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2, dat) 

### Submodel 3. 
m3i.v2 <- lm(density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1*citations_t1 + rsEURm2_t1*woody + citations_t1*woody, data=dat)

### Submodel 4. 
m4i.db <- lm(en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +
               rsUSAm2_t1*density_europe_t1 + yr1_1*woody, 
             data=dat)

### Submodel 5. logit-transformed lm. Best of best model subset
m5i.poly <- lm(expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + woody + woody*rsEURm2_t1, dat)

##### Construct and investigate the SEM #####
p1 <- psem(m1, m2q, m3i.v2, m4i.db, m5i.poly) ## 
summary(p1, .progressBar = F) ## Fisher's C-test: if result is **non-significant** at 5% you aren't missing any relationships.

## Update submodel 2 with possible missing relationships
summary(m2qu <- lm(rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2+ yr1_1 + woody, data=dat))

## Update submodel 5 with possible missing relationships
## Add curved and linear effects of density as wouldn't normally just add cured effect
summary(m5i.polyu <- lm(expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + woody*rsEURm2_t1, data=dat))

## New psem
p2 <- psem(m1, m2qu, m3i.v2, m4i.db, m5i.polyu, citations_t2%~~%expansion_tl, citations_t1%~~%yr1_1) 
summary(p2, .progressBar = F) ## Fisher's C-test: only density_europe_t1%~~%density_europe_t2 covariance remains. Confusing.

### See if coding relationship between naturalised range size and enemy release as covariance affects the model
## Update m4 to allow relationship between naturalised range size and enemy release to be coded as covariance
m4i.db.cov <- lm(en_rel_t1 ~ yr1_1 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody, 
             data=dat)

p3 <- psem(m1, m2qu, m3i.v2, m4i.db.cov, m5i.polyu,
               citations_t2%~~%expansion_tl, citations_t1%~~%yr1_1, en_rel_t1%~~%rsUSAm2_t1, en_rel_t1%~~%rsUSAm2_t2)
summary(p3, .progressBar = F)

## Extract the relationships
## The most typical implementation of standardization is placing the coefficients in units of standard deviations of the mean.
coefs(p1, standardize = "scale") 
coefs(p2, standardize = "scale") ## The final piecewiseSEM
coefs(p3, standardize = "scale") ## Although coding two relationships as covariance changes that model slightly (notably effect size of woody), the expansion model doesn't change.

write.csv(coefs(p2, standardize = "scale"), file="piecewiseParams_nonPhylo.csv", row.names=F)
