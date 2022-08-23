################################################################################### 
##### pSEM submodels #####
##### Construct the individual models that could make up the global SEM #####
##### See code 2b for the submodel including niche expansion #####
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

setwd("E:\\NON_PROJECT\\TEACHING\\CORNWALL\\UG_PROJECTS\\2014-15\\JACK_HARNESS\\PATH_ANALYSIS\\2022")

dat <- read.csv("regan_dataset_analyse2_transformed.csv")

### Use 'woody' or not rather than all three growth forms. 
## There are only 9 shrubs, so little point coding them distinctly. 
## Having a binary variable makes it easier to examine interactions in SEM. 
## See powerpoint for previous analyses with shrub and tree.
## Results not meaningfully changed by replacing three growth forms with two.
dat$woody <- 0
dat$woody[dat$growth %in% c("Shrub","Tree")] <- 1

##### Submodel 1 #####
plot(yr1 ~ rsEURm2_t, data=dat)

### Linear model
summary(m1 <- lm(yr1 ~ rsEURm2_t, data=dat))
autoplot(m1, label.size = 1) ## Seems ok

R2 <- round(summary(m1)$r.squared, 4)
int <- round(coefficients(m1)["(Intercept)"], 4)
slope <- round(coefficients(m1)["rsEURm2_t"], 7)
p <- round(summary(m1)$coefficients[8], 4)

p1 <- ggplot(dat, aes(x=rsEURm2_t, y=yr1)) + geom_point(shape=1) + geom_smooth(method=lm, se=T)
p1 <- p1 + annotate("text", x=600000, y=1740, label=paste0("R^2=", R2)) +
  annotate("text", x=600000, y=1720, label=paste0("intercept=", int)) +
  annotate("text", x=600000, y=1700, label=paste0("slope=",slope)) +
  annotate("text", x=600000, y=1680, label=paste0("p=",p)) ## p-value of range size
plot(p1)

## Quadratic model
summary(m1q <- lm(yr1 ~ rsEURm2_t + I(rsEURm2_t^2), data=dat)) ## Doesn't fit
autoplot(m1q, label.size = 1) ## Seems fine

## Stepwise model
summary(m1s <- step(m1q)) ## Selects model with single linear term + intercept

##### Submodel 2 #####
plot(rsUSAm2_t ~ rsEURm2, data=dat)

### Linear model
summary(m2 <- lm(rsUSAm2_t ~ rsEURm2_t, data=dat))
autoplot(m2, label.size = 1) ## Seems fine

## Quadratic model
summary(m2q <- lm(rsUSAm2_t ~ rsEURm2_t + I(rsEURm2_t^2), data=dat)) ## Linear term sig. quadratic term marginally NS
autoplot(m2q, label.size = 1) ## Seems fine

## Stepwise model
summary(m2s <- step(m2q)) ## Selects model with both linear and quadratic term + intercept

p2 <- ggplot(dat, aes(x=rsEURm2_t)) +
               geom_point(aes(y=rsUSAm2_t), shape=1) + 
               stat_smooth(aes(y=rsUSAm2_t), method="lm", formula=y ~ x + I(x^2), size = 1)

R2 <- round(summary(m2s)$r.squared, 4)
int <- round(coefficients(m2s)["(Intercept)"], 0)
l <- round(coefficients(m2s)["rsEURm2_t"], 4)
q <- round(coefficients(m2s)["I(rsEURm2_t^2)"], 7)

p2 <- p2 + annotate("text", x=200000, y=2000000, label=paste0("R^2=", R2)) +
  annotate("text", x=300000, y=1800000, label=paste0("intercept=", int)) +
  annotate("text", x=300000, y=1600000, label=paste0("linear term=", l)) +
  annotate("text", x=300000, y=1400000, label=paste0("quad term=", q))

plot(p2)

##### Submodel 3 - area corrected #####
par(mfrow=c(2,2))
par(mar=c(4, 4, 1, 1)) # c(bottom, left, top, right)
plot(density_europe_t ~ rsEURm2, data=dat)
plot(density_europe_t ~ citations_t, data=dat)
plot(density_europe_t ~ growth, data=dat)

### Linear model
summary(m3 <- lm(density_europe_t ~ rsEURm2_t + citations_t + woody, data=dat))
autoplot(m3, label.size = 1) ## Seems fine

## Quadratic model
## Using q and I(q^2) is problematic because they will be correlated and correlated variables can cause problems. 
## The use of poly() avoids this by producing orthogonal polynomials
## However, in order to test for interactions between linear effects and to dredge, have to create new variables for first and second order polynomial
dat$rsEURm2_t1 <- poly(dat$rsEURm2_t, 2)[,1]
dat$rsEURm2_t2 <- poly(dat$rsEURm2_t, 2)[,2]
dat$citations_t1 <- poly(dat$citations_t, 2)[,1]
dat$citations_t2 <- poly(dat$citations_t, 2)[,2]

summary(m3q <- lm(density_europe_t ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody, data=dat, na.action=na.fail)) ## all relationships significant. 
autoplot(m3q, label.size = 1) ## Seems fine
summ(m3q) 

## Stepwise model
summary(m3s <- step(m3q)) ## Selects model with both linear and quadratic terms + intercept

## Check for interactions
summary(m3i <- lm(density_europe_t ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody +
                    rsEURm2_t1*citations_t1 + rsEURm2_t1*woody + citations_t1*woody,
                  data=dat, na.action=na.fail)) ## na.action is for if dredge is needed.
summary(m3is <- step(m3i)) ## selects model with all linear and quadratic effects and with interaction between woody and citations_t1, EUR range size and citations, and EUR range size and woody (when all three growth forms were used, no interactions were selected)
autoplot(m3is, label.size = 1) ## Seems fine
summ(m3is) 

## Plot best fitting model
## all non-plotted variables are centered by default, i.e. relationship plotted given mean values of other explanatory variables
effect_plot(m3is, pred=rsEURm2_t1, interval=T, plot.points=T, data=dat, main.title="Observed data")
effect_plot(m3is, pred=rsEURm2_t1, interval=T, plot.points=T, partial.residuals=T, data=dat, main.title="Partial residuals")

effect_plot(m3is, pred=citations_t1, interval=T, plot.points=T, data=dat, main.title="Observed data") 
effect_plot(m3is, pred=citations_t1, interval=T, plot.points=T, partial.residuals=T, data=dat, main.title="Partial residuals")

effect_plot(m3is, pred=woody, interval=T, plot.points=T, data=dat, main.title="Observed data")
effect_plot(m3is, pred=woody, interval=T, plot.points=T, partial.residuals=T, data=dat, main.title="Partial residuals")

## Variance Inflation Factor
vif(m3is) ## all less than 5: OK

## Check final model no different when use first order polynomial as response
dat$density_europe_t1 <- poly(dat$density_europe_t, 2)[,1]
dat$density_europe_t2 <- poly(dat$density_europe_t, 2)[,2]

summary(m3i.v2 <- lm(density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody +
                    rsEURm2_t1*citations_t1 + rsEURm2_t1*woody + citations_t1*woody, data=dat))

plot(summary(m3i)$coefficients[,"Estimate"], summary(m3i.v2)$coefficients[,"Estimate"]) ## Effect sizes remain the same (apart from intercept, which is different at the numerical value of the response variable differs)

##### Submodel 3 - NOT area corrected #####
### This is from when growth form was coded as three variables - now obsolete. 
### Did not repeat with growth coded as 'woody' as decided not to use area correction
par(mfrow=c(2,2))
par(mar=c(4, 4, 1, 1)) # c(bottom, left, top, right)
plot(fungi_eur ~ rsEURm2, data=dat)
plot(fungi_eur ~ citations_t, data=dat)
plot(fungi_eur ~ growth, data=dat)

### Linear model
summary(m33 <- lm(fungi_eur_t ~ rsEURm2_t + citations_t + growth, data=dat)) ## All linear terms significant
autoplot(m33, label.size = 1) ## Seems fine

## Quadratic model
summary(m33q <- lm(fungi_eur_t ~ poly(rsEURm2_t, 2) + poly(citations_t, 2) + growth, data=dat, na.action=na.fail)) ## all relationships significant except quadratic term of rsEURm2_t. 
autoplot(m33q, label.size = 1) ## Seems fine
summ(m33q) 

## Stepwise model
summary(m33s <- step(m33q)) ## Selects model with all linear and quadratic terms + intercept

## Check for interactions
summary(m33i <- lm(fungi_eur_t ~ poly(rsEURm2_t, 2) * poly(citations_t, 2) * growth,
                  data=dat, na.action=na.fail)) ## Native range size and shrub interactions significant. Citations and shrub interactions near significant. 3-way interaction between ntv range size (quad), citations, and shrub significant. Argh!!
summary(m33is <- step(m33i)) ## selects model with all linear and quadratic effects and without interactions. Phew.

## Plot best fitting model
## all non-plotted variables are centered by default, i.e. relationship plotted given mean values of other explanatory variables
effect_plot(m33q, pred=rsEURm2_t, interval=T, plot.points=T, data=dat, main.title="Observed data")
effect_plot(m33q, pred=rsEURm2_t, interval=T, plot.points=T, partial.residuals=T, data=dat, main.title="Partial residuals")

effect_plot(m33q, pred=citations_t, interval=T, plot.points=T, data=dat, main.title="Observed data") 
effect_plot(m33q, pred=citations_t, interval=T, plot.points=T, partial.residuals=T, data=dat, main.title="Partial residuals")

effect_plot(m33q, pred=growth, interval=T, plot.points=T, data=dat, main.title="Observed data")
effect_plot(m33q, pred=growth, interval=T, plot.points=T, partial.residuals=T, data=dat, main.title="Partial residuals")

## Variance Inflation Factor
vif(m33q) ## all less than 5: OK

##### Submodel 4 #####
par(mfrow=c(2,3))
par(mar=c(4, 4, 1, 1)) # c(bottom, left, top, right)
plot(en_rel_t ~ rsUSAm2, data=dat)
plot(en_rel_t ~ yr1, data=dat)
plot(en_rel_t ~ density_europe_t, data=dat)
plot(en_rel_t ~ citations_t, data=dat)
plot(en_rel_t ~ woody, data=dat)

### Linear model
summary(m4 <- lm(en_rel_t ~ rsUSAm2_t + yr1 + density_europe_t + citations_t + woody, data=dat))
autoplot(m4, label.size = 1) ## Seems to have broader variance at mid-levels of enemy release than at high or low ends - caused by a few rogue values in the centre of the distribution of the response variable. Leverage graph is flatlined, so none of these values seem to have undue influence.

### Quadratic model
## Using q and I(q^2) is problematic because they will be correlated and correlated variables can cause problems. 
## The use of poly() avoids this by producing orthogonal polynomials
## However, in order to test for interactions between linear effects and to dredge, have to create new variables for first and second order polynomial
dat$rsUSAm2_t1 <- poly(dat$rsUSAm2_t, 2)[,1]
dat$rsUSAm2_t2 <- poly(dat$rsUSAm2_t, 2)[,2]
dat$yr1_1 <- poly(dat$yr1, 2)[,1]
dat$yr1_2 <- poly(dat$yr1, 2)[,2]

summary(m4q <- lm(en_rel_t ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 +
                    density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody,
                  data=dat, na.action=na.fail)) ## yr1 not significant. Woody significant (GrowthTree was only just significant previously)
autoplot(m4q, label.size = 1) ## Seems ok
summ(m4q) 

### Stepwise model
summary(m4s <- step(m4q)) ## Selects model with all linear and quadratic terms + intercept
summ(m4s)
dredge(m4q) ## Best model that respects marginality contains all linear and quadratic terms 

### Check final model for whether quadratic terms and contentious main effect terms are required:
## Note this investigation was done with three growth forms. Did not repeat as used the model identified by step.
## yr1 and growth could both be discarded and dAIC remaining < 2
AIC(m4q) # -1856.189
summary(m4q.1 <- lm(en_rel_t ~ poly(rsUSAm2_t, 2) + yr1 +
              poly(density_europe_t, 2) + poly(citations_t, 2) + growth,
            data=dat, na.action=na.fail))
AIC(m4q.1) # -1855.348. Can be considered equal to model with quadratic effect of yr1

summary(m4q.2 <- lm(en_rel_t ~ poly(rsUSAm2_t, 2) + 
                      poly(density_europe_t, 2) + poly(citations_t, 2) + growth,
                    data=dat, na.action=na.fail)) 
AIC(m4q.2) # -1855.285 ## Can be considered equal to model containing yr1 with or without quadratic effect.

summary(m4q.3 <- lm(en_rel_t ~ poly(rsUSAm2_t, 2) + poly(yr1, 2) +
                      poly(density_europe_t, 2) + citations_t + growth,
                    data=dat, na.action=na.fail)) 
AIC(m4q.3) # -1850.657. Delta AIC > 2 so reject.

summary(m4q.4 <- lm(en_rel_t ~ poly(rsUSAm2_t, 2) + poly(yr1, 2) +
                      poly(density_europe_t, 2) + growth,
                    data=dat, na.action=na.fail)) 
AIC(m4q.4) # -1840.208. Delta AIC > 2 so reject.

summary(m4q.5 <- lm(en_rel_t ~ poly(rsUSAm2_t, 2) + poly(yr1, 2) +
                      poly(density_europe_t, 2) + poly(citations_t, 2),
                    data=dat, na.action=na.fail))
AIC(m4q.5) # -1854.857 Can be considered equal to model with growth.

summary(m4q.6 <- lm(en_rel_t ~ poly(rsUSAm2_t, 2) + 
                      poly(density_europe_t, 2) + poly(citations_t, 2),
                    data=dat, na.action=na.fail)) 
AIC(m4q.6) # -1854.577. Can be considered equal to model with growth and yr1.

effect_plot(m4q.6, pred=rsUSAm2_t, interval=T, plot.points=T, partial.residuals=F, data=dat, main.title="Observed data")
effect_plot(m4q.6, pred=rsUSAm2_t, interval=T, plot.points=T, partial.residuals=T, data=dat, main.title="Partial residuals")

effect_plot(m4q.6, pred=density_europe_t, interval=T, plot.points=T, partial.residuals=F, data=dat, main.title="Observed data")
effect_plot(m4q.6, pred=density_europe_t, interval=T, plot.points=T, partial.residuals=T, data=dat, main.title="Partial residuals")

effect_plot(m4q.6, pred=citations_t, interval=T, plot.points=T, partial.residuals=F, data=dat, main.title="Observed data")
effect_plot(m4q.6, pred=citations_t, interval=T, plot.points=T, partial.residuals=T, data=dat, main.title="Partial residuals")

## Check for interactions. 
summary(m4q <- lm(en_rel_t ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 +
                    density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody,
                  data=dat, na.action=na.fail)) ## yr1 not significant. Woody significant (GrowthTree was only just significant previously)


summary(m4i <- lm(en_rel_t ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 +
                    density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +
                    rsUSAm2_t1*yr1_1 + rsUSAm2_t1*density_europe_t1 + rsUSAm2_t1*citations_t1 + rsUSAm2_t1*woody +
                    density_europe_t1*yr1_1 + density_europe_t1*citations_t1 + density_europe_t1*woody +
                    citations_t1*yr1_1 + citations_t1*woody +
                    yr1_1*woody,                                                                   , 
                  data=dat, na.action=na.fail)) ## na.action is for if dredge is needed. density_europe (linear and curved), USA range-size (curved) only variables significant at 0.05
summary(m4is <- step(m4i)) ## selects model without quadratic effect of yr1, and with interaction between rsUSA and density_europe, and between yr1 and woody (when growth coded as three forms, yr1*citations was selected and no interaction between growth form and any other variable)
m4i.d <- dredge(m4i) ## best model that respects marginality excludes quadratic effect of yr1, and with interaction between rsUSA and density_europe, and between yr1 and woody (when growth coded as three forms, yr1*citations was selected and no interaction between growth form and any other variable)
m4i.db <- get.models(m4i.d, subset=1)$'270592'
summ(m4i.db)
autoplot(m4i.db)
vif(m4q) ## all fine. VIF won't work when applied to model with interactions but that's fine.