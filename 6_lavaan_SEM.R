################################################################################### 
##### Test SEM from pSEM with ML in Lavaan #####
##### Modified to work with lavaan version 2022 as the update seems to have changed the best model #####
##### Written by Regan Early #####
##### Written on: 6th September 2021 #####
##### Modified by: #####
##### Modified on: 23rd August 2022 #####
##############################################################################

.libPaths("D:/SOFTWARE/R-4.1.2/library")
library(lavaan)
library(lavaanPlot)
library(car) ## logit transformation
library(semTools) ## toolbox for sem including normality and kurtosis tests. 
library(MVN)
library(stringr)

setwd("E:\\NON_PROJECT\\TEACHING\\CORNWALL\\UG_PROJECTS\\2014-15\\JACK_HARNESS\\PATH_ANALYSIS\\2022")

dat <- read.csv("regan_dataset_analyse2_transformed.csv")

##### Code as coded in piecewiseSEM #####
## Quadratic effects need to be coded as a variable from the data frame, otherwise pSEM doesn't recognise them.
## sem fails if (^2) or poly(, 2) is included in teh formulae
# Read this: https://groups.google.com/g/lavaan/c/ZKhKQq2fq1E?pli=1. Basically suggests doing what I already did, and explains how to calculate the coefficient for the linear + quadradic relationship
# Other group chats suggest I did it the right way
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

# ## Interactions
# dat$rsEURm2_t1Xcitations_t1 <- dat$rsEURm2_t1*dat$citations_t1
# dat$rsEURm2_t1Xwoody <- dat$rsEURm2_t1*dat$woody
# dat$citations_t1Xwoody <- dat$citations_t1*dat$woody
# dat$rsUSAm2_t1Xdensity_europe_t1 <- dat$rsUSAm2_t1*dat$density_europe_t1
# dat$yr1_1Xwoody <- dat$yr1_1*dat$woody
# dat$woodyXrsEURm2_t1 <- dat$woody*dat$rsEURm2_t1

##### Test for multivariate normality ##### 
mardiaSkew(dat[,56:68]) ## not normal
mardiaKurtosis(dat[,56:68]) ## not normal
## Use alternative estimators

##### Run SEMs #####
## These regression formulas are similar to the way ordinary linear regression formulas are used
## in R, but they may include latent variables. Interaction terms are currently not supported.
## Include correlated errors (~~)
## Note interactions must be coded with : (not)
submodels.i <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + woody:rsEURm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
submodels.ii <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + rsEURm2_t1:woody + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody 
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
 ## Does not converge if submodels 3 and 5 both contain woody:rsEURm2_t1

## Note this structure dictates the covariance of citations_t1+t2, which requires them to be random effects.
## Because they are exogenous (nothing affects them) lavaan by default wants their covariance to be fixed to the sample vales (i.e. be fixed effects)
## Setting fixed.x=F within sem sets all exogenous variables to be random and prevents an error message, while seemingly not changing the model output. 
## See ?lavOptions ", fixed.x=F"
summary(sem1.i <- sem(submodels.i, dat), standardized=T, rsq=T, fit.measures=T) 
summary(sem1.ii <- sem(submodels.ii, dat), standardized=T, rsq=T, fit.measures=T) 
anova(sem1.i, sem1.ii) ## sem1.i fits significantly better
sem1 <- sem1.i ## remove interaction from submodel 3
submodels <- submodels.i

## The first column (Estimate) contains the (estimated or fixed) parameter value for each model parameter;
## the second column (Std.err) contains the standard error for each estimated parameter;
## the third column (Z-value) contains the Wald statistic (which is simply obtained by dividing the parameter value by its standard error)
## the fourth column (P(>jzj)) contains the p-value for testing the null hypothesis that the parameter equals zero in the population.
## In the fifth column (labeled Std.lv), only the latent variables are standardized.
## If "std.lv", the standardized estimates are on the variances of the (continuous) latent variables only
## In the sixth column (labeled Std.all), both latent and observed variables are standardized.
## If "std.all", the standardized estimates are based on both the ***variances of both (continuous) observed*** and latent variables, including the variances of exogenous covariates.
## The latter is often called the `completely standardized solution'.

## Note that in the Variances: section, there is a dot before the observed variables names.
## This is because they are dependent (or endogenous) variables (predicted by the latent variables),
## and therefore, the value for the variance that is printed in the output is an estimate of the residual variance:
## the left-over variance that is not explained by the predictor(s). 
## By contrast, there is no dot before the latent variable names, because they are exogenous variables
## in this model (there are no single-headed arrows pointing to them).
## The values for the variances here are the estimated total variances of the latent variables.

standardizedSolution(sem1) ## default is Std.all. Gives slightly different results to summary, particularly for growth, where the values are an order of magnitude different (and different from the final solution in piecewiseSEm as well). Why?
coef(sem1)

lavaanPlot(model=sem1, coefs=T)
lavaanPlot(model=sem1, coefs=T, stand=T) ## standardized coefficients

##### Try to improve fit of sem #####
### Investigate magnitude of deviation between the relationships in the data and those implied by the model.
round(inspect(sem1, "sample")$cov[1:9, 1:9], 4)
round(fitted(sem1)$cov[1:9, 1:9], 4)

round(inspect(sem1, "sample")$cov[10:18, 10:18], 4)
round(fitted(sem1)$cov[10:18, 10:18], 4)

### Include woody ~~ citation covariances
submodels2 <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + woody:rsEURm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
citations_t1~~woody
citations_t2~~woody
'
summary(sem2 <- sem(submodels2, dat), standardized=T, rsq=T, fit.measures=T) 

round(inspect(sem2, "sample")$cov[10:18, 10:18], 4)
round(fitted(sem2)$cov[10:18, 10:18], 4)## fitted and observed now match for citations_t1 and citation_t2

### Include causal relationship of yr1 on expansion as well
submodels3 <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + yr1_1 + woody:rsEURm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
citations_t1~~woody
citations_t2~~woody
'
summary(sem3 <- sem(submodels3, dat), standardized=T, rsq=T, fit.measures=T) 

round(inspect(sem3, "sample")$cov[1:9, 1:9], 4)
round(fitted(sem3)$cov[1:9, 1:9], 4) ## fitted and observed now closer for expansion~yr1, but effect of yr1_1 is not significant

### Compare three alternative models
anova(sem1, sem2, sem3) 

## LRtest reveals sig between the models 1 and 2, and AIC supports sem1 as the best model so the simpler model (sem1) is better. 
AIC(sem1); AIC(sem2); AIC(sem3)
## However, neither model appears to be correctly specified (Chi-square value)
## Causal relationship of yr1_ on expansion_tl, and covariance between citations and growth NOT supported

##### Test for multivariate normality in endogenous variables - are alternative estimators valid? #####
endovars <- dat[,c("yr1_1", "rsUSAm2_t1", "density_europe_t1", "en_rel_t1", "expansion_tl")]
mvn(endovars, mvnTest="mardia") ## data skewed and has kurtosis
## For individual variables, skewness values approaching 2 and kurtosis values over 7 can mean more significant problems with non-normality, presumably to the extent where regressions inappropriate?
## Non of my variables do that, so hopefully piecewiseSEM OK. 

mvn(endovars, multivariatePlot="qq")
mvn(endovars, multivariateOutlierMethod="quan")
## There is evidence of multivariate non-normality

### Test alternative estimators, as ML requires multivariate normality in endogenous variables ###
summary(sem1.1 <- sem(submodels, dat, estimator="mlm"), standardized=T, rsq=T, fit.measures=T) ## Chi-sq & other three tests now all within preferred values
summary(sem2.1 <- sem(submodels2, dat, estimator="mlm"), standardized=T, rsq=T, fit.measures=T) ## Chi-sq & other three tests now all within preferred values
summary(sem3.1 <- sem(submodels3, dat, estimator="mlm"), standardized=T, rsq=T, fit.measures=T) ## Chi-sq not sig.
# write.csv(parameterEstimates(sem1.1), file="semparams1_1_mlm.csv", row.names=F)
# Better to use 2nd version of robust CFI etc... min 14 in https://www.youtube.com/watch?v=HvYW_GeHpD8. Savelei 2018.

summary(sem1.2 <- sem(submodels, dat, test="bollen.stine", se="boot", bootstrap=1000), standardized=T, rsq=T, fit.measures=T) ## Chi-sq not sig but RMSEA and CFI outside of preferred range.
summary(sem2.2 <- sem(submodels2, dat, test="bollen.stine", se="boot", bootstrap=1000), standardized=T, rsq=T, fit.measures=T) ## Chi-sq & other three tests now all within preferred values. 
summary(sem3.2 <- sem(submodels3, dat, test="bollen.stine", se="boot", bootstrap=1000), standardized=T, rsq=T, fit.measures=T) ## Chi-sq & other three tests now all within preferred values. 
# write.csv(parameterEstimates(sem1.2), file="semparams1_2_bs.csv", row.names=F) ## converted this and above to semparams1_mlm_vs_bollerstine.xls

##### Test removal of NS paths, with sem1 as a basis #####

### Remove submodel 1
submodels1.updates1 <- '
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1*citations_t1 + citations_t1*woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1*density_europe_t1 + yr1_1*woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + woody*rsEURm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
summary(sem1.updates1 <- sem(submodels1.updates1, dat, estimator="mlm"), standardized=T, rsq=T)

### Remove NS paths in submodel 2
## rsEURm2_t2
submodels1.updates2a <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + woody:rsEURm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
summary(sem1.updates2a <- sem(submodels1.updates2a, dat, estimator="mlm"), standardized=T, rsq=T) 

## yr1_1
submodels1.updates2b <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + woody:rsEURm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
summary(sem1.updates2b <- sem(submodels1.updates2b, dat, estimator="mlm"), standardized=T, rsq=T) 

## rsEURm2_t2 & yr1_1
submodels1.updates2c <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + woody:rsEURm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
summary(sem1.updates2b <- sem(submodels1.updates2c, dat, estimator="mlm"), standardized=T, rsq=T) 

### Remove NS paths in submodel 3
## citations_t2
submodels1.updates3a <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + woody:rsEURm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
summary(sem1.updates3a <- sem(submodels1.updates3a, dat, estimator="mlm"), standardized=T, rsq=T) 

## citations_t2
submodels1.updates3b <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + woody + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + woody:rsEURm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
summary(sem1.updates3b <- sem(submodels1.updates3b, dat, estimator="mlm"), standardized=T, rsq=T) 

## Remove NS paths in submodel 4
## yr1_2
submodels1.updates4a <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + woody:rsEURm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
summary(sem1.updates4a <- sem(submodels1.updates4a, dat, estimator="mlm"), standardized=T, rsq=T) 

##yr1_1 and yr1_2 and interactions
submodels1.updates4b <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + woody:rsEURm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
summary(sem1.updates4b <- sem(submodels1.updates4b, dat, estimator="mlm"), standardized=T, rsq=T) 

# ## citations_t2
# submodels1.updates4c <- '
# yr1_1 ~ rsEURm2_t1
# rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
# density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
# en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
# expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + woody:rsEURm2_t1
# citations_t2~~expansion_tl
# citations_t1~~yr1_1
# '
# summary(sem1.updates4c <- sem(submodels1.updates4c, dat, estimator="mlm"), standardized=T, rsq=T) 
# 
# ## woody and interactions
# submodels1.updates4d <- '
# yr1_1 ~ rsEURm2_t1
# rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
# density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
# en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + rsUSAm2_t1:density_europe_t1
# expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + woody:rsEURm2_t1
# citations_t2~~expansion_tl
# citations_t1~~yr1_1
# '
# summary(sem1.updates4d <- sem(submodels1.updates4d, dat, estimator="mlm"), standardized=T, rsq=T) 

## all
submodels1.updates4e <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + density_europe_t1 + density_europe_t2 + citations_t1 + rsUSAm2_t1:density_europe_t1
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + woody:rsEURm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
summary(sem1.updates4e <- sem(submodels1.updates4e, dat, estimator="mlm"), standardized=T, rsq=T) 

### Remove NS paths in submodel 5
# ## rsEURm2_t2
# submodels1.updates5a <- '
# yr1_1 ~ rsEURm2_t1
# rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
# density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
# en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
# expansion_tl ~ rsEURm2_t1 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody + woody:rsEURm2_t1
# citations_t2~~expansion_tl
# citations_t1~~yr1_1
# '
# summary(sem1.updates5a <- sem(submodels1.updates5a, dat, estimator="mlm"), standardized=T, rsq=T) 

# ## rsEURm2_t1 & rsEURm2_t2 & interaction
# submodels1.updates5b <- '
# yr1_1 ~ rsEURm2_t1
# rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
# density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
# en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
# expansion_tl ~ rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2 + woody
# citations_t2~~expansion_tl
# citations_t1~~yr1_1
# '
# summary(sem1.updates5b <- sem(submodels1.updates5b, dat, estimator="mlm"), standardized=T, rsq=T) 

## en_rel_t1
submodels1.updates5c <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + density_europe_t1 + density_europe_t2 + woody + woody:rsEURm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
summary(sem1.updates5c <- sem(submodels1.updates5c, dat, estimator="mlm"), standardized=T, rsq=T) 

# ## density_europe_t2
# submodels1.updates5d <- '
# yr1_1 ~ rsEURm2_t1
# rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
# density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
# en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
# expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + woody + woody:rsEURm2_t1
# citations_t2~~expansion_tl
# citations_t1~~yr1_1
# '
# summary(sem1.updates5d <- sem(submodels1.updates5d, dat, estimator="mlm"), standardized=T, rsq=T) 

## density_europe_t2 & density_europe_t1 
submodels1.updates5e <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + woody + woody:rsEURm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
summary(sem1.updates5e <- sem(submodels1.updates5e, dat, estimator="mlm"), standardized=T, rsq=T) 

## woody & interactions
submodels1.updates5f <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + density_europe_t1 + density_europe_t2
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
summary(sem1.updates5f <- sem(submodels1.updates5f, dat, estimator="mlm"), standardized=T, rsq=T) 

## all
submodels1.updates5g <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
summary(sem1.updates5g <- sem(submodels1.updates5g, dat, estimator="mlm"), standardized=T, rsq=T) 

### Compare models with paths removed
sem1.mlm <- sem(submodels, dat, estimator="mlm") 
## lrtest
# anova(sem1.mlm,sem1.updates1,sem1.updates2a,sem1.updates2b,sem1.updates3a,sem1.updates4a,sem1.updates4b,sem1.updates4e,sem1.updates5a,sem1.updates5b,sem1.updates5c,sem1.updates5d,sem1.updates5e,sem1.updates5f,sem1.updates5g)
anova(sem1.mlm,sem1.updates1,sem1.updates2a,sem1.updates2b,sem1.updates3a,sem1.updates4a,sem1.updates4b,sem1.updates4e,sem1.updates5c,sem1.updates5e,sem1.updates5f,sem1.updates5g)

## AIC
x <- as.data.frame(AIC(sem1.mlm,sem1.updates1,sem1.updates2a,sem1.updates2b,sem1.updates3a,sem1.updates4a,sem1.updates4b,sem1.updates4e,sem1.updates5c,sem1.updates5e,sem1.updates5f,sem1.updates5g))
x[order(x$AIC),] 
## sem1.updates4a AIC -2938.619 with delta AIC of ~1
## removes yr1_2 from submodel 4 
## sem1.updates5d is next best, which removes just density_europe_t2 

### Test some combinations of dropped variables
## density_europe_t2 & density_europe_t2 (submodel 5e) & yr1_2 (submodel 4a)
submodels1.updates6a <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1 + yr1_1:woody
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + woody + woody:rsEURm2_t1
citations_t2~~expansion_tl
citations_t1~~yr1_1
'
summary(sem1.updates6a <- sem(submodels1.updates6a, dat, estimator="mlm"), standardized=T, rsq=T) 

# ## density_europe_t2 & density_europe_t2 (submodel 5) & yr1_1 and yr1_2 and interactions  (submodel 4)
# submodels1.updates6b <- '
# yr1_1 ~ rsEURm2_t1
# rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
# density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
# en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1
# expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + woody + woody:rsEURm2_t1
# citations_t2~~expansion_tl
# citations_t1~~yr1_1
# '
# summary(sem1.updates6b <- sem(submodels1.updates6b, dat, estimator="mlm"), standardized=T, rsq=T) 

# ## density_europe_t2 & density_europe_t2 & woody & interactions(submodel 5) & yr1_1 and yr1_2 and interactions (submodel 4)
# submodels1.updates6c <- '
# yr1_1 ~ rsEURm2_t1
# rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + woody
# density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + woody + rsEURm2_t1:citations_t1 + citations_t1:woody
# en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + woody +  rsUSAm2_t1:density_europe_t1
# expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1
# citations_t2~~expansion_tl
# citations_t1~~yr1_1
# '
# summary(sem1.updates6c <- sem(submodels1.updates6c, dat, estimator="mlm"), standardized=T, rsq=T) 

xx <- as.data.frame(AIC(sem1.updates6a,sem1.mlm,sem1.updates1,sem1.updates2a,sem1.updates2b,sem1.updates3a,sem1.updates4a,sem1.updates4b,sem1.updates4e,sem1.updates5c,sem1.updates5e,sem1.updates5f,sem1.updates5g))
xx[order(xx$AIC),] ## 4a still better

##### Compare fit with MLM and bollerstine #####
summary(sem1.updates4a.bs <- sem(submodels1.updates4a, dat, test="bollen.stine", se="boot"), standardized=T, rsq=T) 

# write.csv(parameterEstimates(sem1.updates4a), file="semparams1_updates4a_mlm.csv", row.names=F) 
# write.csv(parameterEstimates(sem1.updates4a.bs), file="semparams1_updates4a_bs.csv", row.names=F) 

## Use standardizedSolution rather than parameterEstimates to obtain SEs and test statistics for standardized estimates
write.csv(standardizedsolution(sem1.updates4a), file="semparams1_updates4a_mlmS.csv", row.names=F) 
write.csv(standardizedsolution(sem1.updates4a.bs), file="semparams1_updates4a_bsS.csv", row.names=F) 

### Model nesting: http://web.pdx.edu/~newsomj/semclass/ho_nested.pdf
## A nested model is a model that uses the same variables (and cases!) as another model but
## specifies at least one additional parameter to be estimated. The model with fewer restrictions or more free
## parameters (i.e., fewer degrees of freedom), which could be called a reduced model, is nested within the
## more restricted model, which could be called the full model

## lavaan WARNING: some models are based on a different set of observed variables.
## --> THis seems to just mean the variables in the model, which will differ between models unless you're just altering the relationships. https://rdrr.io/cran/lavaan/src/R/lav_test_LRT.R
## lavaan WARNING: some restricted models fit better than less restricted models; either these models are not nested, or the less restricted model failed to reach a global optimum.
## --> The Chisq of sem2b is less than the Chisq of sem2a. The function is supposed to order the models in order of degrees of freedom, so the earlier model should have fewer df and therfore be the mroe restricted(??)
## --> The test is then based on the order of variables. However, if two models have the same df, the original order is kept. So you don't get the error when swap sems in the function, and I don't think it's meaningful. https://rdrr.io/cran/lavaan/src/R/lav_test_LRT.R
## lavaan WARNING: some models have the same degrees of freedom
## --> sem2a and sem2b have the same df. 

##### Change some causal relationships for covariances #####
# Humped effect of natd range size on enemy release is counter-intuitive. Perhaps covariance rather than effect?
submodels1e.cov <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + growthTree + growthShrub
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + growthTree + growthShrub
en_rel_t1 ~ yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + growthTree + growthShrub + rsUSAm2_t1:density_europe_t1
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + rsUSAm2_t1 + en_rel_t1 + growthTree + growthShrub
citations_t2~~expansion_tl
citations_t1~~yr1_1
en_rel_t1~~rsUSAm2_t1
en_rel_t1~~rsUSAm2_t2
'

summary(sem1e.cov <- sem(submodels1e.cov, dat, test="bollen.stine", se="boot"), standardized=T, rsq=T, fit.measures=T) ## CFI too low but Chi-sq & other two tests all within preferred values
write.csv(parameterEstimates(sem1e.cov), file="semparams1e_cov.csv", row.names=F) ## converted this and above to semparams1_mlm_vs_bollerstine.xls

AIC(sem1e) ## -1525.667
AIC(sem1e.cov) ## -1787.048 Much lower

lavTestLRT(sem1e, sem1e.cov) ## no sig. diff. sem1e is more parsimonious (df=40 vs 46)

### Naturalised range size and expansion could be responding to the same unmeasured variable, rather than being in a causal relationship
submodels1e.cov1 <- '
yr1_1 ~ rsEURm2_t1
rsUSAm2_t1 ~ rsEURm2_t1 + rsEURm2_t2 + yr1_1 + growthTree + growthShrub
density_europe_t1 ~ rsEURm2_t1 + rsEURm2_t2 + citations_t1 + citations_t2 + growthTree + growthShrub
en_rel_t1 ~ rsUSAm2_t1 + rsUSAm2_t2 + yr1_1 + yr1_2 + density_europe_t1 + density_europe_t2 + citations_t1 + citations_t2 + growthTree + growthShrub + rsUSAm2_t1:density_europe_t1
expansion_tl ~ rsEURm2_t1 + rsEURm2_t2 + growthTree + growthShrub
citations_t2~~expansion_tl
citations_t1~~yr1_1
expansion_tl~~rsUSAm2_t1
'
# summary(sem2 <- sem(submodels2, dat), standardized=T, rsq=T, fit.measures=T) ## 
summary(sem1e.cov1 <- sem(submodels1e.cov1, dat, test="bollen.stine", se="boot"), standardized=T, rsq=T, fit.measures=T) ## DFI and RMSEA out of range. Chi-sq & other test now all within preferred values
write.csv(parameterEstimates(sem1e.cov1), file="semparams1e_cov1.csv", row.names=F) 

AIC(sem1e) ## -1525.667
AIC(sem1e.cov) ## -1787.048 Much lower
AIC(sem1e.cov1) ## -1523.31 worse than others, so causal relationship should remain

lavTestLRT(sem1e, sem1e.cov1) ## no sig. diff. Same df.


