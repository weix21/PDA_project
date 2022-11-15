## Xin and Chichun
## Model
## Date: 11/14/2022

# Load used packages
library(lme4)
library(MASS)
library(AER)

## Load Data
load("Processed_data_11_14.RData")

## Do Poisson regression
m1 <- glm(nCollection_Quarter ~ Species+Collection_Season+Region+Species:Collection_Season+Species:Region+Collection_Season:Region, family="poisson", data=df_2000)

## Do ANOVA to test the significance
m1_spec <- glm(nCollection_Quarter ~ Collection_Season+Region+Collection_Season:Region, family="poisson", data=df_2000)
m1_seas <- glm(nCollection_Quarter ~ Species+Region+Species:Region, family="poisson", data=df_2000)
m1_reg <- glm(nCollection_Quarter ~ Species+Collection_Season+Species:Collection_Season, family="poisson", data=df_2000)
m1_all <- glm(nCollection_Quarter ~ Species+Collection_Season+Region, family="poisson", data=df_2000)

## Test if there's overdispersion
dispersiontest(m1)

## Do quasi-Poisson regression
m2 <- glm(nCollection_Quarter ~ Species+Collection_Season+Region+Species:Collection_Season+Species:Region+Collection_Season:Region, family="quasipoisson", data=df_2000)

## Do Negative Binomial regression
m3 <- glm.nb(nCollection_Quarter ~ Species+Collection_Season+Region+Species:Collection_Season+Species:Region+Collection_Season:Region, data=df_2000)

## Do likelihood ratio test
pchisq(2 * (logLik(m3) - logLik(m1)), df = 1, lower.tail = FALSE)

# Check colinearlity
vif(m2)
vif(m3)