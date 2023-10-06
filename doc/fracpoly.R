## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  echo=TRUE,
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup, results="hide", message=FALSE, warning=FALSE----------------------
library(BristolTAG)
library(tidyverse)
library(R2jags)
library(mcmcplots)
library(doBy)
library(survival)
library(RColorBrewer)

## -----------------------------------------------------------------------------
# Load the data to use
#df <- read.csv("Analysis/NS1/fracpoly/NS1 PFS IPD.csv")
df <- NS1_PFS # (example here loads it from package to allow vignette to run)

loglog_plot(df)

## -----------------------------------------------------------------------------
# Apply function
anova <- anova_data(timepoints=seq(0,72, by=1), 
                      df=df)

## -----------------------------------------------------------------------------
jagsdat <- fp_data(anova, trtnames=c("BCP", "ABCP", "PPCT", "PCT"))

## -----------------------------------------------------------------------------
inits <- fp_geninits(ns=4, polyorder=1, seed=890421)

## ---- results="hide"----------------------------------------------------------
modseq <- sequence_fpoly(jagsdat, powers=c(-1,-0.5), polyorder=1, inits=inits)

## ---- results="hide"----------------------------------------------------------
# Set fractional polynomial power
jagsdat$P1 <- -0.5

# Model does not include (or monitor) rankings - can easily be done in R
jagsmod <- jags(data=jagsdat, inits=inits,
                parameters.to.save=c("d", "mu", "dev", "totresdev"),
                model.file=system.file("JAGSmodels", "FE_1st_order_model.jags", package="BristolTAG"),
                jags.seed=890421,
                n.iter=15000, n.burnin=2000, n.thin=1
)

# Set attributes to match those of jagsdat
attr(jagsmod, "trtnames") <- attr(jagsdat, "trtnames")
attr(jagsmod, "studynames") <- attr(jagsdat, "studynames")
attr(jagsmod, "ipd") <- attr(jagsdat, "ipd")

## ---- eval=FALSE--------------------------------------------------------------
#  # Note that eval=FALSE so this command won't overwrite files document is rendered
#  saveRDS(jagsmod, file=paste0("Analysis/NS1/JAGS_Results/pfs_fpoly1_", jagsdat$P1))

## -----------------------------------------------------------------------------
summary(modseq)

## ---- eval=FALSE--------------------------------------------------------------
#  mcmcplot(jagsmod, parms = c("d"))

## -----------------------------------------------------------------------------
# Looking at model fit of a single model
devplot(jagsmod)

## -----------------------------------------------------------------------------
# Dev-dev plot for comparing models
# We'll compare the two models run previously using sequence_fpoly()
devplot(jagsmod1=modseq$`FP_-1`, jagsmod2=modseq$`FP_-0.5`)

## -----------------------------------------------------------------------------
hazardrats <- hrcalc(jagsmod, times=seq(1,60, length.out=100), eform=TRUE)

## -----------------------------------------------------------------------------
plot(hazardrats, reftrt = "BCP")

## ---- results="hide"----------------------------------------------------------
# Generate survival quantities
surv <- survcalc(jagsmod, refstudy="KEYNOTE-189")

## -----------------------------------------------------------------------------
head(surv$summary$S)

## -----------------------------------------------------------------------------
plot(surv, quantity="S", treats=c("ABCP", "PPCT"))

## -----------------------------------------------------------------------------
auc.surv <- auc(surv, quantity="S", tint=c(1,50))
print(auc.surv)

## ---- warning=FALSE-----------------------------------------------------------
plot(surv, treats=c("PCT", "PPCT"), overlay.km = TRUE)

