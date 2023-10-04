## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  echo=TRUE,
  collapse = TRUE,
  comment = "#>"
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

# Add treatment as a factor variable 
df$treatment <- factor(df$txCode, level=1:4,
          label=c("PPCT", "PCT", "ABCP", "BCP"))

# Select time points for aggregating data - here we use 1 week intervals
timepoints <- seq(1,72, by=1)

# Time points including zero
timepoints2 <- c(0, timepoints)

# Apply function
anova <- anova_data(timepoints=timepoints, timepoints2=timepoints2, ref.study=3,
                      df=df)

## -----------------------------------------------------------------------------
jagsdat <- fp_data(data=anova, trtnames=c("BCP", "ABCP", "PPCT", "PCT"),
                   polyorder=1)

## -----------------------------------------------------------------------------
inits <- fp_geninits(ns=4, polyorder=1, seed=890421)

## ---- results="hide"----------------------------------------------------------
# Set fractional polynomial power
jagsdat$P1 <- -0.5

# Model does not include (or monitor) rankings - can easily be done in R
jagsmod <- jags(data=jagsdat, inits=inits,
                parameters.to.save=c("d", "mu", "dev", "totresdev"),
                model.file=system.file("JAGSmodels", "FE_1st_order_model.jags", package="BristolTAG"),
                jags.seed=890421,
                n.iter=30000, n.burnin=2000, n.thin=1
)

## ---- eval=FALSE--------------------------------------------------------------
#  # Note that eval=FALSE so this command won't overwrite files document is rendered
#  saveRDS(jagsmod, file=paste0("Analysis/NS1/JAGS_Results/pfs_fpoly1_", jagsdat$P1))

## ---- results="hide"----------------------------------------------------------
modseq <- sequence_fpoly(jagsdat, powers=c(-1,-0.5))

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
summary(modseq)

## -----------------------------------------------------------------------------
hazardrats <- hrcalc(jagsmod, times=seq(1,60, length.out=100), eform=TRUE)

## -----------------------------------------------------------------------------
plot(hazardrats, reftrt = 1)

## ---- results="hide"----------------------------------------------------------
# Generate survival quantities
surv <- survcalc(jagsmod, refstudy=1)

## -----------------------------------------------------------------------------
head(surv$summary$S)

## -----------------------------------------------------------------------------
plot(surv, quantity="S", treats=c(3,4))

## -----------------------------------------------------------------------------
auc.surv <- auc(surv, quantity="S", tint=c(1,50))
print(auc.surv)

## ---- results="hide"----------------------------------------------------------
# Loop over studies
glist <- list()
for (i in 1:4) {
  surv <- survcalc(jagsmod, refstudy=i, n.mcmc=500)
  glist[[length(glist)+1]] <- plot(surv, quantity="S", treats=c(1,3))
}

## -----------------------------------------------------------------------------
ggpubr::ggarrange(glist[[1]], glist[[2]], glist[[3]], glist[[4]])

## -----------------------------------------------------------------------------
ipd1 <- readRDS(system.file("extdata", "PPCT PFS ITT PDL1 LT50pct KEYNOTE189.rds", package="BristolTAG"))
ipd2 <- readRDS(system.file("extdata", "PCT PFS ITT PDL1 LT50pct KEYNOTE189.rds", package="BristolTAG"))

## ---- warning=FALSE-----------------------------------------------------------
plot(surv, treats=3:4, kmdat = list("3"=ipd1, "4"=ipd2))

