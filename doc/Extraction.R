## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, results="hide", warning=FALSE, message=FALSE----------------------
library(BristolTAG)
library(survival)
library(readr)
library(tidyverse)

## ---- results="hide", warning=FALSE, message=FALSE----------------------------
#### Load the .csv file ####

# In this example we read from the data within the package
km1 <- read_csv(system.file("extdata", "IMpower150_ABCP_OS_TCIC123.csv", package="BristolTAG"), 
                col_names=FALSE)

# Clean the data
km1 <- clean.digitised(km1)

## ---- results="hide", warning=FALSE, message=FALSE----------------------------
ipd1 <- guyot.method(km1$time, km1$survival, 
                     t=seq(0,52,4),
                     r=c(143,126,109,94,85,72,62,54,47,28,13,5,1,0),
                     tot.events=NA
                     )

## ---- results="hide", warning=FALSE, message=FALSE----------------------------
#### TCIC123 ####

# We already have ipd1


#### TCIC0 ####
km2 <- read_csv(system.file("extdata", "IMpower150_ABCP_OS_TCIC0.csv", package="BristolTAG"), 
                col_names=FALSE)
km2 <- clean.digitised(km2)
ipd2 <- guyot.method(km2$time, km2$survival, 
                     t=seq(0,44,4),
                     r=c(120,105,91,75,60,48,39,30,25,14,6,2),
                     tot.events=NA
                     )

#### TCIC3 ####
km3 <- read_csv(system.file("extdata", "IMpower150_ABCP_OS_TCIC3.csv", package="BristolTAG"), 
                col_names=FALSE)
km3 <- clean.digitised(km3)
ipd3 <- guyot.method(km3$time, km3$survival, 
                     t=seq(0,52,4),
                     r=c(50,44,36,34,31,29,27,25,22,13,4,2,1,0),
                     tot.events=NA
                     )

## ---- results="hide", warning=FALSE, message=FALSE----------------------------
ipd12 <- combine.ipd.guyot(ipd1, ipd2)

## ---- results="hide", warning=FALSE, message=FALSE----------------------------
ipd.subgroup <- sub.ipd.guyot(ipd12, ipd3)

## -----------------------------------------------------------------------------
survival_fit <- survfit(Surv(ipd.subgroup$patient$time, ipd.subgroup$patient$event) ~ 1)
plot(survival_fit)

