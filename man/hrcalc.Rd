% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fracpoly_functions.R
\name{hrcalc}
\alias{hrcalc}
\alias{plot.hazard.ratios}
\title{Calculate time-varying hazard ratios (HR)}
\usage{
hrcalc(jagsmod, times = seq(1, 60, length.out = 100), eform = TRUE)

\method{plot}{hazard.ratios}(hr, reftrt, treats = attr(hr, "trtnames"))
}
\arguments{
\item{jagsmod}{an object of class `"rjags"` that contains the results of a
fractional polynomial model.}

\item{times}{a vector of times at which to estimate time-varying HRs}

\item{eform}{whether results should be on log (`FALSE`) or exponentiated onto
natural (`TRUE`) scale.}

\item{hr}{An object of class `"hazard.ratios"` containing the output of `hrcalc()`}

\item{reftrt}{The name of the reference treatment against which HRs should be plotted. Must match
a name in `attr(hr, "trtnames")`}

\item{treats}{treatments to plot. Default plots them all}
}
\description{
Calculate time-varying hazard ratios (HR)
}
\section{Functions}{
\itemize{
\item \code{plot(hazard.ratios)}: Plot time-varying HRs

}}
