
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BristolTAG

<!-- badges: start -->
<!-- badges: end -->

The goal of BristolTAG is to provide a set of R functions for running
various evidence synthesis analyses needed by the University of
Bristol’s Technology Assessment Group. At the moment this just involves
running fractional polynomial Network Meta-Analyses of survival data…but
more may be added in the future.

Note that the fractional polynomial models are based on the work of
Jeroen Jansen (Network meta-analysis of survival data with fractional
polynomials; BMC Medical Research Methodology; 2011, 11:61) and BUGS
code has been adapted from Freeman, Cooper, Sutton, Crowther, Carpenter
and Hawkins (Challenges of modelling approaches for network
meta-analysis of time-to-event outcomes in the presence of
non-proportional hazards to aid decision making: Application to a
melanoma network; Statistical Methods in Medical Research; 2022, 31(5)).
I also used some bit from Roche’s `gemtcPlus` package
(<https://github.com/Roche/gemtcPlus>).

## Installation

You can install the development version of BristolTAG from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hugaped/BristolTAG")
```

## Getting started

For getting started with fractional polynomial NMA of survival data,
check out the `fracpoly` vignette.

## Citing the package

You probably won’t use this package if you’re not part of the Bristol
TAG, but if you do then please do cite it. We’re academics so this helps
us get that extra snippet of fame and glory (if not money). There’s no
publication, and we’re very unlikely to put the package on CRAN, but
just write somewhere that you used it and include a reference like this:

*BristolTAG R package; Hugo Pedder; Bristol Technology Assessment Group,
University of Bristol; 2023. Available at:
<https://github.com/hugaped/BristolTAG>.*
