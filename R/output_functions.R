

#####################################################################
###### Functions for outputting results from models #######
#####################################################################


#' Output results from fractional polynomial NMA for use in health economic model
#'
#' @inheritParams hrcalc
#' @param format Can take `"mvn"` to output parameters for multivariate normal
#' distribution, or `"coda"` for a full set of MCMC samples (given as a matrix)
#' @param params A vector of parameter names corresponding to those in the model
#' to extract values for.
#' @param refstudy Name of the study to use as reference when making predictions.
#' `refstudy` must either match a name in `attr(jagsmod, "studynames")` or can be
#' left as `NULL` (the default) for parameters for the reference curve not to be
#' included in the coda.
#'
#'
#' @export
output_coda_fp <- function(jagsmod, format="mvn", refstudy=NULL,
                           treatments=attr(jagsmod, "trtnames")) {
  if (!"rjags" %in% class(jagsmod)) {
    stop("jagsmod must be an object of class rjags")
  }

  sims.list <- jagsmod$BUGSoutput$sims.list

  trtnames <- attr(jagsmod, "trtnames")
  studynames <- attr(jagsmod, "studynames")

  if (!is.null(refstudy)) {
    mu <- sims.list$mu[,which(studynames %in% refstudy),]

    # Get correct parameter estimates into matrix form
    mu <- apply(mu, MARGIN=2, cbind)
  }

  d <- sims.list$d

  # Loop over FP parameters
  dmat <- matrix(nrow=dim(d)[1])
  cols <- vector()
  for (i in 1:dim(d)[3]) {
    dmat <- cbind(dmat, apply(d[,which(trtnames %in% treatments),i], MARGIN=2, cbind))
    cols <- append(cols, paste(trtnames[which(trtnames %in% treatments)], i, sep="_"))
  }
  dmat <- dmat[,-1]

  colnames(dmat) <- cols

  if (!is.null(refstudy)) {
    colnames(mu) <- paste(refstudy, 1:ncol(mu), sep="_")
    outmat <- cbind(mu, dmat)
  } else {
    outmat <- dmat
  }

  if (format=="mvn") {

    # Means
    param.means <- apply(outmat, MARGIN=2, mean)

    # Covariance - Should reference treatment d's be dropped?
    param.covar <- cov(outmat)

    # Param names
    temp <- paste0("d[", which(trtnames %in% treatments), ",")
    temp <- as.vector(t(sapply(temp, FUN=function(x) {paste0(x, 1:dim(d)[3], "]")})))

    if (!is.null(refstudy)) {
      param.symbol <- c(paste0("mu[", which(studynames %in% refstudy), ",", 1:ncol(mu), "]"),
                        temp)
    } else {
      param.symbol <- temp
    }

    out <- list(mean=param.means, covar=param.covar, params=param.symbol)

  } else if (format=="coda") {

    out <- outmat

  } else {
    stop("format must be 'mvn' or 'coda'")
  }

  return(out)
}




#' Estimate survival quantiles (e.g. median survival)
#'
#' @inheritParams survcalc
#' @param quantile A quantile (between 0 and 1) at which to estimate the survival
#' time - e.g. for median survival `quantile=0.5`
#' @param treatments A character vector of treatments names for which to estimate
#' survival quantiles
#' @param interval Two numbers representing the limits over which to evaluate the
#' solver that calculates the survival quantile. Note that this **must** be wide
#' enough to allow for all MCMC survival quantiles or the function will throw
#' an error. However, a more narrow interval will speed up the solver.
#'
#' @export
survquantile <- function(jagsmod, refstudy, refmod=jagsmod,
                         quantile=0.5, n.mcmc=NULL,
                         interval=c(0.1, 300),
                         treatments=attr(jagsmod, "trtnames")) {

  if (!"rjags" %in% class(jagsmod)) {
    stop("jagsmod must be an object of class rjags")
  }
  if (!"rjags" %in% class(refmod)) {
    stop("refmod must be an object of class rjags")
  }
  if (!all(c("d", "mu") %in% jagsmod$parameters.to.save)) {
    stop("d and mu must be monitored in jagsmod")
  }
  if (dim(jagsmod$BUGSoutput$sims.list$d)[3]==2) {
    polyorder <- 1
  } else if (dim(jagsmod$BUGSoutput$sims.list$d)[3]==3) {
    polyorder <- 2
  }

  jagsdat <- jagsmod$model$data()
  sims.list <- jagsmod$BUGSoutput$sims.list

  # Check treatments are all in jagsmod
  if (!all(treatments %in% attr(jagsmod, "trtnames"))) {
    stop("Named 'treatments' do not match those in jagsmod")
  }
  trt.ind <- which(attr(jagsmod, "trtnames") %in% treatments)

  # Speeds up computation
  # And ensures the largest MCMC object samples are fully used
  if (!is.null(n.mcmc)) {
    mcmc.index <- sample(1:jagsmod$BUGSoutput$n.sims, size=n.mcmc)
    matsize <- n.mcmc

    mcmc.index.jags <- mcmc.index
    mcmc.index.ref <- mcmc.index

  } else {
    mcmc.index.jags <- 1:jagsmod$BUGSoutput$n.sims
    mcmc.index.ref <- 1:refmod$BUGSoutput$n.sims

    if (jagsmod$BUGSoutput$n.sims > refmod$BUGSoutput$n.sims) {
      mcmc.index.ref <- c(mcmc.index.ref,
                          sample(1:refmod$BUGSoutput$n.sims,
                                 size=jagsmod$BUGSoutput$n.sims - refmod$BUGSoutput$n.sims)
      )
      matsize <- jagsmod$BUGSoutput$n.sims
    } else if (jagsmod$BUGSoutput$n.sims < refmod$BUGSoutput$n.sims) {
      mcmc.index.jags <- c(mcmc.index.jags,
                           sample(1:jagsmod$BUGSoutput$n.sims,
                                  size=refmod$BUGSoutput$n.sims - jagsmod$BUGSoutput$n.sims)
      )
      matsize <- refmod$BUGSoutput$n.sims
    } else {
      matsize <- jagsmod$BUGSoutput$n.sims
    }
    n.sims <- max(jagsmod$BUGSoutput$n.sims, refmod$BUGSoutput$n.sims)
  }

  # Get index of reference study
  if (!identical(jagsmod, refmod)) {
    message("Model for reference curve is different to the treatment effect model")
  }

  # Get index of reference study
  refstudy.ind <- which(attr(refmod, "studynames") %in% refstudy)
  if (length(refstudy.ind)==0) {
    if (dim(refmod$BUGSoutput$sims.list$mu)[2]==1) {
      refstudy.ind <- 1
    } else {
      stop("refstudy is not a named study in attr(refmod, 'studynames')")
    }
  }

  reftrt <- jagsdat$t[refstudy.ind,1]

  haz.df <- data.frame()
  cumhaz.df <- data.frame()
  mort.df <- data.frame()
  S.df <- data.frame()

  mu <- refmod$BUGSoutput$sims.list$mu[, refstudy.ind, ]
  d  <- sims.list$d

  exponents <- jagsdat$P1
  if (!is.null(jagsdat$P2)) {
    exponents <- c(exponents, jagsdat$P2)
  }

  outlist <- list()
  for (k in seq_along(trt.ind)) {
    beta <- mu[mcmc.index.ref,] +
      (d[mcmc.index.jags,trt.ind[k],] - d[mcmc.index.jags,reftrt,])

    outlist[[treatments[k]]] <- apply(beta, MARGIN=1,
                                    FUN=function(x) {
                                      survquant_chunk(params=x,
                                                      exponents=exponents,
                                                      quantile=quantile,
                                                      interval=interval)
                                    })
  }
  class(outlist) <- "survquantile"
  return(outlist)
}


#' Calculate hazard function
hazard_function <- function(t, params, exponents) {
  exp(get_fp(t, params=params, exponents=exponents))
}


#' Calculate cumulative hazard function
#'
#' @param method Can take either `"stats"` or `"pracma"`
cumulative_hazard_function <- function(t, lower=0.1, params, exponents,
                                       method="pracma") {

  if (method=="stats") {
    integrate(hazard_function, params=params, exponents=exponents,
              rel.tol=.Machine$double.eps^0.1,
              lower = lower, upper = t)$value
  } else if (method=="pracma") {
    pracma::quadgk(hazard_function, a = lower, b = t,
                 params=params, exponents=exponents,
                 tol=.Machine$double.eps^0.1)
  }
}



#' Solve for cumulative hazard function at desired survival quantile
survquant_chunk <- function(params, exponents, quantile=0.5, interval=c(0.1, 1000), ...) {

  quantile <- -log(quantile)

  uniroot(function(t) {cumulative_hazard_function(t,
                                                  params=params,
                                                  exponents=exponents) - quantile},
          interval = interval,
          ...
          )$root
}




#' Prints survival quantiles
#'
#' @param x an object of class `"survquantile"`
#' @param ... Arguments for `quantile()`
#'
#' @export
print.survquantile <- function(x, ...) {
  outmat <- matrix(ncol=3, nrow=length(x))
  for (i in seq_along(x)) {
    outmat[i,] <- quantile(x[[i]], probs=c(0.025, 0.5, 0.975), ...)
  }
  df <- data.frame(names(x))
  df <- cbind(df, outmat)
  names(df) <- c("Treatment", "2.5%", "50%", "97.5%")
  return(df)
}
