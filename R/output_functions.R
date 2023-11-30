

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
#'
#'
#' @export
output_coda_fp <- function(jagsmod, format="mvn", refstudy="KEYNOTE189",
                           treatments=attr(jagsmod, "trtnames")) {
  if (!"rjags" %in% class(jagsmod)) {
    stop("jagsmod must be an object of class rjags")
  }

  sims.list <- jagsmod$BUGSoutput$sims.list

  trtnames <- attr(jagsmod, "trtnames")
  studynames <- attr(jagsmod, "studynames")

  mu <- sims.list$mu[,which(studynames %in% refstudy),]
  d <- sims.list$d

  # Get correct parameter estimates into matrix form
  mu <- apply(mu, MARGIN=2, cbind)

  # Loop over FP parameters
  dmat <- matrix(nrow=nrow(mu))
  cols <- vector()
  for (i in 1:dim(d)[3]) {
    dmat <- cbind(dmat, apply(d[,which(trtnames %in% treatments),i], MARGIN=2, cbind))
    cols <- append(cols, paste(trtnames[which(trtnames %in% treatments)], i, sep="_"))
  }
  dmat <- dmat[,-1]

  colnames(dmat) <- cols
  colnames(mu) <- paste(refstudy, 1:ncol(mu), sep="_")

  outmat <- cbind(mu, dmat)

  if (format=="mvn") {

    # Means
    param.means <- apply(outmat, MARGIN=2, mean)

    # Covariance - Should reference treatment d's be dropped?
    param.covar <- cov(outmat)

    # Param names
    temp <- paste0("d[", which(trtnames %in% treatments), ",")
    temp <- as.vector(t(sapply(temp, FUN=function(x) {paste0(x, 1:dim(d)[3], "]")})))
    param.symbol <- c(paste0("mu[", which(studynames %in% refstudy), ",", 1:ncol(mu), "]"),
                      temp)

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
#' @param quantiles
#'
#' @export
survquantile <- function(jagsmod, refstudy, quantile=c(0.5), n.mcmc=NULL,
                         interval=c(0.1, 300)) {

  trtnames <- attr(jagsmod, "trtnames")

  if (!"rjags" %in% class(jagsmod)) {
    stop("jagsmod must be an object of class rjags")
  }
  if (!all(c("d", "mu") %in% jagsmod$parameters.to.save)) {
    stop("d and mu must be monitored in jagsmod")
  }
  if (dim(jagsmod$BUGSoutput$sims.list$d)[3]==2) {
    polyorder <- 1
  } else if (dim(jagsmod$BUGSoutput$sims.list$d)[3]==3) {
    polyorder <- 2
  }

  # Speeds up computation
  if (!is.null(n.mcmc)) {
    mcmc.index <- sample(1:jagsmod$BUGSoutput$n.sims, size=n.mcmc)
    matsize <- n.mcmc
  } else {
    mcmc.index <- 1:jagsmod$BUGSoutput$n.sims
    matsize <- jagsmod$BUGSoutput$n.sims
  }

  jagsdat <- jagsmod$model$data()
  sims.list <- jagsmod$BUGSoutput$sims.list

  # Get index of reference study
  refstudy.ind <- which(attr(jagsmod, "studynames") %in% refstudy)

  reftrt <- jagsdat$t[refstudy.ind,1]

  haz.df <- data.frame()
  cumhaz.df <- data.frame()
  mort.df <- data.frame()
  S.df <- data.frame()

  mu <- sims.list$mu[, refstudy.ind, ]
  d  <- sims.list$d

  exponents <- jagsdat$P1
  if (!is.null(jagsdat$P2)) {
    exponents <- c(exponents, jagsdat$P2)
  }

  outlist <- list()
  for (k in 1:jagsdat$nt) {
    beta <- mu[mcmc.index,] + (d[mcmc.index,k,] - d[mcmc.index,reftrt,])

    outlist[[trtnames[k]]] <- apply(beta, MARGIN=1,
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
cumulative_hazard_function <- function(t, params, exponents) {
  integrate(hazard_function, params=params, exponents=exponents, lower = 0.1, upper = t)$value
}



#' Solve for cumulative hazard function at desired survival quantile
survquant_chunk <- function(params, exponents, quantile=0.5, interval=c(0.1, 1000)) {

  quantile <- -log(quantile)

  uniroot(function(t) {cumulative_hazard_function(t,
                                                  params=params,
                                                  exponents=exponents) - quantile},
          interval = interval
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
