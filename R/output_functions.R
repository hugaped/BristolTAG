

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
