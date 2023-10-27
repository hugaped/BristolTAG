# Extraction functions



#' Go from KM curve data from published figures into individual-level data
#'
#' This function is a simple adaptation of code given by Guyot et al., taken from https://github.com/certara/survivalnma
#'
#' @param time X axis coordinates taken from the digitised Kaplan-Meier curve
#' @param survival Y axis coordinates taken from the digitised Kaplan-Meier curve; same length as `x`
#' @param tint times at which number at risk `rint` is given
#' @param rint number at risk; same length as `tint`
#' @param tot.events total number of events reported (integer) (optional - can be left as `NA`)
#' @export
#'
#' @examples
#' #Make up some data:
#' ipd_curve1 <- guyot.method(time = seq(0, 100), survival = 1-pexp(seq(0, 100), rate = 1/50),
#'              tint=c(0, 10, 50), rint = c(1000, 800, 250))
#' ipd_curve2 <- guyot.method(time = seq(0, 100), survival = 1-pexp(seq(0, 100), rate = 1/30),
#'              tint=c(0, 10, 50), rint = c(1000, 700, 100))
#' library(survival)
#' ipd_curve1$patient$treatment <- "active"
#' ipd_curve2$patient$treatment <- "chemo"
#' ipd <- rbind(ipd_curve1$patient, ipd_curve2$patient)
#' survival_fit <- survfit(Surv(time, event) ~ treatment, data = ipd)
#' survival_fit
#' plot(survival_fit, col = c("blue", "red"))
#'
#' @references Guyot, Patricia, AE Ades, Mario JNM Ouwens, and Nicky J. Welton.
#'             “Enhanced Secondary Analysis of Survival Data: Reconstructing the
#'             Data from Published Kaplan-Meier Survival Curves.”
#'             BMC Medical Research Methodology 12, no. 1 (February 1, 2012): 9.
#'             https://doi.org/10.1186/1471-2288-12-9.
guyot.method <- function(time, survival, tint, rint,
                         tot.events=NA) {

  x <- time
  y <- survival
  t <- tint
  r <- rint

  if(max(y) > 1)
    stop("Survival (y coordinates) can't be greater than 1. Please scale your y to be in [0, 1] range.")
  if(min(x) > 0){
    message("Adding coordinate (x=0, y=1) to your data.")
    x <- c(0, x)
    y <- c(1, y)
  }
  if(min(t) > 0) {
    warning("Number at risk at time 0 is needed. The calculation might not work correctly otherwise.")
  }
  if (length(t)!=length(r)) {
   stop("length(tint) should equal length(rint)")
  }
  if(max(t) > max(x)) {
    warning("Times for number at risk (t) where t > x (x = coordinates of the curve) are not used.")
    r <- r[t < max(x)]
    t <- t[t < max(x)]
  }
  # warn <- FALSE
  # for (i in 1:(length(x)-1)) {
  #
  #   if (x[i]==x[i+1]) {
  #     warn <- TRUE
  #     x[i+1] <- x[i+1]+0.0000001
  #   }
  # }
  # if (warn==TRUE) {
  #   warning("Adding 0.0000001 to identical time values in t")
  # }

  # Parameter names used by the script below
  t.S<-x; S<-y; n.risk<-r; t.risk<-t

  # Need to grab `lower` and `upper`, to "match" vector x against vector t
  # See Guyot Table 2 for visual explanation of this
  lower <- sapply(t, function(current_t) min(which(x >= current_t)))
  upper <- sapply(c(t[-1], Inf), function(current_t) max(which(x < current_t)))

  n.int<-length(n.risk)
  n.t<- upper[n.int]
  #Initialise vectors
  n.censor<- rep(0,(n.int-1))
  n.hat<-rep(n.risk[1]+1,n.t)
  cen<-rep(0,n.t)
  d<-rep(0,n.t)
  KM.hat<-rep(1,n.t)
  last.i<-rep(1,n.int)
  sumdL<-0
  if (n.int > 1){
    #Time intervals 1,...,(n.int-1)
    for (i in 1:(n.int-1)){
    #for (i in 1:16){
      #First approximation of no. censored on interval i
      n.censor[i]<- round(n.risk[i]*S[lower[i+1]]/S[lower[i]]- n.risk[i+1])
      #Adjust tot. no. censored until n.hat = n.risk at start of interval (i+1)
      while((n.hat[lower[i+1]]>n.risk[i+1])||((n.hat[lower[i+1]]<n.risk[i+1])&&(n.censor[i]>0))){
        # if ((n.hat[lower[i+1]]>n.risk[i+1])==FALSE) {
        #   stop("ME")
        # }

        if (n.censor[i]<=0){
          cen[lower[i]:upper[i]]<-0
          n.censor[i]<-0
        }
        if (n.censor[i]>0){
          cen.t<-rep(0,n.censor[i])
          for (j in 1:n.censor[i]){
            cen.t[j]<- t.S[lower[i]] +
              j*(t.S[lower[(i+1)]]-t.S[lower[i]])/(n.censor[i]+1)
          }
          #Distribute censored observations evenly over time. Find no. censored on each time interval.
          cen[lower[i]:upper[i]]<-graphics::hist(cen.t,breaks=t.S[lower[i]:lower[(i+1)]],
                                                 plot=F)$counts
        }
        #Find no. events and no. at risk on each interval to agree with K-M estimates read from curves
        n.hat[lower[i]]<-n.risk[i]
        last<-last.i[i]


        temp.low <- min(c(lower[i], upper[i]))
        temp.high <- max(c(lower[i], upper[i]))
        for (k in temp.low:temp.high){
        #for (k in lower[i]:upper[i]){
          if (i==1 & k==lower[i]){
            d[k]<-0
            KM.hat[k]<-1
          }
          else {
            d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))
            KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
          }
          n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
          if (d[k] != 0) last<-k
        }
        n.censor[i]<- n.censor[i]+(n.hat[lower[i+1]]-n.risk[i+1])
      }
      if (n.hat[lower[i+1]]<n.risk[i+1]) {
        n.risk[i+1]<-n.hat[lower[i+1]]
        }
      last.i[(i+1)]<-last
    }
  }
  #Time interval n.int.
  if (n.int>1){
    #Assume same censor rate as average over previous time intervals.
    n.censor[n.int]<- min(round(sum(n.censor[1:(n.int-1)])*(t.S[upper[n.int]]-
                                                              t.S[lower[n.int]])/(t.S[upper[(n.int-1)]]-t.S[lower[1]])), n.risk[n.int])
  }
  if (n.int==1){n.censor[n.int]<-0}
  if (n.censor[n.int] <= 0){
    cen[lower[n.int]:(upper[n.int]-1)]<-0
    n.censor[n.int]<-0
  }
  if (n.censor[n.int]>0){
    cen.t<-rep(0,n.censor[n.int])
    for (j in 1:n.censor[n.int]){
      cen.t[j]<- t.S[lower[n.int]] +
        j*(t.S[upper[n.int]]-t.S[lower[n.int]])/(n.censor[n.int]+1)
    }
    cen[lower[n.int]:(upper[n.int]-1)]<-graphics::hist(cen.t,breaks=t.S[lower[n.int]:upper[n.int]],
                                                       plot=F)$counts
  }
  #Find no. events and no. at risk on each interval to agree with K-M estimates read from curves
  n.hat[lower[n.int]]<-n.risk[n.int]
  last<-last.i[n.int]
  for (k in lower[n.int]:upper[n.int]){
    if(KM.hat[last] !=0){
      d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))} else {d[k]<-0}
    KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
    n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
    #No. at risk cannot be negative
    if (n.hat[k+1] < 0) {
      n.hat[k+1]<-0
      cen[k]<-n.hat[k] - d[k]
    }
    if (d[k] != 0) last<-k
  }
  #If total no. of events reported, adjust no. censored so that total no. of events agrees.
  if (!is.na(tot.events)){
    if (n.int>1){
      sumdL<-sum(d[1:upper[(n.int-1)]])
      #If total no. events already too big, then set events and censoring = 0 on all further time intervals
      if (sumdL >= tot.events){
        d[lower[n.int]:upper[n.int]]<- rep(0,(upper[n.int]-lower[n.int]+1))
        cen[lower[n.int]:(upper[n.int])]<- rep(0,(upper[n.int]-lower[n.int]+1))
        n.hat[(lower[n.int]+1):(upper[n.int]+1)]<- rep(n.risk[n.int],(upper[n.int]+1-lower[n.int]))
      }
    }
    #Otherwise adjust no. censored to give correct total no. events
    if ((sumdL < tot.events)|| (n.int==1)){
      sumd<-sum(d[1:upper[n.int]])
      while ((sumd > tot.events)||((sumd< tot.events)&&(n.censor[n.int]>0))){
        n.censor[n.int]<- n.censor[n.int] + (sumd - tot.events)
        if (n.censor[n.int]<=0){
          cen[lower[n.int]:(upper[n.int]-1)]<-0
          n.censor[n.int]<-0
        }
        if (n.censor[n.int]>0){
          cen.t<-rep(0,n.censor[n.int])
          for (j in 1:n.censor[n.int]){
            cen.t[j]<- t.S[lower[n.int]] +
              j*(t.S[upper[n.int]]-t.S[lower[n.int]])/(n.censor[n.int]+1)
          }
          cen[lower[n.int]:(upper[n.int]-1)]<-graphics::hist(cen.t,breaks=t.S[lower[n.int]:upper[n.int]],
                                                             plot=F)$counts
        }
        n.hat[lower[n.int]]<-n.risk[n.int]
        last<-last.i[n.int]
        for (k in lower[n.int]:upper[n.int]){
          d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))
          KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
          if (k != upper[n.int]){
            n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
            #No. at risk cannot be negative
            if (n.hat[k+1] < 0) {
              n.hat[k+1]<-0
              cen[k]<-n.hat[k] - d[k]
            }
          }
          if (d[k] != 0) last<-k
        }
        sumd<- sum(d[1:upper[n.int]])
      }
    }
  }
  wt1 <- matrix(c(t.S,n.hat[1:n.t],d,cen),ncol=4,byrow=F)

  ### Now form IPD ###
  #Initialise vectors
  t.IPD<-rep(t.S[n.t],n.risk[1])
  event.IPD<-rep(0,n.risk[1])
  #Write event time and event indicator (=1) for each event, as separate row in t.IPD and event.IPD
  k=1
  for (j in 1:n.t){
    if(d[j]!=0){
      t.IPD[k:(k+d[j]-1)]<- rep(t.S[j],abs(d[j]))
      event.IPD[k:(k+d[j]-1)]<- rep(1,abs(d[j]))
      k<-k+d[j]
    }
  }
  #Write censor time and event indicator (=0) for each censor, as separate row in t.IPD and event.IPD
  cen <- abs(cen)
  for (j in 1:(n.t-1)){
    if(cen[j]!=0){
      t.IPD[k:(k+cen[j]-1)]<- rep(((t.S[j]+t.S[j+1])/2),cen[j])
      event.IPD[k:(k+cen[j]-1)]<- rep(0,cen[j])
      k<-k+cen[j]
    }
  }

  #Output IPD
  IPD<-matrix(c(t.IPD,event.IPD),ncol=2,byrow=F)

  out <- list("curve"=as.data.frame(wt1),
              "patient"=data.frame(time = t.IPD, event = event.IPD))

  class(out) <- "ipd.guyot"
  return(out)
}







#' Take a subset of digitised IPD
#'
#' @param x An object of class `"ipd.guyot"` generated from digitised IPD containing
#' the larger set of IPD
#' @param y An object of class `"ipd.guyot"` generated from digitised IPD that should be
#' removed from `x`
#'
#' @details
#' Note that `x` must include patients in `y`. The function then removes IPD patients specified in `y`
#' from IPD in `x` using nearest-neighbour.
#'
#' @export
sub.ipd.guyot <- function(x, y) {

  if (nrow(x$patient)<nrow(y$patient)) {
    stop("y is larger than x - y must be a subset of x")
  }

  out <- x$patient

  for (i in seq_along(y$patient[,1])) {
    # Subset for event
    sub <- out[out$event==y$patient$event[i],]

    sub$tdif <- abs(sub$time - y$patient$time[i])
    sub <- arrange(sub, tdif)

    tmatch <- sub$time[1]

    drop <- which(out$time==tmatch & out$event==y$patient$event[i])

    #print(tmatch)

    out <- out[-drop[1],]
  }

  out <- list(patient=out)
  class(out) <- "ipd.guyot"

  return(out)
}



#' Combines digitised IPD
#'
#' @param x An object of class `"ipd.guyot"` generated from digitised IPD
#' @param y An object of class `"ipd.guyot"` generated from digitised IPD
#'
#' @return Returns a list containing a data frame of IPD that combines data from `x$patient` and
#' `y$patient`
#'
#' @export
combine.ipd.guyot <- function(x, y) {

  # Combining the patient data frame
  pat <- arrange(rbind(x$patient, y$patient), time, desc(event))

  out <- list("patient"=pat)
  class(out) <- "ipd.guyot"

  return(out)
}






#' Cleans digitised KM data
#'
#' @param km a data frame with two columns that correspond to locations
#' on the x- (1st column) and y- (2nd column) axes of a digitised Kaplan Meier
#' curve
#'
#' @export
clean.digitised <- function(km) {

  km <- km[order(km[[1]]),] # order the data if it has not need ordered in the webplot tool

  digitised <- data.frame(1:nrow(km), km[,1:2])

  names(digitised) <- c("k","time", "survival")

  ## If survival probability is greater than 1 then divided by 100 ##
  if(digitised[1,"survival"] > 10) {
    digitised[,"survival"] <- digitised[,"survival"]/100
  }
  head(digitised)


  ## clean up the data, time is increasing, survival is decreasing ##

  digitised$time[1] <- ifelse(digitised$time[1] < 0 | digitised$time[1] > 10, 0 , digitised$time[1])
  digitised$survival[1] <- ifelse(digitised$survival[1] > 1 | digitised$survival[1] < 0, 1 , digitised$survival[1])

  for(i in 2:nrow(digitised)){
    if(digitised[i,"time"] < digitised[i-1,"time"])
    {
      digitised[i,"time"] <- digitised[i-1,"time"]
    }
    else if(digitised[i,"survival"] > digitised[i-1,"survival"] | digitised[i,"survival"] < 0)
    {
      digitised[i,"survival"] <- digitised[i-1,"survival"]
    }
  }

  return(digitised)
}





#' Checks digitised KM data
#'
#' Checks a data frame of digitised data to ensure survival and times are
#' correctly specified (i.e. survival at each time-point is <= survival at previous
#' time-point and time >= previous time)
#'
#' @param digitised A data frame containing variables `time` and `survival`, ordered by
#' `time`
#'
check.digitised <- function(digitised) {
  for (i in 2:nrow(digitised)) {
    if (digitised$time[i]<digitised$time[i-1]) {
      warning(paste0("Row ", digitised$k[i], ": ", digitised$time[i], " < ", digitised$time[i-1]))
      digitised$time[i] <- digitised$time[i-1]
    }
  }

  for (i in 2:nrow(digitised)) {
    if (digitised$survival[i]>digitised$survival[i-1]) {
      warning(paste0("Row ", digitised$k[i], ": ", digitised$survival[i], " > ", digitised$survival[i-1]))
      digitised$survival[i] <- digitised$survival[i-1]
    }
  }
  return(digitised)
}








