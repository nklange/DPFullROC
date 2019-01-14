#' Dual process MLE: Make start parameters
#'
#' This function estimates start parameters for the DP univariate model fit, with starting parameters for the
#' first 40% of iterations sampled from around the observed data, and 60% sampled from a uniform distribution.
#'
#' @author Nicholas Lange, \email{lange.nk@gmail.com}
#' @param freqFA frequency of FA respones in category bins, Lure to Target order
#' @param freqHit frequency of Hit responses in category bins, Lure to Target order
#' @param iterations the number of runs
#' @return The function returns a dataframe:
#' \item{(parameters)}{Starting parameters for all parameters that versions of the DP Recog-Source model could vary to date}
#' @export

makeMLEstartparameters<-function(freqFA,freqHit,iterations){
  nBins <- length(freqHit)

  ratingMid <- length(freqFA)/2


  p_hit <- (sum(freqHit[(ratingMid+1):nBins]) + 0.5)/(sum(freqHit) + 1)
  p_fa <- (sum(freqFA[(ratingMid+1):nBins]) + 0.5)/(sum(freqFA) + 1)

  #mu
  e.dpri <- stats::qnorm(p_hit) - stats::qnorm(p_fa)


  # criteria
  ProbFA <- vector()
  ProbH <- vector()

  # calc response probabilities from data
  for (i in c(1:nBins)) {
    ProbFA[i] <-  (freqFA[i] + 0.5)/(sum(freqFA) + 1)

    ProbH[i] <-  (freqHit[i] + 0.5)/(sum(freqHit) + 1)
  }

  # Get cumulative probabilities
  CumProbH <- rev(cumsum(rev(ProbH)))
  CumProbH <- CumProbH[2:nBins]
  CumProbFA <- rev(cumsum(rev(ProbFA)))
  CumProbFA <- CumProbFA[2:nBins]

  for (i in c(1:length(CumProbFA))){
    if (CumProbFA[i] >= 1) {CumProbFA[i] <- .99}
    if (CumProbH[i] >= 1) {CumProbH[i] <- .99}
  }

  # Recognition Criteria
  e.TMPCs <- c()
  for (i in c(1:(nBins-1))){
    if (i < 2){
      e.TMPCs[i] <- stats::qnorm(1-CumProbFA[i], sd = 1) # c values are relative to new
    }

    else {
      e.TMPCs[i] <- stats::qnorm(1-CumProbFA[i], sd = 1) - stats::qnorm(1-CumProbFA[i-1], sd = 1) # diff - c2 will = c1 + dc2. dc2 = c2 - c1

    }
  }

  # Sigf Old

  zCumProbN <- stats::qnorm(CumProbFA)
  zCumProbO <- stats::qnorm(CumProbH)

  zROCSlope <- round(stats::lm(zCumProbO ~ zCumProbN)$coef[2], digits = 2)

  e.sd_target   <- 1/zROCSlope # if unequal variance


  close<-round((2/5)*iterations) #some runs drawn from close around observed mean
  random<-round((3/5)*iterations) #some runs with starting parameters drawn at relatively random

  dpri<-log(c(truncnorm::rtruncnorm(close,mean=e.dpri,sd=0.3,a=0.1), stats::runif(random,min=0.1,max=5)))
  sd_target<-log(c(truncnorm::rtruncnorm(close,mean=e.sd_target,sd=0.3,a=0), stats::runif(random,min=1,max=2)))





  Cs <- data.frame(matrix(nrow = iterations, ncol = nBins-1))
  tempCsLabels<- NULL
  for (i in c(1:(nBins-1))){
    tempCsLabels[i] <- paste0("Cs",i)
  }
  colnames(Cs) <- tempCsLabels
  for (i in c(1:(nBins-1))){
    if (i < 2){
      Cs[,i] <- c(stats::rnorm(close,mean=e.TMPCs[i],sd=.02), stats::runif(random,min=-3,max=0))
    }
    else{
      Cs[,i]<-log(c(truncnorm::rtruncnorm(close,mean=e.TMPCs[i],sd=.02,a=0),stats::runif(random,min=0,max=1))) # Cs2 - Csi formally dCs(i), i.e. increase on top of Cs(i-1)
    }
  }

  # Recollection parameters not estimated from the data!
  RTraw <- c(truncnorm::rtruncnorm(close,mean=0.2,sd=0.2,a=0, b = 1), stats::runif(random,min=0,max=1))
  RLraw <- c(truncnorm::rtruncnorm(close,mean=0.2,sd=0.2,a=0, b = 1), stats::runif(random,min=0,max=1))

  rt <- log(RTraw/(1-RTraw))
  rl <- log(RLraw/(1-RLraw))


  startparss<-cbind(rt, rl, dpri, sd_target,
                    Cs)
  return(startparss)
}
