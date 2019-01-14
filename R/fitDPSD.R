#' Estimation of recollection and familiarity by fitting recognition memory data with the Dual Process Signal Detection (DPSD) model
#'
#' This function allows to estimate recollection and familiarity for recognition memory data by fitting data to the DPSD model.
#' The optimization is attempted by minimizing the summed log-likelihood using the Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm
#' in \code{\link{optim}}.
#' The function uses random start values on each iteration in order to find the set of parameters,
#' which fit the data best by returning the values with the lowest negative log likelihood.
#' Recollection and Familiarity for the lure distribution (new items) are set to 0.
#' Optional arguments in the function allow the user to specify an equal-variance model.
#' Recollection is bounded to be between 0 and 1, Familiarity and the standard deviation of the target distribution to be positive.
#' Criteria are ordered.
#'
#' @author Nicholas Lange, \email{lange.nk@gmail.com}
#' @param falseAlarms A vector containing the number of false alarms per recognition category rating.
#' @param hit A vector containing the hit number of hits per recognition category rating.
#' @param iterations A numeric value specifying the number of iterations. Default is set to 200.
#' @param eqVar A boolean value specifying if the standard deviation of the target distribution is equal to that of the lure distribution (i.e. = 1) (TRUE) or estimated separately (FALSE). Default is set to TRUE.
#' @return The function returns a dataframe with components:
#' \item{(parameters)}{The estimated parameters (recollection_target, recollection_lure = 0, familiarity, sd_target, criteria) for the iteration with the lowest SumSquareError}
#' \item{SSE}{Minimum sum square error}
#' @references Yonelinas, A. P. (1999). The Contribution of Recollection and Familiarity to Recognition and Source-Memory Judgments: A Formal Dual-Process Model and an Analysis of Receiver Operating Characteristics. Journal of Experimental Psychology: Learning, Memory, and Cognition, 25(6), 1415 - 1434. http://doi.org/10.1037//0278-7393.25.6.1415
#' @keywords ROC recollection familiarity DPSD
#' @export

fitDPSD <- function(falseAlarms, hit, iterations = 200, eqVar = TRUE){

  if (length(falseAlarms) != length(hit)) ('Vectors containing hit and false alarm frequencies do not have the same length')

  parameters            <- c()
  results               <- c()
  value                 <- c()

  # Function calculating total squared prediction error for hit and false alarm rates
  solver  <- function(x) {

    if (eqVar == TRUE) {
      rt <- exp(x[1]) / (1 + exp(x[1]))
      rl <- 0
      dpri <- exp(x[2])
      sd_target <- 1
      crit <- c()
      for (i in c(1:length(falseAlarms))) {
        crit[i] <- x[2 + i]
      }
    } else if (eqVar == FALSE) {
      rt <- exp(x[1]) / (1 + exp(x[1]))
      rl <- 0
      dpri <- exp(x[2])
      sd_target <- exp(x[3])
      crit <- c()
      for (i in c(1:length(falseAlarms))) {
        crit[i] <- x[3 + i]
      }
    }

    I <- c(-Inf,crit,Inf)

    # Lure items: never recollected
    predFA <- vector()
    for (i in 1:length(falseAlarms)) {
      predFA[i] <- log(stats::pnorm(I[i+1],mean=0,sd=1)-stats::pnorm(I[i],mean=0,sd=1) )
    }

    # Old items (1 -5) [not recollected]
    predHit <- vector()
    for (i in 1:(length(hit) - 1)) {
      predHit[i] <- log( (1-rt)*(stats::pnorm(I[i+1],mean=dpri,sd=sd_target)-stats::pnorm(I[i],mean=dpri,sd=sd_target)) )
    }

    # Old items (6) [Recollected or not recollected]
    predHit[6] <- log( rt + ((1-rt)*(1-stats::pnorm(I[6],mean=dpri,sd=sd_target))) )

    hitLik <- hit[i] * predHit[i]
    faLik <- falseAlarms[i] * predFA[i]

    return( -sum(c(hitLik,faLik)) )
  }

  # estimate starting parameters from data
  rstart <- makeMLEstartparameters(falseAlarms,hit,iterations = iterations)
  x0 <- NULL
  for (i in 1:iterations) {

    criteria <- rstart[,grepl( "Cs" , names(rstart))]
    if (eqVar == TRUE) {

      modelpara <- rstart[,c("rt","dpri")]
      x0 <- cbind(modelpara,criteria)
    } else if (eqVar == FALSE) {

      modelpara <- rstart[,c("rt","dpri","sd_target")]
      x0 <- cbind(modelpara,criteria)

    }


    cat('\rProgress: |',rep('=',floor((i/iterations)*50)),rep(' ',50 - floor((i/iterations)*50)),'|', sep = '')

    # Optimize
    control <- list('maxit', 10000000, 'reltol', 0.0000000001)
    temp    <- try(stats::optim(x0[i,], solver, method = "BFGS", control = control), silent = TRUE)

    # Move to next iteration if it crashes out of one
    if (class(temp) == "try-error") {
      parameters[i]        <- NA
      value[i]             <- NA
    } else {
      parameters           <- rbind(parameters, temp$par)
      value                <- rbind(value, temp$value)
    }

  }
  # Safe-guard against wonky runs
  parameters <- parameters[stats::complete.cases(parameters),]
  value <- value[stats::complete.cases(value),]
  #Identify run with min negLL
  Best <- parameters[which(value == min(value, na.rm = TRUE)),]

  #Prepare output
  Bestcolumns <- c("recollection_target","recollection_lure","familiarity","sd_target")
  fanames<-NULL
  for (i in c(1:length(falseAlarms))){
    fanames[i] <- paste0("c",i)
  }
  resultscolnames<-c(Bestcolumns,fanames,"negLL")

  tempresult <- NULL
  if (eqVar == TRUE) {
    tempresult <-
      c(exp(Best[1]) / (1 + exp(Best[1])),
        0,
        exp(Best[2]),
        1,
        Best[3:length(Best)],
        min(value))

  } else if (eqVar == FALSE) {
    tempresult <-
      c(exp(Best[1]) / (1 + exp(Best[1])),
        0,
        exp(Best[2]),
        exp(Best[3]),
        Best[4:length(Best)],
        min(value)
      )
  }

  results <- as.data.frame(matrix(tempresult,nrow=1, dimnames = list(NULL, resultscolnames)))
  cat('\n')
  cat('\n')
  return(results)
}
