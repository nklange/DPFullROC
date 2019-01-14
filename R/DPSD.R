#' Wrapper for model fitting by calculation of cumulative hit and false alarm rates
#'
#' This function allows to estimate recollection and familiarity for recognition and source memory data by fitting
#' data to the DPSD model.
#' The optimization is attempted by minimizing the total squared difference between observed and
#' predicted cumulative hit and false alarm rates using the Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm
#' in \code{\link{optim}}. The function uses random start values on each iteration in order to find the set of parameters,
#' which fit the data best by returning the values with the lowest total squared difference.
#' Optional arguments in the function allow the user to specify if they are fitting recognition or source memory data,
#' an equal-variance model and if recollection is to be estimated as a separate parameter for both target and lure items,
#' when fitting source memory data.
#' Recollection is bounded to be between 0 and 1, Familiarity and the standard deviation of the target distribution to be positive.
#' Criteria are unbounded.
#' @author Nicholas Lange, \email{lange.nk@gmail.com}
#' @param responseScale An vector containing  possible levels of confidence rating responses ordered from highest to lowest (e.g. 6:1).
#' @param confidenceRatings An vector containing participant confidence rating responses to individual items.
#' @param TargetLure An vector coding containing information whether each individual item is a Target or Lure item (e.g. Old, New, Old, Old)
#' @param targetLabel A string/integer designating the label for the target source in TargetLure
#' @param lureLabel A string/integer designating the label for the lure source in TargetLure
#' @param iterations A numeric value specifying the number of iterations. Default is set to 200.
#' @param eqRecollection A boolean value specifying if recollection is set equal for the target and lure source (TRUE) or is estimates separately for both sources (FALSE). Default is set to FALSE.
#' @param eqVar A boolean value specifying if the standard deviation of the target distribution is equal to that of the lure distribution (i.e. = 1) (TRUE) or estimated separately (FALSE). Default is set to TRUE.
#' @param fitSource A boolean value specifying if the model is being fitted to source memory data (TRUE) or recognition memory data (FALSE). Default is set to TRUE
#' @param fitROC A boolean specifying if the model being fitted is fitted SSE to ROCs or MLE to raw category ratings
#' @return The function returns a dataframe with components:
#' \item{(parameters)}{The estimated parameters (recollection_target, recollection_lure, familiarity, sd_target, criteria) for the iteration with the lowest SumSquareError}
#' \item{SSE}{Minimum sum square error}
#' @references Yonelinas, A. P. (1999). The Contribution of Recollection and Familiarity to Recognition and Source-Memory Judgments: A Formal Dual-Process Model and an Analysis of Receiver Operating Characteristics. Journal of Experimental Psychology: Learning, Memory, and Cognition, 25(6), 1415 - 1434. http://doi.org/10.1037//0278-7393.25.6.1415
#' @keywords ROC recollection familiarity DPSD
#' @export
#'
DPSD <- function(responseScale, confidenceRatings, TargetLure, targetLabel, lureLabel,
                 iterations = 200, eqVar = TRUE, eqRecollection = FALSE, fitSource = TRUE, fitROC = FALSE){

  rates   <- cumRates(responseScale, confidenceRatings, TargetLure, targetLabel, lureLabel)
  frequencies   <- ratingFreq(responseScale, confidenceRatings, TargetLure, targetLabel, lureLabel)

  if (fitSource == TRUE & fitROC == TRUE) {
    results <- fitDPSDROCsource(rates$fa, rates$hit, iterations = iterations, eqVar = eqVar, eqRecollection = eqRecollection)
  }
  else if (fitSource == FALSE & fitROC == TRUE){
    results <- fitDPSDROC(rates$fa, rates$hit, iterations = iterations, eqVar = eqVar)
  }
  else if (fitSource == FALSE & fitROC == FALSE){
    results <- fitDPSD(frequencies$fa,frequencies$hit, iterations = iterations, eqVar = eqVar)
  }
  else if (fitSource == TRUE & fitROC == FALSE){
    results <- fitDPSDsource(frequencies$fa,frequencies$hit, iterations = iterations, eqVar = eqVar, eqRecollection = eqRecollection)
  }


  return(results)
}
