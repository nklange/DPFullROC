#' Extraction of cumulative hit and false alarm rates for memory ROC analysis
#'
#' This function allows you to extract cumulative hit and false alarm rates for memory ROC analysis.
#' @param responseScale An vector containing  possible levels of confidence rating responses ordered from highest to lowest (e.g. 6:1).
#' @param confidenceRatings An vector containing participant confidence rating responses to individual items.
#' @param TargetLure An vector coding containing information whether each individual item is a Target or Lure item (e.g. Old, New, Old, Old).
#' @param targetLabel A string/integer designating the label for the target source in TargetLure.
#' @param lureLabel A string/integer designating the label for the lure source in TargetLure.
#' @return The function returns a data frame with components:
#' \item{fa}{The extracted frequency of hit responses per confidence rating.}
#' \item{hit}{The extracted frequency of fa responses per confidence rating.}
#' @export

ratingFreq <- function(responseScale, confidenceRatings, TargetLure, targetLabel, lureLabel){
  scaleLength <- length(responseScale)
  falseAlarm  <- c()
  hit         <- c()
  results     <- c()
  lureStimuli <- confidenceRatings[which(TargetLure == lureLabel)]
  targetStimuli <- confidenceRatings[which(TargetLure == targetLabel)]
  for(i in 1:(scaleLength)){
      hit[i]        <-
        sum(targetStimuli == rev(responseScale)[i])

      falseAlarm[i] <-
        sum(lureStimuli == rev(responseScale)[i])


    }
  results$fa         <- falseAlarm
  results$hit        <- hit

  return(data.frame(results))
}
