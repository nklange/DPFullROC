#' Extraction of cumulative hit and false alarm rates for memory ROC analysis
#'
#' This function allows you to extract cumulative hit and false alarm rates for memory ROC analysis.
#' @param responseScale An vector containing  possible levels of confidence rating responses ordered from highest to lowest (e.g. 6:1).
#' @param confidenceRatings An vector containing participant confidence rating responses to individual items.
#' @param TargetLure An vector coding containing information whether each individual item is a Target or Lure item (e.g. Old, New, Old, Old).
#' @param targetLabel A string/integer designating the label for the target source in TargetLure.
#' @param lureLabel A string/integer designating the label for the lure source in TargetLure.
#' @return The function returns a data frame with components:
#' \item{fa}{The extracted cumulative false alarm rate.}
#' \item{hit}{The extracted cumulative hit rate.}
#' \item{zfa}{The z-transformed cumulative false alarm rate.}
#' \item{zhit}{The z-transformed cumulative hit rate.}
#' @keywords ROC
#' @export

cumRates <- function(responseScale, confidenceRatings, TargetLure, targetLabel, lureLabel){
  scaleLength <- length(responseScale)
  falseAlarm  <- c()
  hit         <- c()
  zhit        <- c()
  zfalseAlarm <- c()
  results     <- c()
  lureStimuli <- confidenceRatings[which(TargetLure == lureLabel)]
  targetStimuli <- confidenceRatings[which(TargetLure == targetLabel)]
  for(i in 1:(scaleLength - 1)){
    if(i < 2){
      hit[i]        <-
        sum(targetStimuli == responseScale[i]) / length(targetStimuli)

      falseAlarm[i] <-
        sum(lureStimuli == responseScale[i]) / length(lureStimuli)

      zhit[i]       <- if (hit[i] == 1) {
        stats::qnorm(1 - 1 / (2*length(targetStimuli)))
      } else if (hit[i] == 0) {
        stats::qnorm(1 / (2*length(targetStimuli)))
      } else {
        stats::qnorm(hit[i])
      }

      zfalseAlarm[i] <- if (falseAlarm[i] == 1) {
        stats::qnorm(1 - 1 / (2*length(targetStimuli)))
      } else if (falseAlarm[i] == 0) {
        stats::qnorm(1 / (2*length(targetStimuli)))
      } else {
        stats::qnorm(falseAlarm[i])
      }

    } else {
      hit[i]        <- hit[i - 1] + sum(targetStimuli == responseScale[i]) / length(targetStimuli)
      falseAlarm[i] <- falseAlarm[i - 1] + sum(lureStimuli == responseScale[i]) / length(lureStimuli)
      zhit[i]       <- if (hit[i] == 1) {
        stats::qnorm(1 - 1 / (2*length(targetStimuli)))
      } else if (hit[i] == 0) {
        stats::qnorm(1 / (2*length(targetStimuli)))
      } else {
        stats::qnorm(hit[i])
      }

      zfalseAlarm[i] <- if (falseAlarm[i] == 1) {
        stats::qnorm(1 - 1 / (2*length(targetStimuli)))
      } else if (falseAlarm[i] == 0) {
        stats::qnorm(1 / (2*length(targetStimuli)))
      } else {
        stats::qnorm(falseAlarm[i])
      }
    }
  }
  results$fa         <- falseAlarm
  results$hit        <- hit
  results$zfa        <- zfalseAlarm
  results$zhit       <- zhit

  return(data.frame(results))
}
