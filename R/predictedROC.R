#' Returning false alarm and hit rates for the Dual Process Signal Detection (DPSD) for recognition and source memory
#'
#' This function allows you to get the correspending false alarm and hit rates for a given set of parameters (recollection
#' for both distributions, familiarity, standard deviation for the target distribution). A boolean
#' allows you to fit a source ROC when set TRUE, with target and lure distributions representing designated target source
#' and lure source, and a recognition ROC when set FALSE (with recollection_lure and familiarity set to 0 for the lure
#' distribution).
#'
#' @author Nicholas Lange, \email{lange.nk@gmail.com}
#' @param recollection_target A value representing the recollection rate for target items.
#' @param recollection_lure A value representing the recollection rate for lure items.
#' @param familiarity A value representing the familiarity rate
#' @param sd_target A value representing the standard deviation of the target distribution (default: sd_target = 1)
#' @param fitSource A boolean value representing whether a recognition or source memory ROC is fitted (default: fitSource = TRUE)
#' @return The function returns a data frame with components:
#' \item{fa}{The fitted cumulative false alarm rate.}
#' \item{hit}{The fitted cumulative hit rate.}
#' @references Yonelinas, A. P. (1999). The Contribution of Recollection and Familiarity to Recognition and Source-Memory Judgments: A Formal Dual-Process Model and an Analysis of Receiver Operating Characteristics. Journal of Experimental Psychology: Learning, Memory, and Cognition, 25(6), 1415 - 1434. http://doi.org/10.1037//0278-7393.25.6.1415
#' @keywords memory ROC
#' @export

predictedROC       <- function(recollection_target,recollection_lure, familiarity,
                                  sd_target = 1, fitSource = TRUE){

  rt <- recollection_target
  rl <- recollection_lure
  sdt <- sd_target
  sdl <- 1

  if (fitSource == FALSE) {

    rl <- 0
    mut <- familiarity
    mul <- 0

  }
  else {

    mut <- familiarity/2
    mul <- - familiarity/2
  }

  results             <- c()
  results$fa  <- c(seq(0.001, 0.99999, 0.01), 0.99999)
  criterion           <- stats::qnorm(results$fa/(1 - rl), - mul, sdl)
  criterion[which(results$fa/(1 - rl) > 1)] <- 1000
  results$hit        <- results$fa + rt + (1 - rt) * stats::pnorm(criterion, - mut, sdt) - (1 - rl) * stats::pnorm(criterion, - mul, sdl)
  results$hit[results$hit > 1] <- 1
  return(data.frame(results))
}
