# fullDPROC: Rpackage

* author: Nicholas Lange
* contact: lange.nk@gmail.com
* cran link: tba
* Source: [fullDPROC_0.1.0.tar.gz](fullDPROC_0.1.0.tar.gz)

This package fits the univariate dual-process model (Yonelinas, 1999) to recognition or source memory data by fitting the model to the ROC curve (minimizing SSE) or to the frequencies of responses by minimizing log-likelihood.

All fitting can be done by using the [DPSD](R/DPSD.R) wrapper function and specifying data and model/parameter restrictions in the function arguments.
The predicted ROC curve can be produced using the [predictedROC](R/predictedROC.R) function.

## [DPDS.R](R/DPSD.R) function

This is the main wrapper function. The package requires individual item-level confidence rating responses from recognition judgments or source judgments (for experiments with a single source dimension, i.e. a single rating scale ranging from Source A to Source B).

### Arguments

* responseScale: Vector containing  possible levels of confidence rating responses ordered from highest to lowest (e.g. 6:1).
* confidenceRatings: Vector containing participant confidence rating responses to individual items.
* TargetLure: Vector coding containing information whether each individual item is a Target or Lure item (e.g. Old, New, Old, Old). 
* targetLabel: A string/integer designating the label for the target items in TargetLure. If fitSource = TRUE, one source has to be designated target, the other lure.
* lureLabel: A string/integer designating the label for the lure items in TargetLure. If fitSource = TRUE, one source has to be designated target, the other lure.
* iterations (default = 200) A numeric value specifying the number of iterations.
* eqRecollection (default = FALSE) A boolean value specifying if recollection is set equal for the target and lure source (TRUE) or is estimates separately for both sources (FALSE).
* eqVar (default = TRUE) A boolean value specifying if the standard deviation of the target distribution is equal to that of the lure distribution (i.e. = 1) (TRUE) or estimated separately (FALSE).
* fitSource (default = TRUE) A boolean value specifying if the model is being fitted to source memory data (TRUE) or recognition memory data (FALSE)
* fitROC (default = FALSE) A boolean specifying if the model is fitted by minimizing SSE to ROC (hit/fa rates) or by MLE to rating frequencies

### Parameter restrictions

Recollection is bounded to be between 0 and 1, Familiarity to be positive. Criteria are ordered.

* eqVar = FALSE: sd_target is estimated, bounded to be positive.
* fitSource = FALSE: recollection_lure = 0 (overrides eqRecollection argument)

### Optimization

* fit with BFGS in optim
* Starting parameters:
	* for MLE: for 40% of iterations, starting parameters are estimated (mostly) from data and 60% are drawn randomly from uniform distributions
	* for SSE minimization: drawn randomly from uniform distributions with reasonable range

### Output

* Parameter estimates (recollection_target, recollection_lure, familiarity, sd_target, criteria) for the iteration with the lowest SSE or likelihood)
	* if eqVar = TRUE -> sd_target = 1
	* if eqRecollection = TRUE -> recollection_lure = recollection_target
	* if fitSource = FALSE -> recollection_lure = 0
* minimized value:
	* SSE if fitROC = TRUE
	* negative log-likelihood if fitROC = FALSE

## [predictedROC](R/predictedROC.R) function

This function provides the output for the predicted ROC curve from parameter values (for graphs and such).

### Arguments

* recollection_target
* recollection_lure
* familiarity
* sd_target (default = 1)
* fitSource (default = TRUE)

### Output

* dataframe with false alarms (fa) in first column for approximately continuous 0 - 1, and corresponding hit rates (hit) in second column

## Auxilliary functions

These are functions that are called from the wrapper function depending on arguments, but do not need to be called manually.

* [fitDPSD.R](R/fitDPSD.R): Function that fits recognition data (fitSource = FALSE) using MLE on rating frequencies
* [fitDPSDsource.R](R/fitDPSDsource.R): Function that fits source data (fitSource = TRUE) using MLE on rating frequencies
* [fitDPSDROC.R](R/fitDPSDROC.R): Function that fits recognition data (fitSource = FALSE) by minimizing SSE on cumulative rates
* [fitDPSDROCsource.R](R/fitDPSDROCsource.R): Function that fits source data (fitSource = TRUE) by minimizing SSE on cumulative rates
* [cumRates.R](R/cumRates.R): Function that transforms the raw rating data into cumulative hit/fa rates
* [ratingFreq.R](R/ratingFreq.R): Function that transforms the raw rating data into frequencies of ratings per rating bin and target/lure
* [makeMLEstartparameters.R](R/makeMLEstartparameters.R): Function that estimates starting parameters for MLE optimization. For 2/5 of iterations, starting parameters are chosen from normal distributions centred on parameter estimates taken from the observed data; for remaining iterations, starting parameters are drawn from uniform distributions.

## To Do:

* Integrate Bayesian estimation from [Hierarchical Bayes version](https://github.com/nklange/HierBayes/tree/master/Univariate)
* Tidy for cran
	* vignette
	* example data or simulation function for some data
* Also accept already summarized cumrates or frequency data rather than only item-level data
* Make predictedROC talk directly to output from DPSD