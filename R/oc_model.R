#' Model to estimate over-coverage on population
#'
#' The `oc_model()` function is able to get an estimate of over-coverage in a
#'   population based on register data and log-linear models.
#'
#' @param model_formula Model formula to be used in the log-linear model.
#' @param freq_table Frequency table with all observational data.
#' @param censored Indexes of all individuals who are not observed in any of
#'  the registers.
#' @param nsample Number of posterior draws in the MCMC estimation process.
#' @param null.move.prob Parameter to control model selection algorithm. See
#'   `conting::bict()` for more information.
#' @param n.burnin Number of burnin samples to be discarded from the MCMC
#'   algorithm.
#' @param thin Thinning parameter in the MCMC algorithm.
#' @param prob_level Probability level to be used when calculating the credible
#'  interval for the overcoverage level.
#' @param ... Additional arguments to be passed to `conting::bict()` function,
#'  which is slightly changed here to function properly in the latest versions
#'  of R.
#'
#' @return A list with overcoverage estimates and summaries of the number of
#'  false positives based on estimates of the model. The object model itself is
#'  also returned.
#' @export
#'
#' @examples
#' \dontrun{model_oc <- oc_model(
#'  qty ~ bin1 * list1 + cat1 * list2,
#'  freq_table,
#'  cens_ind)}
#'
oc_model <- function(model_formula,
                     freq_table,
                     censored,
                     nsample = 2000,
                     null.move.prob = 1,
                     n.burnin = 1000,
                     thin = 1,
                     prob_level = 0.95, ...){

  model <- bict2(
    formula = model_formula,
    data = freq_table,
    cens = cens_ind,
    n.sample = 2000,
    null.move.prob = null.move.prob, ...)

  if (n.burnin > 0) {
    YY0 <- matrix(model$Y0[-(1:n.burnin), ],
                  ncol = dim(model$Y0)[2])
  } else YY0 <- model$Y0

  obs_z <- sum(model$maximal.mod$y[-model$missing2])

  fp_sample <- apply(YY0, 1, sum)
  fp_means <- apply(YY0, 2, mean)
  fp_lower <- apply(YY0, 2, quantile, (1 - prob_level)/2)
  fp_upper <- apply(YY0, 2, quantile, 1 - (1 - prob_level)/2)

  oc_estimates <- (sum(freq_table$count) -  (obs_z + fp_sample)) /
    sum(freq_table$count)

  list(oc_estimates = oc_estimates,
       fp_sample = fp_sample,
       fp_means = fp_means,
       fp_lower = fp_lower,
       fp_upper = fp_upper,
       model = model)
}

