#' Creation of list of people that deregister from the population.
#'
#' Given a population and a list of presences, `create_list_out()` creates a
#'   list of individuals who deregister.
#'
#' @param pop A baseline population from which one will select individuals to
#'   arrive and leave in a number of years.
#' @param presences A matrix of presences of all the individuals in the
#'   population for a number of years.
#' @param departures A vector of departure times of all the individuals in the
#'   population for a number of years.
#' @param formula_prob A linear predictor formula to calculate the probability
#'   of deregistering.
#' @param coef_values The values for each coefficient to be used in the linear
#'   predictor defined in `formula_prob`. The coefficients will be used
#'   considering a logit link function.
#'
#' @return A matrix with the information of the deregistering (1) or non
#'   non deregistering (0) for each individual of the population considering
#'   a number of years of observation.
#' @export

create_list_out <- function(pop,
                            presences,
                            departures,
                            formula_prob,
                            coef_values){

  n <- nrow(pop)

  pred_vars <- stats::model.matrix(formula_prob, data = pop)
  lin_predictor <-  as.numeric(pred_vars %*%
    matrix(coef_values, ncol = 1))

  prob_selection <- exp(lin_predictor) /
    (1 + exp(lin_predictor))

  aux_zeros <- rep(0, n)
  sapply(1:ncol(presences), function(a){
    departed <- which(departures == a)

    n_departed <- length(departed)
    if (length(n_departed) > 0){
      informed <- departed[stats::runif(n_departed) < prob_selection[departed]]
      aux_zeros[informed] <- 1
    }
    aux_zeros
  })
}
