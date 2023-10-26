#' Create matrix of presences and absences.
#'
#' Based on a given population, `create_presences()` creates a matrix of
#'   presences and absences based on a probability function to remain in the
#'   country.
#'
#' @param pop A baseline population from which one will select individuals to
#'   arrive and leave in a number of years.
#' @param varying_arrival A logical parameter that controls if all individuals
#'   arrive in the first year or not. Default is FALSE.
#' @param const_rate_arrival In case individuals are arriving in different
#'   years, then if `TRUE` individuals will arrive at a constant rate each year.
#'   If `FALSE` then will arrive randomly in the different years.
#' @param formula_phi A linear predictor formula to calculate the probability of
#'   leaving the country.
#' @param coef_values The values for each coefficient to be used in the linear
#'   predictor defined in `formula_phi`. The coefficients will be used
#'   considering a logit link function.
#' @param years Number of years to define this matrix of absences and presences.
#'
#' @return A matrix with the information of the presence (1) or absence (0) for
#'   each individual of the population considering a number of years of
#'   observation.
#' @export
#'
#' @examples
#' main_pop <- create_population(
#' size = 500,
#' n_cat_var = 3,
#' prob_bin = c(0.5))
#'
#' # example with all individuals arriving in the first year.
#' presences <- create_presences(main_pop,
#'   formula_phi = ~ bin1,
#'   coef_values = c(2, -1),
#'   years = 3)
#' colSums(presences)
#'
#' # example with constant rate of arrival of individuals in each year.
#' presences <- create_presences(main_pop,
#'   formula_phi = ~ bin1,
#'   coef_values = c(2, -1),
#'   varying_arrival = TRUE,
#'   years = 3)
#' colSums(presences)
#'
#' # example with varying rate of arrival of individuals in each year.
#' presences <- create_presences(main_pop,
#'   formula_phi = ~ bin1,
#'   coef_values = c(2, -1),
#'   varying_arrival = TRUE,
#'   const_rate_arrival = FALSE,
#'   years = 3)
#' colSums(presences)

create_presences <- function(pop,
                             varying_arrival = FALSE,
                             const_rate_arrival = TRUE,
                             formula_phi,
                             coef_values,
                             years = 2){
  n <- nrow(pop)

  pred_vars <- stats::model.matrix(formula_phi, data = pop)

  lin_predictor <- as.numeric(pred_vars %*% matrix(coef_values, ncol = 1))

  prob_selection <- exp(lin_predictor) /
    (1 + exp(lin_predictor))

  if (!varying_arrival){
    z <- matrix(0, nrow = n, ncol = years)
    z[, 1] <- 1
  }
  else {
    if (const_rate_arrival){
      size_arrival_year <- n / (years - 1)
      z <- kronecker(diag(years - 1),
                     rep(1, size_arrival_year))
      z <- cbind(z, 0)
    }
    else {
      arrival_year <- sample(1:(years - 1),
                             size = n,
                             replace = TRUE)
      z <- matrix(0, nrow = n, ncol = years)
      for (k in 1:n) z[k, arrival_year[k]] <- 1
    }
  }
  for (j in 2:years){
    present <- which(z[, j - 1] == 1)
    z[present, j] <- as.numeric(
      stats::runif(length(present)) < prob_selection[present])
  }
  z
}
