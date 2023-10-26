#' Creating observations of lists in a population
#'
#' The function `create_list_presences` creates a matrix of presences based on
#'  a population of individuals and the matrix of presences.
#'
#' @param pop A baseline population from which one will select individuals to
#'   arrive and leave in a number of years.
#' @param presences A matrix of presences of all the individuals in the
#'   population for a number of years.
#' @param formula_prob A linear predictor formula to calculate the probability of
#'   appearing in such a list.
#' @param coef_values The values for each coefficient to be used in the linear
#'   predictor defined in `formula_prob`. The coefficients will be used
#'   considering a logit link function.
#'
#' @return A matrix with the information of the detection (1) or non detection
#'   (0) for each individual of the population considering a number of years of
#'   observation.
#' @export
#'
#' @examples
#' set.seed(111)
#'
#' # creating main population
#' main_pop <- create_population(
#'   size = 500,
#'   n_cat_var = 3,
#'   prob_bin = c(0.5))
#'
#' # creating matrix of presences.
#' presences <- create_presences(main_pop,
#'                               formula_phi = ~ bin1,
#'                               coef_values = c(2, -1),
#'                               varying_arrival = TRUE,
#'                               const_rate_arrival = FALSE,
#'                               years = 3)
#'
#' # creating list1 observation as a function of bin1
#' list1 <- create_list_presences(main_pop, presences,
#'                                ~ bin1,
#'                                c(1.5, -0.5))
#'
#' # checking the conditional distribution of list1 given bin1
#' prop.table(table(main_pop$bin1, list1[,1]), 1)
#'
#' # creating list1 observation as a function of cat1, that needs two
#' # coefficients, because cat1 has 3 categories.
#' list2 <- create_list_presences(main_pop, presences,
#'                                ~ cat1,
#'                                c(-0.5, -0.5, 1.0))
#'
#' # checking the conditional distribution of list2 given bin1
#' prop.table(table(main_pop$bin1, list2[, 1]), 1)
#'
#' # checking the conditional distribution of list2 given cat1
#' prop.table(table(main_pop$cat1, list2[, 1]), 1)


create_list_presences <- function(pop, presences, formula_prob, coef_values){

  n <- nrow(pop)

  pred_vars <- stats::model.matrix(formula_prob, data = pop)

  lin_predictor <- as.numeric(pred_vars %*% matrix(coef_values, ncol = 1))

  prob_selection <- exp(lin_predictor) / (1 + exp(lin_predictor))

  aux_zeros <- rep(0, n)
  apply(presences, 2, function(a){
    present <- which(a == 1)
    n_pres <- length(present)
    selected <- present[stats::runif(n_pres) < prob_selection[present]]
    aux_zeros[selected] <- 1
    aux_zeros
  })
}
