#' Create a simulated population to test models for overcoverage.
#'
#' `create_population()` creates a simulated population with a defined number
#'   number of binary, continuous and categorical variables.
#'
#' @param size Size of population to be created.
#' @param n_bin_var Number of binary variables to be created in the population.
#' @param n_cont_var Number of continuous variables to be created in the
#'   population. Continuous variables are generated according to a Normal
#'   distribution.
#' @param n_cat_var Number of categorical variables to be created in the
#'   population. Categorical variables are always created with three categories:
#'   A, B and C.
#' @param prob_bin Vector with probabilities of each binary variable to be
#'   created. In case this does not have the same length of `n_bin_var` then
#'   the first value of the vector is repeated `n_bin_var` times.
#' @param prob_cat Vector with probabilities of categorical variables to be
#'   created. This must have length 3 and by default the values are 0.5, 0.3 and
#'   0.2.
#' @param names_bin_var Vector of names for binary variables.
#' @param names_cont Vector of names for continuous variables.
#' @param names_cat_var Vector of names for categorical variables.
#'
#' @return A `data.frame` with all the population with their respective
#'  variables.
#' @export
#'
#' @examples
#' # basic usage of create_population function
#' main_pop <- create_population(
#'   size = 1e6,
#'   n_cont_var = 1,
#'   n_cat_var = 3,
#'   prob_bin = c(0.5))

create_population <- function(size,
                              n_bin_var = 1,
                              n_cont_var = 1,
                              n_cat_var = 1,
                              prob_bin,
                              prob_cat = c(0.5, 0.3, 0.2),
                              names_bin_var = NULL,
                              names_cont = NULL,
                              names_cat_var = NULL){
  id <- replicate(size,
                  paste0(sample(c(letters, 0:9),
                                size = 8,
                                replace = TRUE), collapse = ""))

  if (length(prob_bin) != n_bin_var) prob_bin[1:n_bin_var] <- prob_bin[1]

  bin_var <- sapply(1:n_bin_var, function(a){
    sample(0:1, size,
           replace = TRUE,
           prob = c(1 - prob_bin[a], prob_bin[a]))
  })
  if (!is.null(names_bin_var)) colnames(bin_var) <- names_bin_var
  else colnames(bin_var) <- paste0("bin", 1:n_bin_var)

  cont_var <- sapply(1:n_cont_var, function(a){
    stats::rnorm(size)
  })
  if (!is.null(names_cont)) colnames(cont_var) <- names_cont
  else colnames(cont_var) <- paste0("cont", 1:n_cont_var)

  if (length(prob_cat) != 3) stop("prob_cat must have length 3.")

  cat_var <- sapply(1:n_cat_var, function(a){
    sample(c("A", "B", "C"), size,
           replace = TRUE,
           prob = prob_cat)
  })
  if (!is.null(names_cat_var)) colnames(cat_var) <- names_cat_var
  else colnames(cat_var) <- paste0("cat", 1:n_cat_var)

  data.frame(id, bin_var, cont_var, cat_var)
}
