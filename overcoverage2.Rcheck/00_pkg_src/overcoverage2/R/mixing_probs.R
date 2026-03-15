#' Compute mixing probabilities for a partition
#'
#' Computes the mixing probabilities for each individual using the compiled
#' C++ routine \code{mixing_probs}. This is packaged as a wrapper function that
#' accepts inputs directly and optionally saves the results.
#'
#' @param y Observation matrix (individuals x time).
#' @param covariates Covariate matrix (individuals x 11).
#' @param age Age indicator array (individuals x 3 x time).
#' @param tin Time-in indicator array (individuals x 2 x time).
#' @param combins Register combination matrix.
#' @param estimates Parameter estimates matrix (bootstraps x parameters).
#' @param L Number of registers.
#' @param init Initial state distribution.
#' @param threads Number of threads for RcppParallel.
#' @param cpp_file Optional path to a C++ file to source via Rcpp::sourceCpp.
#' @param save_path Optional path to save the mixing weights as an RDS file.
#'
#' @return A 3D array of mixing weights (individuals x 2 x bootstraps).
#' @export
compute_mixing_weights <- function(
  y,
  covariates,
  age,
  tin,
  combins,
  estimates,
  L = 8,
  init = c(1, 0, 0, 0, 0, 0, 0, 0),
  threads = 8,
  cpp_file = NULL,
  save_path = NULL
) {

  if (!is.null(cpp_file)) {
    Rcpp::sourceCpp(cpp_file)
  }

  covariates <- as.matrix(covariates)
  age2 <- age[, 1, ]
  age3 <- age[, 2, ]
  age4 <- age[, 3, ]
  tin2 <- tin[, 1, ]
  tin3 <- tin[, 2, ]

  age2 <- apply(age2, 2, as.numeric)
  age3 <- apply(age3, 2, as.numeric)
  age4 <- apply(age4, 2, as.numeric)
  tin2 <- apply(tin2, 2, as.numeric)
  tin3 <- apply(tin3, 2, as.numeric)

  # need to find the first detection for all individuals
  get.first <- function(x) min(which(x != 0))
  first <- apply(y, 1, get.first)
  first <- as.vector(first)

  # creation of design matrix for observations (6 registers + kjoenn + citizenship)
  l <- rep(list(1:0), L)
  combinations <- expand.grid(l)
  X <- matrix(unlist(combinations), byrow = FALSE, ncol = L)
  # X contains possible observation combinations
  # incorporating country of birth
  n_cob <- 3
  dummies <- matrix(0, nrow = nrow(X), ncol = n_cob)
  for (i in 1:n_cob) {
    mat <- matrix(0, nrow = nrow(X), ncol = n_cob)
    mat[, i] <- 1
    dummies <- rbind(dummies, mat)
  }
  X_extended <- X[rep(1:nrow(X), n_cob + 1), ]
  X <- cbind(X_extended, dummies)
  # incorporating age
  n_age <- 3
  dummies <- matrix(0, nrow = nrow(X), ncol = n_age)
  for (i in 1:n_age) {
    mat <- matrix(0, nrow = nrow(X), ncol = n_age)
    mat[, i] <- 1
    dummies <- rbind(dummies, mat)
  }
  X_extended <- X[rep(1:nrow(X), n_age + 1), ]
  X <- cbind(X_extended, dummies)

  # all two way interactions between lists
  for (p in 1:5) {
    for (q in (p + 1):6) {
      X <- cbind(X, X[, p] * X[, q])
    }
  }
  # interactions between lists and covariates
  for (p in 1:6) {
    for (q in 7:14) {
      X <- cbind(X, X[, p] * X[, q])
    }
  }

  RcppParallel::setThreadOptions(numThreads = threads)

  num_bootstraps <- nrow(estimates)
  weights <- vector("list", num_bootstraps)
  for (boot in 1:num_bootstraps) {
    weights[[boot]] <- mixing_probs(
      y, covariates, age2, age3, age4, tin2, tin3,
      X, init, first, L, combins, estimates[boot, ]
    )
  }

  mixing_weights <- array(NA, dim = c(nrow(y), 2, num_bootstraps))
  for (boot in 1:num_bootstraps) {
    mixing_weights[, , boot] <- weights[[boot]]
  }

  if (!is.null(save_path)) {
    saveRDS(mixing_weights, file = save_path)
  }

  mixing_weights
}
