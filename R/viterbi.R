#' Compute Viterbi paths for each individual
#'
#' Runs the Viterbi algorithm using the compiled C++ routine \code{viterbi}.
#' Inputs are provided directly; any file reads should happen outside this
#' wrapper.
#'
#' @param y Observation matrix (individuals x time).
#' @param covariates Covariate matrix (individuals x 11).
#' @param age Age indicator array (individuals x 2 x time).
#' @param tin Time-in indicator array (individuals x 2 x time).
#' @param combins Register combination matrix.
#' @param estimates Parameter estimates matrix (bootstraps x parameters).
#' @param mixing_weights Mixing weights array (individuals x 2 x bootstraps).
#' @param L Number of registers.
#' @param init Initial state distribution.
#' @param threads Number of threads for RcppParallel.
#' @param cpp_file Optional path to a C++ file to source via Rcpp::sourceCpp.
#' @param save_path Optional path to save the optimal paths as an RDS file.
#'
#' @return A 3D array of optimal paths (individuals x states x bootstraps).
#' @export
compute_viterbi_paths <- function(
  y,
  covariates,
  age,
  tin,
  combins,
  estimates,
  mixing_weights,
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
  if (length(dim(age)) != 3 || dim(age)[2] != 2) {
    stop("age must be a 3D array with 2 levels (individuals x 2 x time).")
  }
  age2 <- age[, 1, ]
  age3 <- age[, 2, ]
  tin2 <- tin[, 1, ]
  tin3 <- tin[, 2, ]

  age2 <- apply(age2, 2, as.numeric)
  age3 <- apply(age3, 2, as.numeric)
  tin2 <- apply(tin2, 2, as.numeric)
  tin3 <- apply(tin3, 2, as.numeric)

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
  n_age <- 2
  dummies <- matrix(0, nrow = nrow(X), ncol = n_age)
  for (i in 1:n_age) {
    mat <- matrix(0, nrow = nrow(X), ncol = n_age)
    mat[, i] <- 1
    dummies <- rbind(dummies, mat)
  }
  X_extended <- X[rep(1:nrow(X), n_age + 1), ]
  X <- cbind(X_extended, dummies)
  # all two way interactions between lists
  for (p in 1:(L - 1)) {
    for (q in (p + 1):L) {
      X <- cbind(X, X[, p] * X[, q])
    }
  }
  # interactions between lists and covariates
  covariate_start <- L + 1
  covariate_end <- L + n_cob + n_age
  for (p in 1:L) {
    for (q in covariate_start:covariate_end) {
      X <- cbind(X, X[, p] * X[, q])
    }
  }

  # need to find the first detection for all individuals
  get.first <- function(x) min(which(x != 0))
  first <- apply(y, 1, get.first)
  first <- as.vector(first)

  RcppParallel::setThreadOptions(numThreads = threads)

  num_bootstraps <- nrow(estimates)
  optimal <- array(NA, dim = c(nrow(y), 17, num_bootstraps))

  for (boot in 1:num_bootstraps) {
    optimal[, , boot] <- viterbi(
      y, covariates, age2, age3, tin2, tin3,
      combins, X, init, first, L, estimates[boot, ],
      mixing_weights[, , boot]
    )
  }

  if (!is.null(save_path)) {
    saveRDS(optimal, file = save_path)
  }

  optimal
}
