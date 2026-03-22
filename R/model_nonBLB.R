#' Fit the non-BLB model
#'
#' Fits the model without the Bag of Little Bootstraps using the compiled
#' C++ log-likelihood routine \code{loglikelihood_nonBLB_parallel}. Inputs are
#' provided directly; any file reads should happen outside this wrapper.
#'
#' @param y Observation matrix (individuals x time).
#' @param covariates Covariate matrix (individuals x 11).
#' @param age Age indicator array (individuals x 2 x time). If a 3-level array
#'   is supplied, the third level is used; if a 2-level array is supplied, the
#'   third level is treated as all zeros.
#' @param tin Time-in indicator array (individuals x 2 x time).
#' @param combins Register combination matrix.
#' @param L Number of registers.
#' @param init Initial state distribution.
#' @param sample_n Optional sample size to subsample individuals.
#' @param seed Optional RNG seed used when sampling or generating init params.
#' @param init_params Initial parameter vector. If NULL, a random vector is used.
#' @param num_betas Number of beta parameters (default 78).
#' @param num_deltas Number of delta parameters (default 112).
#' @param threads Number of threads for RcppParallel.
#' @param cpp_file Optional path to a C++ file to source via Rcpp::sourceCpp.
#' @param save_path Optional path to save the parameter estimates as an RDS file.
#'
#' @return The optimization result object from \code{nlm}.
#' @export
model_nonBLB <- function(
  y,
  covariates,
  age,
  tin,
  combins,
  L = 8,
  init = c(1, 0, 0, 0, 0, 0, 0, 0),
  sample_n = NULL,
  seed = NULL,
  init_params = NULL,
  num_betas = 78,
  num_deltas = 112,
  threads = 8,
  cpp_file = NULL,
  save_path = NULL
) {

  if (!is.null(cpp_file)) {
    Rcpp::sourceCpp(cpp_file)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.null(sample_n)) {
    samp <- sample(seq_len(nrow(y)), sample_n, replace = FALSE)
    y <- y[samp, ]
    covariates <- covariates[samp, ]
    age <- age[samp, , ]
    tin <- tin[samp, , ]
    combins <- combins[samp, ]
  }

  covariates <- as.matrix(covariates)
  if (length(dim(age)) != 3) {
    stop("age must be a 3D array (individuals x age levels x time).")
  }
  if (dim(age)[2] == 2) {
    age2 <- age[, 1, ]
    age3 <- age[, 2, ]
    age4 <- matrix(0, nrow = dim(age)[1], ncol = dim(age)[3])
  } else if (dim(age)[2] == 3) {
    age2 <- age[, 1, ]
    age3 <- age[, 2, ]
    age4 <- age[, 3, ]
  } else {
    stop("age must have 2 or 3 levels.")
  }
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

  if (is.null(init_params)) {
    logit <- function(x) log(x / (1 - x))
    init_params <- c(
      logit(runif(1, 0, 1)),
      rnorm(num_betas, 0, 1),
      rnorm(num_deltas, 0, 1)
    )
  }

  objective <- function(params) {
    loglikelihood_nonBLB_parallel(y, covariates, age2, age3, age4, tin2,
                                  tin3, X, init, first, L, combins, params)
  }

  RcppParallel::setThreadOptions(numThreads = threads)

  result <- nlm(f = objective,
                p = init_params,
                gradtol = 1e-06,
                iterlim = 200)

  if (!is.null(save_path)) {
    saveRDS(result$estimate, file = save_path)
  }

  result
}
