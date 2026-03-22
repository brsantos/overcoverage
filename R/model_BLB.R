#' Run the BLB model
#'
#' Runs the BLB log-likelihood optimization loop for the model described in the
#' paper. For confidentiality, any data reads are commented out; inputs must be
#' provided directly.
#'
#' @param part Partition index used for file naming when reading or saving.
#' @param y Observation matrix (individuals x time).
#' @param covariates Covariate matrix (individuals x 11).
#' @param age Age indicator array (individuals x 2 x time). If a 3-level array
#'   is supplied, the third level is used; if a 2-level array is supplied, the
#'   third level is treated as all zeros.
#' @param tin Time-in indicator array (individuals x 2 x time).
#' @param combins Register combination matrix.
#' @param init_params Initial parameter vector.
#' @param init Initial state distribution.
#' @param N Population size. Defaults to nrow(y).
#' @param L Number of registers.
#' @param num_bootstraps Number of bootstrap runs.
#' @param boot_start First bootstrap index to run.
#' @param threads Number of threads for RcppParallel.
#' @param progress_path Optional path for saving incremental results.
#' @param results_path Optional path for saving final estimates.
#' @param cpp_file Optional path to a C++ file to source via Rcpp::sourceCpp.
#'
#' @return A list with elements \code{estimates}, \code{means}, \code{sds},
#'   and \code{results}.
#' @export
#'
#' @examples
#' # See model_BLB_simulated_example() for a simulated run
model_BLB <- function(
  part = 20,
  y = NULL,
  covariates = NULL,
  age = NULL,
  tin = NULL,
  combins = NULL,
  init_params = NULL,
  init = c(1, 0, 0, 0, 0, 0, 0, 0),
  N = NULL,
  L = 8,
  num_bootstraps = 100,
  boot_start = 1,
  threads = 8,
  progress_path = NULL,
  results_path = NULL,
  cpp_file = NULL
) {

  # Example confidential loads (commented out):
  # y <- readRDS(file = paste0("/ssb/stamme02/stud590/wk24/lyb/partitions/y/y", part, ".rds"))
  # covariates <- readRDS(file = paste0("/ssb/stamme02/stud590/wk24/lyb/partitions/covariates/covariates", part, ".rds"))
  # age <- readRDS(file = paste0("/ssb/stamme02/stud590/wk24/lyb/partitions/age/age", part, ".rds"))
  # tin <- readRDS(file = paste0("/ssb/stamme02/stud590/wk24/lyb/partitions/tin/tin", part, ".rds"))
  # combins <- readRDS(file = paste0("/ssb/stamme02/stud590/wk24/lyb/partitions/combinations/combination", part, ".rds"))
  # init_params <- readRDS(file = "/ssb/stamme02/stud590/prog/lyb/init_params.rds")

  if (!is.null(cpp_file)) {
    Rcpp::sourceCpp(cpp_file)
  }

  missing_inputs <- c()
  if (is.null(y)) missing_inputs <- c(missing_inputs, "y")
  if (is.null(covariates)) missing_inputs <- c(missing_inputs, "covariates")
  if (is.null(age)) missing_inputs <- c(missing_inputs, "age")
  if (is.null(tin)) missing_inputs <- c(missing_inputs, "tin")
  if (is.null(combins)) missing_inputs <- c(missing_inputs, "combins")
  if (is.null(init_params)) missing_inputs <- c(missing_inputs, "init_params")
  if (length(missing_inputs) > 0) {
    stop("Missing required inputs: ", paste(missing_inputs, collapse = ", "))
  }

  if (is.null(N)) {
    N <- nrow(y)
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

  b <- dim(y)[1]
  probs <- rep(1 / b, b)

  # bootstrapping function done in sequence using lapply
  bootstrap <- function(y, covariates, age2, age3, age4, tin2, tin3, X, init,
                        first, L, combins, N, probs, init_params, boot) {

    cat("Running bootstrap", boot, "\n")

    if (!is.null(cpp_file)) {
      Rcpp::sourceCpp(cpp_file)
    }

    RcppParallel::setThreadOptions(numThreads = threads)

    objective <- function(params) {
      loglikelihood_BLB_parallel(y, covariates, age2, age3, age4, tin2, tin3,
                                 X, init, first, L, combins, N, probs, params)
    }

    result <- optim(par = init_params,
                    fn = objective,
                    method = "BFGS",
                    control = list(reltol = 1e-05))

    # return everything so we can ensure convergence
    return(result)
  }

  results <- vector("list", num_bootstraps)
  if (!is.null(progress_path) && file.exists(progress_path)) {
    results <- readRDS(progress_path)
  }

  for (boot in boot_start:num_bootstraps) {
    result <- bootstrap(y, covariates, age2, age3, age4, tin2, tin3, X,
                        init, first, L, combins, N, probs, init_params, boot)
    results[[boot]] <- result$par
    if (!is.null(progress_path)) {
      saveRDS(results, file = progress_path)
    }
  }

  estimates <- matrix(NA, nrow = num_bootstraps, ncol = length(init_params))
  for (boot in 1:num_bootstraps) {
    estimates[boot, ] <- results[[boot]]
  }

  means <- numeric(length(init_params))
  sds <- numeric(length(init_params))
  for (par in 1:length(init_params)) {
    means[par] <- mean(estimates[, par], na.rm = TRUE)
    sds[par] <- sd(estimates[, par], na.rm = TRUE)
  }

  if (!is.null(results_path)) {
    saveRDS(estimates, file = results_path)
  }

  list(
    estimates = estimates,
    means = means,
    sds = sds,
    results = results
  )
}
