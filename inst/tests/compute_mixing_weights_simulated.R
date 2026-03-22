compute_mixing_weights_simulated_example <- function(
  N = 200,
  T = 6,
  L = 8,
  seed = 123,
  run = FALSE,
  cpp_file = NULL
) {
  set.seed(seed)

  y <- matrix(sample(1:10, N * T, replace = TRUE), N, T)
  covariates <- matrix(rnorm(N * 11), N, 11)
  age <- array(rbinom(N * 2 * T, 1, 0.5), dim = c(N, 2, T))
  tin <- array(rbinom(N * 2 * T, 1, 0.5), dim = c(N, 2, T))
  combins <- matrix(sample(1:4, N * T, replace = TRUE), N, T)

  # Simulated parameter estimates (bootstraps x params)
  estimates <- matrix(rnorm(2 * 184), nrow = 2, ncol = 184)

  if (!run) {
    return(list(
      y = y,
      covariates = covariates,
      age = age,
      tin = tin,
      combins = combins,
      estimates = estimates,
      L = L
    ))
  }

  compute_mixing_weights(
    y = y,
    covariates = covariates,
    age = age,
    tin = tin,
    combins = combins,
    estimates = estimates,
    L = L,
    cpp_file = cpp_file
  )
}
