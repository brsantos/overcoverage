model_nonBLB_simulated_example <- function(
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
  age <- array(rbinom(N * 3 * T, 1, 0.5), dim = c(N, 3, T))
  tin <- array(rbinom(N * 2 * T, 1, 0.5), dim = c(N, 2, T))
  combins <- matrix(sample(1:4, N * T, replace = TRUE), N, T)

  init_params <- rnorm(191)

  if (!run) {
    return(list(
      y = y,
      covariates = covariates,
      age = age,
      tin = tin,
      combins = combins,
      init_params = init_params,
      L = L
    ))
  }

  model_nonBLB(
    y = y,
    covariates = covariates,
    age = age,
    tin = tin,
    combins = combins,
    L = L,
    init_params = init_params,
    cpp_file = cpp_file
  )
}
