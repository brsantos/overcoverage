library(microbenchmark)

set.seed(123)

N <- 200
T <- 6
L <- 3

y <- matrix(sample(1:10, N * T, replace = TRUE), N, T)
covariates <- matrix(rnorm(N * 11), N, 11)
age2 <- matrix(rbinom(N * T, 1, .5), N, T)
age3 <- matrix(rbinom(N * T, 1, .5), N, T)
age4 <- matrix(rbinom(N * T, 1, .5), N, T)
tin2 <- matrix(rbinom(N * T, 1, .5), N, T)
tin3 <- matrix(rbinom(N * T, 1, .5), N, T)
combins <- matrix(sample(1:(64 * 2^(L - 2)), N * T, replace = TRUE), N, T)
first <- rep(1, N)
init <- c(1, 0, 0, 0, 0, 0, 0, 0)

X <- matrix(rnorm((64 * 2^(L - 2)) * 77), 64 * 2^(L - 2), 77)

params <- rnorm(200)

microbenchmark(
  fast = loglikelihood_nonBLB_parallel(
    y, covariates, age2, age3, age4, tin2, tin3,
    X, init, first, L, combins, params
  ),
  legacy = loglikelihood_nonBLB_parallel_legacy(
    y, covariates, age2, age3, age4, tin2, tin3,
    X, init, first, L, combins, params
  ),
  times = 100
)

fast_val <- loglikelihood_nonBLB_parallel(
  y, covariates, age2, age3, age4, tin2, tin3,
  X, init, first, L, combins, params
)
legacy_val <- loglikelihood_nonBLB_parallel_legacy(
  y, covariates, age2, age3, age4, tin2, tin3,
  X, init, first, L, combins, params
)

stopifnot(isTRUE(all.equal(fast_val, legacy_val, tolerance = 1e-8)))
