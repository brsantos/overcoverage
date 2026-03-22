library(microbenchmark)

set.seed(123)

N <- 200
T <- 6
L <- 3

y <- matrix(sample(1:10, N*T, replace = TRUE), N, T)
covariates <- matrix(rnorm(N * 11), N, 11)
age2 <- matrix(rbinom(N*T,1,.5), N, T)
age3 <- matrix(rbinom(N*T,1,.5), N, T)
age4 <- matrix(0, N, T)
tin2 <- matrix(rbinom(N*T,1,.5), N, T)
tin3 <- matrix(rbinom(N*T,1,.5), N, T)
combins <- matrix(sample(1:(64 * 2^(L - 2)), N * T, replace = TRUE), N, T)
first <- rep(1, N)
init <- c(1,0,0,0,0,0,0,0)

X <- matrix(rnorm((64 * 2^(L - 2)) * 77), 64 * 2^(L - 2), 77)
proportions <- rep(1 / N, N)
n <- N

params <- rnorm(200)

microbenchmark(
  fast = {
    set.seed(123)
    loglikelihood_BLB_parallel(
      y, covariates, age2, age3, age4, tin2, tin3,
      X, init, first, L, combins, n, proportions, params
    )
  },
  legacy = {
    set.seed(123)
    loglikelihood_BLB_parallel_legacy(
      y, covariates, age2, age3, age4, tin2, tin3,
      X, init, first, L, combins, n, proportions, params
    )
  },
  times = 100
)

set.seed(123)
fast_val <- loglikelihood_BLB_parallel(
  y, covariates, age2, age3, age4, tin2, tin3,
  X, init, first, L, combins, n, proportions, params
)
set.seed(123)
legacy_val <- loglikelihood_BLB_parallel_legacy(
  y, covariates, age2, age3, age4, tin2, tin3,
  X, init, first, L, combins, n, proportions, params
)

stopifnot(isTRUE(all.equal(fast_val, legacy_val, tolerance = 1e-8)))
