# Fit the non-BLB model

Fits the model without the Bag of Little Bootstraps using the compiled
C++ log-likelihood routine `loglikelihood_nonBLB_parallel`. Inputs are
provided directly; any file reads should happen outside this wrapper.

## Usage

``` r
model_nonBLB(
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
)
```

## Arguments

- y:

  Observation matrix (individuals x time).

- covariates:

  Covariate matrix (individuals x 11).

- age:

  Age indicator array (individuals x 2 x time). If a 3-level array is
  supplied, the third level is used; if a 2-level array is supplied, the
  third level is treated as all zeros.

- tin:

  Time-in indicator array (individuals x 2 x time).

- combins:

  Register combination matrix.

- L:

  Number of registers.

- init:

  Initial state distribution.

- sample_n:

  Optional sample size to subsample individuals.

- seed:

  Optional RNG seed used when sampling or generating init params.

- init_params:

  Initial parameter vector. If NULL, a random vector is used.

- num_betas:

  Number of beta parameters (default 78).

- num_deltas:

  Number of delta parameters (default 112).

- threads:

  Number of threads for RcppParallel.

- cpp_file:

  Optional path to a C++ file to source via Rcpp::sourceCpp.

- save_path:

  Optional path to save the parameter estimates as an RDS file.

## Value

The optimization result object from `nlm`.
