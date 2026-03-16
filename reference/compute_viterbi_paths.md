# Compute Viterbi paths for each individual

Runs the Viterbi algorithm using the compiled C++ routine `viterbi`.
Inputs are provided directly; any file reads should happen outside this
wrapper.

## Usage

``` r
compute_viterbi_paths(y, covariates, age, tin, combins, estimates,
  mixing_weights, L = 8, init = c(1, 0, 0, 0, 0, 0, 0, 0),
  threads = 8, cpp_file = NULL, save_path = NULL)
```

## Arguments

- y:

  Observation matrix (individuals x time).

- covariates:

  Covariate matrix (individuals x 11).

- age:

  Age indicator array (individuals x 3 x time).

- tin:

  Time-in indicator array (individuals x 2 x time).

- combins:

  Register combination matrix.

- estimates:

  Parameter estimates matrix (bootstraps x parameters).

- mixing_weights:

  Mixing weights array (individuals x 2 x bootstraps).

- L:

  Number of registers.

- init:

  Initial state distribution.

- threads:

  Number of threads for RcppParallel.

- cpp_file:

  Optional path to a C++ file to source via Rcpp::sourceCpp.

- save_path:

  Optional path to save the optimal paths as an RDS file.

## Value

A 3D array of optimal paths (individuals x states x bootstraps).
