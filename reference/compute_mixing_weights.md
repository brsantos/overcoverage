# Compute mixing probabilities for a partition

Computes the mixing probabilities for each individual using the compiled
C++ routine `mixing_probs`. This is packaged as a wrapper function that
accepts inputs directly and optionally saves the results.

## Usage

``` r
compute_mixing_weights(
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

- estimates:

  Parameter estimates matrix (bootstraps x parameters).

- L:

  Number of registers.

- init:

  Initial state distribution.

- threads:

  Number of threads for RcppParallel.

- cpp_file:

  Optional path to a C++ file to source via Rcpp::sourceCpp.

- save_path:

  Optional path to save the mixing weights as an RDS file.

## Value

A 3D array of mixing weights (individuals x 2 x bootstraps).
