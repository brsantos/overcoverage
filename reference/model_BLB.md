# Run the BLB model

Runs the BLB log-likelihood optimization loop for the model described in
the paper. For confidentiality, any data reads are commented out; inputs
must be provided directly.

## Usage

``` r
model_BLB(
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
)
```

## Arguments

- part:

  Partition index used for file naming when reading or saving.

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

- init_params:

  Initial parameter vector.

- init:

  Initial state distribution.

- N:

  Population size. Defaults to nrow(y).

- L:

  Number of registers.

- num_bootstraps:

  Number of bootstrap runs.

- boot_start:

  First bootstrap index to run.

- threads:

  Number of threads for RcppParallel.

- progress_path:

  Optional path for saving incremental results.

- results_path:

  Optional path for saving final estimates.

- cpp_file:

  Optional path to a C++ file to source via Rcpp::sourceCpp.

## Value

A list with elements `estimates`, `means`, `sds`, and `results`.

## Examples

``` r
# See model_BLB_simulated_example() for a simulated run
```
