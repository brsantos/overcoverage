<img src="man/figures/overcoverage_small.png" align="right" width="140" />

# overcoverage

Capture-Recapture Hidden Markov Models for register-based population inference.
The package includes reusable checks to detect inconsistencies in register
data before modeling.

## Installation

Install from GitHub:

```r
install.packages("remotes")
remotes::install_github("brsantos/overcoverage")
```

On Linux, you may need TBB for RcppParallel before installation, for example:

```
sudo apt-get install -y libtbb-dev
```

## Quick start: model_BLB

Below is a minimal example using simulated inputs. Replace these with your data and starting values.

```r
library(overcoverage)

set.seed(123)
N <- 200
T <- 6
L <- 3

y <- matrix(sample(1:10, N * T, replace = TRUE), N, T)
covariates <- matrix(rnorm(N * 11), N, 11)
age <- array(rbinom(N * 2 * T, 1, .5), dim = c(N, 2, T))
tin <- array(rbinom(N * 2 * T, 1, .5), dim = c(N, 2, T))
combins <- matrix(sample(1:(64 * 2^(L - 2)), N * T, replace = TRUE), N, T)
init_params <- rnorm(184)

fit <- model_BLB(
  y = y,
  covariates = covariates,
  age = age,
  tin = tin,
  L = L,
  combins = combins,
  init_params = init_params,
  num_bootstraps = 2,
  threads = 1
)

fit
```
