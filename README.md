# overcoverage2

Capture-Recapture Hidden Markov Models for register-based population inference.
The package includes reusable checks to detect inconsistencies in register
data before modeling.

## Installation

Install from GitHub:

```r
install.packages("remotes")
remotes::install_github("brsantos/overcoverage2")
```

On Linux, you may need TBB for RcppParallel before installation, for example:

```
sudo apt-get install -y libtbb-dev
```

## Quick start: model_BLB

Below is a minimal example using simulated inputs. Replace these with your data and starting values.

```r
library(overcoverage2)

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
initial <- c(1, 0, 0, 0, 0, 0, 0, 0)

X <- matrix(rnorm((64 * 2^(L - 2)) * 77), 64 * 2^(L - 2), 77)

params <- rnorm(200)

fit <- model_BLB(
  y = y,
  covariates = covariates,
  age2 = age2,
  age3 = age3,
  age4 = age4,
  tin2 = tin2,
  tin3 = tin3,
  X = X,
  initial = initial,
  first = first,
  L = L,
  combins = combins,
  params = params,
  n = N,
  proportions = rep(1 / N, N)
)

fit
```
