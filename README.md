
<!-- README.md is generated from README.Rmd. Please edit that file -->

# overcoverage <img src="inst/figures/overcoverage.png" align="right" />

A package to estimate overcoverage on register based data considering
Multiple System Estimation (MSE) models, based on the package `conting`.

The discussion of this method is available on the paper:

- Mussino, E., Santos, B., Monti, A. et al. Multiple systems estimation
  for studying over-coverage and its heterogeneity in population
  registers. **Quality & Quantity** (2023).
  <https://doi.org/10.1007/s11135-023-01757-x>

## Prerequisites

Before using this package, you need to install the archived package
`conting`. Because the package is archived, we need to install it in a
different way. Using the package `devtools`, we can use the following
code

``` r
devtools::install_version("conting",
                          version = "1.7")
```

After installing `conting`, you can install our package, also using
`devtools` with the following lines

``` r
devtools::install_github("brsantos/overcoverage")
```

## Creating a population
