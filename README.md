
<!-- README.md is generated from README.Rmd. Please edit that file -->

# overcoverage

![](inst/figures/overcoverage.png)

# overcoverage <img src="inst/figures/overcoverage.png" align="right" />

A package to estimate overcoverage on register based data considering
Multiple System Estimation (MSE) models, based on package `conting`.

## Prerequisites

Before using our package, one needs to install the archived package
`conting`. Because the package is archived, we need to install the
package in a different way. Using the package `devtools`, we can use the
following code

``` r
devtools::install_version("conting",
                          version = "1.7")
```

Then one can install our package, also using `devtools` with the
following lines

``` r
devtools::install_github("brsantos/overcoverage")
```

## Creating a population
