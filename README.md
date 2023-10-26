
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

Because this type of data is sensitive and not publicly available, we
create functions that are able to recreate similar scenarios to showcase
the use of this package.

First, we can create a general database, with a similar purpose as the
Register of Total Population (RTB) in Sweden. We assume that if someone
enter the country, they are definitely in this dataset. For this use the
function `create_population`.

We can define how many people enter the country each year, for the
purpose of comparison with different countries, for instance. We can
create dummy variables, that could represent information such as sex
(male or female), employed (yes or no), among others. We can create
numerical variables that could represent income, age, time in the
country, etc. And lastly, we can create factor variables with three
levels. These could represent age groups, income groups or any other
type of grouping.

In case there is interest in another type of variable, adaptation for
this function are easy to implement. The binary variables are able to
receive the proportion of successful cases, for instance. We can easily
create a binary variable with equal proportion in the population, such
as sex, or we could create a variable with the proportion of people
employed, with values higher than 0.5 for example.

In the following example, we show how to create a population of 200,000
individuals. We create 2 binary variables, with names `sex` and
`higher_educated`. One categorical variable named `age_groups` and
numerical variable called `scaled_income`. By default, the categories
for the factor variable are sampled with probabilities 0.5, 0.3 and 0.2,
but this can be changed.

``` r
set.seed(42)

library(overcoverage)

main_pop <- create_population(
   size = 2e5,
   n_bin = 2, 
   n_cont_var = 1,
   n_cat_var = 1,
   prob_bin = c(0.5, 0.7), 
   names_bin = c("sex", "higher_educated"), 
   names_cont = "scaled_income", 
   names_cat_var = "age_groups")
```

We can make plots to look at their distribution across the population.
Because each variable is independent from the other, it seems that are
bars have the same height, but there are slightly variations in each one
of the bars.

``` r
library(ggplot2)
library(patchwork)

g <- ggplot(main_pop) + theme_minimal()
g1 <- g + geom_bar(aes(fill = factor(sex), 
                       x = factor(higher_educated)), 
                   position = "fill") + 
  labs(x = "higher_educated") + 
  scale_fill_viridis_d(name = "sex")
g2 <- g + geom_bar(aes(fill = factor(higher_educated), 
                       x = factor(sex)), 
                   position = "fill") +
  labs(x = "sex") + 
  scale_fill_viridis_d(name = "higher_educated")
g3 <- g + geom_bar(aes(fill = factor(age_groups), 
                       x = factor(higher_educated)), 
                   position = "fill") +
  labs(x = "higher_educated") + 
  scale_fill_viridis_d(name = "age_groups")
g4 <- g + geom_histogram(aes(x = scaled_income), 
                         fill = "royalblue", color = "grey75")

g1 + g2 + g3 + g4
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

## Presence in the country

Given a fixed population, we can create a matrix with the information
about when each one of these individuals were present in the country. We
can assume either that all of them arrived in the first year or that
they have arrived throughout a number of years. For the second option,
we can select whether they arrived constantly or not.

We control the probability of each individual leaving the country given
a linear predictor, which is used within a logistic regression model.
For instance, in the following example we can say that the probability
of leaving the country, $\phi$, is a function of age groups only with
the following equation

$$\log \begin{pmatrix}\frac{\phi_i}{1 - \phi_i} \end{pmatrix} = \beta_0 + \beta_1 X_B + \beta_2 X_C,$$

where $X_i$ is a indicator variable that is equal to 1, if `age_groups`
is equal to $i$ and 0 otherwise. In the following example, we set
$\beta_0 = 2$, $\beta_1 = -0.5$ and $\beta_2 = -1$ and we set the
arrivals to happen constantly over 4 years.

``` r
presences <- create_presences(main_pop,
   formula_phi = ~ age_groups,
   coef_values = c(2.5, -0.5, -1),
   years = 5, varying_arrival = TRUE)

# checking how many presences each year
colSums(presences)
#> [1]  50000  94522 134216 169668 151283
```

## Creating lists

We can create the
