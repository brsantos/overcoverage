# Build Observation Inputs for model_BLB

Create observation matrices and combinations for the BLB model from
prepared register data. This helper is intended for internal use.

## Usage

``` r
oc2_prepare_model_BLB_inputs(
  data,
  id_col = "id",
  year_col = "year",
  register_cols,
  covariate_cols,
  reimmig_col = "reimmig",
  death_col = "death",
  emig_col = "emig",
  year_beginning,
  final_year,
  combins = NULL
)
```

## Arguments

- data:

  A data.frame in long (person-year) format.

- id_col:

  Column name for individual identifiers.

- year_col:

  Column name for calendar year.

- register_cols:

  Character vector of register indicator columns.

- covariate_cols:

  Character vector of covariate columns to include in observation
  combinations.

- reimmig_col:

  Column name for re-immigration indicator (0/1).

- death_col:

  Column name for death indicator (0/1).

- emig_col:

  Column name for emigration indicator (0/1).

- year_beginning:

  First year in the study window.

- final_year:

  Last year in the study window.

- combins:

  Optional matrix of combination indices used for unobserved years.

## Value

A list with elements \`y_matrix\`, \`X\`, and \`num_combos\`.
