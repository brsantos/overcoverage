# Prepare register variables for model input

Converts numeric registers to binary (optional), replaces missing
values, and zeroes registers in years with emigration or death.

## Usage

``` r
oc2_prepare_register_data(data, register_cols, emig_col = "emig",
  death_col = "death", binary_rules = NULL, na_to_zero = TRUE,
  zero_on_emig = TRUE, zero_on_death = TRUE)
```

## Arguments

- data:

  A data.frame in long format with one row per person-year.

- register_cols:

  Character vector of register column names.

- emig_col:

  Column name for emigration indicator.

- death_col:

  Column name for death indicator.

- binary_rules:

  Optional named list defining binary conversion rules.

- na_to_zero:

  Logical; replace NA in registers with 0.

- zero_on_emig:

  Logical; zero out registers when emig == 1.

- zero_on_death:

  Logical; zero out registers when death == 1.

## Value

A data.frame with updated register columns.
