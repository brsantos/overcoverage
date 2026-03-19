# Prepare Register Variables for Modeling

Standardize and clean register indicators (e.g., convert to binary,
replace missing values, and zero-out registers on emigration/death).

## Usage

``` r
oc2_prepare_register_data(
  data,
  register_cols,
  emig_col = "emig",
  death_col = "death",
  binary_rules = NULL,
  na_to_zero = TRUE,
  zero_on_emig = TRUE,
  zero_on_death = TRUE
)
```

## Arguments

- data:

  A data.frame in long (person-year) format.

- register_cols:

  Character vector of register indicator columns.

- emig_col:

  Column name for emigration indicator (0/1).

- death_col:

  Column name for death indicator (0/1).

- binary_rules:

  Optional named list describing binary conversions. Each entry should
  include \`source\` and optionally \`threshold\`.

- na_to_zero:

  Logical; if TRUE, replace NA with 0 in register columns.

- zero_on_emig:

  Logical; if TRUE, set register columns to 0 on emigration.

- zero_on_death:

  Logical; if TRUE, set register columns to 0 on death.

## Value

A data.frame with cleaned register columns.
