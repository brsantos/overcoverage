# Check and clean register data for model preparation

Applies a set of consistency checks to register data and removes
individuals that violate the assumptions used in model preparation.

## Usage

``` r
oc2_check_register_data(data, id_col = "id", year_col = "year",
  firstimmig_col = "firstimmig", death_col = "death", emig_col = "emig",
  immig_col = "immig", reimmig_col = "reimmig", year_beginning,
  final_year)
```

## Arguments

- data:

  A data.frame in long format with one row per person-year.

- id_col:

  Column name for individual identifiers.

- year_col:

  Column name for calendar year.

- firstimmig_col:

  Column name for first immigration year.

- death_col:

  Column name for death indicator.

- emig_col:

  Column name for emigration indicator.

- immig_col:

  Column name for immigration indicator.

- reimmig_col:

  Column name for re-immigration indicator.

- year_beginning:

  First year of the study period.

- final_year:

  Last year of the study period.

## Value

A list with `data` (cleaned data.frame) and `removed` (IDs removed by
each check).
