# Check Register Data Consistency

Apply a set of consistency checks to person-year register data and
return the cleaned data alongside the IDs removed by each check.

## Usage

``` r
oc2_check_register_data(
  data,
  id_col = "id",
  year_col = "year",
  firstimmig_col = "firstimmig",
  death_col = "death",
  emig_col = "emig",
  immig_col = "immig",
  reimmig_col = "reimmig",
  year_beginning,
  final_year
)
```

## Arguments

- data:

  A data.frame in long (person-year) format.

- id_col:

  Column name for individual identifiers.

- year_col:

  Column name for calendar year.

- firstimmig_col:

  Column name for first immigration year.

- death_col:

  Column name for death indicator (0/1).

- emig_col:

  Column name for emigration indicator (0/1).

- immig_col:

  Column name for immigration indicator (0/1).

- reimmig_col:

  Column name for re-immigration indicator (0/1).

- year_beginning:

  First year in the study window.

- final_year:

  Last year in the study window.

## Value

A list with components:

- data:

  Filtered data.frame after removing inconsistent IDs.

- removed:

  Named list of ID vectors removed at each check.
