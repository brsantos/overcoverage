# overcoverage

## Register data preparation

The package provides reusable checks and preparation helpers for
register data. These functions help identify inconsistencies (for
example, emigration without re-immigration, or activity after death),
and standardize register variables before creating the model inputs.

### Consistency checks

``` r
library(overcoverage)

# data_long should be a person-year data.frame
# with columns like id, year, firstimmig, death, emig, immig, reimmig, ...
#
# checks <- oc2_check_register_data(
#   data_long,
#   year_beginning = 2002,
#   final_year = 2022
# )
# cleaned_data <- checks$data
# removed_ids <- checks$removed
```

### Register preparation

``` r
# register_cols <- c("married", "divorced", "amf", "studies", "intmove",
#                    "child", "pension", "job", "social", "faminc_b")
#
# prepared <- oc2_prepare_register_data(
#   cleaned_data,
#   register_cols = register_cols,
#   binary_rules = list(
#     pension = list(source = "aldpens", threshold = 0),
#     job = list(source = "forvers", threshold = 0),
#     social = list(source = "socink", threshold = 0)
#   )
# )
```

These helpers are designed to be adapted to other countries with
different register definitions. You can supply your own register lists,
covariates, and binary conversion rules while keeping the same checks.
