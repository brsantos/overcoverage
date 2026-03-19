\name{oc2_prepare_register_data}
\alias{oc2_prepare_register_data}
\title{Prepare register variables for model input}
\usage{
oc2_prepare_register_data(data, register_cols, emig_col = "emig",
  death_col = "death", binary_rules = NULL, na_to_zero = TRUE,
  zero_on_emig = TRUE, zero_on_death = TRUE)
}
\arguments{
  \item{data}{A data.frame in long format with one row per person-year.}
  \item{register_cols}{Character vector of register column names.}
  \item{emig_col}{Column name for emigration indicator.}
  \item{death_col}{Column name for death indicator.}
  \item{binary_rules}{Optional named list defining binary conversion rules.}
  \item{na_to_zero}{Logical; replace NA in registers with 0.}
  \item{zero_on_emig}{Logical; zero out registers when emig == 1.}
  \item{zero_on_death}{Logical; zero out registers when death == 1.}
}
\value{
  A data.frame with updated register columns.
}
\description{
  Converts numeric registers to binary (optional), replaces missing values,
  and zeroes registers in years with emigration or death.
}
