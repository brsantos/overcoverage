\name{oc2_check_register_data}
\alias{oc2_check_register_data}
\title{Check and clean register data for model preparation}
\usage{
oc2_check_register_data(data, id_col = "id", year_col = "year",
  firstimmig_col = "firstimmig", death_col = "death", emig_col = "emig",
  immig_col = "immig", reimmig_col = "reimmig", year_beginning,
  final_year)
}
\arguments{
  \item{data}{A data.frame in long format with one row per person-year.}
  \item{id_col}{Column name for individual identifiers.}
  \item{year_col}{Column name for calendar year.}
  \item{firstimmig_col}{Column name for first immigration year.}
  \item{death_col}{Column name for death indicator.}
  \item{emig_col}{Column name for emigration indicator.}
  \item{immig_col}{Column name for immigration indicator.}
  \item{reimmig_col}{Column name for re-immigration indicator.}
  \item{year_beginning}{First year of the study period.}
  \item{final_year}{Last year of the study period.}
}
\value{
  A list with \code{data} (cleaned data.frame) and \code{removed}
  (IDs removed by each check).
}
\description{
  Applies a set of consistency checks to register data and removes
  individuals that violate the assumptions used in model preparation.
}
