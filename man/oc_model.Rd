% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oc_model.R
\name{oc_model}
\alias{oc_model}
\title{Model to estimate over-coverage on population}
\usage{
oc_model(
  model_formula,
  freq_table,
  censored,
  nsample = 2000,
  null.move.prob = 1,
  n.burnin = 1000,
  thin = 1,
  prob_level = 0.95,
  ...
)
}
\arguments{
\item{model_formula}{Model formula to be used in the log-linear model.}

\item{freq_table}{Frequency table with all observational data.}

\item{censored}{Indexes of all individuals who are not observed in any of
the registers.}

\item{nsample}{Number of posterior draws in the MCMC estimation process.}

\item{null.move.prob}{Parameter to control model selection algorithm. See
`conting::bict()` for more information.}

\item{n.burnin}{Number of burnin samples to be discarded from the MCMC
algorithm.}

\item{thin}{Thinning parameter in the MCMC algorithm.}

\item{prob_level}{Probability level to be used when calculating the credible
interval for the overcoverage level.}

\item{...}{Additional arguments to be passed to `conting::bict()` function,
which is slightly changed here to function properly in the latest versions
of R.}
}
\value{
A list with overcoverage estimates and summaries of the number of
 false positives based on estimates of the model. The object model itself is
 also returned.
}
\description{
The `oc_model()` function is able to get an estimate of over-coverage in a
  population based on register data and log-linear models.
}
\examples{
\dontrun{
model_oc <- oc_model(
 qty ~ bin1 * list1 + cat1 * list2,
 freq_table,
 cens_ind)}

}
