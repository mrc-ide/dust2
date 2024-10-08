% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface-likelihood.R
\name{dust_likelihood_rng_state}
\alias{dust_likelihood_rng_state}
\alias{dust_likelihood_set_rng_state}
\title{Get filter RNG state}
\usage{
dust_likelihood_rng_state(obj)

dust_likelihood_set_rng_state(obj, rng_state)
}
\arguments{
\item{obj}{A \code{dust_filter} object, created by
\link{dust_filter_create} or a \code{dust_unfilter} object created by
\link{dust_unfilter_create}}

\item{rng_state}{A raw vector of random number generator state,
returned by \code{dust_likelihood_rng_state}}
}
\value{
A raw vector, this could be quite long.  Later we will
describe how you might reseed a filter or system with this state.
}
\description{
Get random number generator (RNG) state from the particle filter.
}
