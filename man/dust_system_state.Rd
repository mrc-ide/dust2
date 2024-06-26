% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface.R
\name{dust_system_state}
\alias{dust_system_state}
\title{Extract system state}
\usage{
dust_system_state(sys)
}
\arguments{
\item{sys}{A \code{dust_system} object}
}
\value{
An array of system state.  If your system is ungrouped, then
this has two dimensions (state, particle).  If grouped, this has
three dimensions (state, particle, group)
}
\description{
Extract system state
}
\seealso{
\code{\link[=dust_system_set_state]{dust_system_set_state()}} for setting state and
\code{\link[=dust_system_set_state_initial]{dust_system_set_state_initial()}} for setting state to the
system-specific initial conditions.
}
