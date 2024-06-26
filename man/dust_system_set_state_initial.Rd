% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface.R
\name{dust_system_set_state_initial}
\alias{dust_system_set_state_initial}
\title{Set system state to initial conditions}
\usage{
dust_system_set_state_initial(sys)
}
\arguments{
\item{sys}{A \code{dust_system} object}
}
\value{
Nothing, called for side effects only
}
\description{
Set system state from a system's initial conditions.  This may depend
on the current time.
}
