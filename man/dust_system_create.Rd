% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface.R
\name{dust_system_create}
\alias{dust_system_create}
\title{Create a dust object}
\usage{
dust_system_create(
  generator,
  pars,
  n_particles,
  n_groups = 0,
  time = 0,
  dt = NULL,
  ode_control = NULL,
  seed = NULL,
  deterministic = FALSE
)
}
\arguments{
\item{generator}{A system generator object, with class
\code{dust_system_generator}}

\item{pars}{A list of parameters.  The format of this will depend on the
system.  If \code{n_groups} is 1 or more, then this must be a list of length
\code{n_groups} where each element is a list of parameters for your system.}

\item{n_particles}{The number of particles to create.}

\item{n_groups}{Optionally, the number of parameter groups}

\item{time}{The initial time, defaults to 0}

\item{dt}{The time step for discrete time systems, defaults to 1
if not given.  It is an error to provide a non-NULL argument
with continuous-time systems.}

\item{ode_control}{The ODE integration control for continuous time
systems.  Defaults to the default return of \link{dust_ode_control}.
It is an error to provide this with discrete-time systems.}

\item{seed}{Optionally, a seed.  Otherwise we respond to R's RNG seed on
initialisation.}

\item{deterministic}{Logical, indicating if the system should be
allocated in deterministic mode.}
}
\value{
A \code{dust_system} object, with opaque format.
}
\description{
Create a dust system object from a system generator.  This allocates a
system and sets an initial set of parameters.  Once created you can use
other dust functions to interact with it.
}
