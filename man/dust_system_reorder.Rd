% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface.R
\name{dust_system_reorder}
\alias{dust_system_reorder}
\title{Reorder states}
\usage{
dust_system_reorder(sys, index)
}
\arguments{
\item{sys}{A \code{dust_system} object}

\item{index}{The parameter ordering.  For an ungrouped system this
is a vector where each element is the parameter index (if
element \code{i} is \code{j} then after reordering the \code{i}th particle will
have the state previously used by \code{j}).  All elements must lie
in \verb{[1, n_particles]}, repetition of an index is allowed (so
that many new particles may have the state as one old particle).
If the system is grouped, \code{index} must be a matrix with
\code{n_particles} rows and \code{n_groups} columns, with each column
corresponding to the reordering for a group.}
}
\value{
Nothing, called for side effects only.
}
\description{
Reorder states within a system.  This function is primarily used
for debugging and may be removed from the interface if it is not
generally useful.
}
