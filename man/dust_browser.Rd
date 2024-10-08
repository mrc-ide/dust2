% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browser.R
\name{dust_browser}
\alias{dust_browser}
\alias{dust_browser_enabled}
\alias{dust_browser_verbosity}
\alias{dust_browser_continue}
\title{The dust debugger}
\usage{
dust_browser_enabled(value = TRUE)

dust_browser_verbosity(level)

dust_browser_continue()
}
\arguments{
\item{value}{Logical, \code{TRUE} for where the debugger should be
enabled, \code{FALSE} otherwise.}

\item{level}{The verbosity level, as a string.  This must be one
of the values \code{quiet} (prevents informational messages),
\code{normal} (prints a single line on entry) and \code{verbose} (prints
several informational messages on entry).  The default is \code{normal}.}
}
\value{
Both \code{dust_browser_enabled} and \code{dust_browser_verbosity}
return the previous value of the option they are setting.
}
\description{
Control browser-based debugging of dust models.  This help page documents
three funtions that can be used to control if and how the browser
is enabled.  You can't enter the debugger from any of these
functions; it is only enabled if present in your C++ code (or if
using \code{odin2} if you have enabled it).
}
\details{
dust2 includes an extremely simple debugging system, and if you
are reading this message, there's a good chance you are inside it.
It is built on top of R's \code{\link[=browser]{browser()}} and so all the usual tips,
tricks and issues for working with this apply.  We recommend
setting the R option \code{browserNLdisabled = TRUE} to avoid surprises
from presssing \verb{<enter>}.
\itemize{
\item You can press \code{n} or \code{c} to proceed to the next enabled iteration
\item You can press \code{Q} to quit the browser (this will end up as an
error by the time you have control back)
}

These commands are established by \code{browser}, and can't be
disabled.  This means that if you have an variable called \code{n} you
will need to work with it as \code{(n)} (i.e. in parentheses).  This
applies to all of \code{browser}'s command variables (\code{c}, \code{f}, \code{n}, \code{s},
\code{r} and \code{Q}); please see \link{browser} for more information.

By the time the environment has been created, some variables from
your model will have been copied into the environment; you can see
these by running \code{ls()} and write expressions involving these
objects.  Changes that you make in R are \emph{not} (currently)
propagated back into the running system.

If you enable the debugger, you may have very many iterations to
get through before control is returned back to the console.  You
can run \code{dust_debug_continue()} to prevent entry into the debugger
until control is passed back to you; this means the time
series will run to completion and then the next time you run the
system the debugger will be triggered again.  Alternatively, you
can run \code{dust_debug_enabled(FALSE)} to disable all calls to the
debugger.
}
