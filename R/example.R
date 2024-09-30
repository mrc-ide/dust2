##' Load example generators from dust2.  These generators exist
##' primarily for the examples and documentation and are not (yet)
##' very interesting.  The examples will likely change as the package
##' evolves and some may be removed.
##'
##' All models exist as source code in the package; to view the `sir`
##' model you could write:
##'
##' ```r
##' file.show(system.file("examples/sir.cpp", package = "dust2"))
##' ```
##'
##' # `sir`
##'
##' A simple SIR (Susceptible-Infected-Recovered) compartmental model.
##' This model has parameters:
##'
##' * `N`: total population size
##' * `I0`: initial infected population size (when using
##'   [dust_system_set_state_initial])
##' * `beta`: per-contact rate of infection
##' * `gamma`: rate of recovery
##' * `exp_noise`: noise parameter used in the comparison to data
##'
##' The system will have compartments `S`, `I`, `R`, `cases_cumul` and
##' `cases_inc`
##'
##' # `sirode`
##'
##' The same model as `sir` but in continuous time, deterministically
##'
##' # `walk`
##'
##' A random walk in discrete time with Gaussian increments.  This
##' model has parameters:
##'
##' * `sd`: The standard deviation of the Gaussian update (per unit time)
##' * `len`: The number of independent walks
##' * random_initial`: Boolean, indicating if the initial position
##'   should be random (changes how [dust_system_set_state_initial]
##'   would initialise the system)
##'
##' @title Example generators
##'
##' @param name The name of the generator as a string; one of `sir`,
##'   `sirode` or `walk`.  See Details.
##'
##' @return A `dust_generator` object, which you might pass into
##'   `dust_system_create`
##'
##' @export
##' @examples
##' walk <- dust_example("walk")
##' walk()
##'
##' sys <- dust_system_create(walk(), list(sd = 1), 20)
##' y <- dust_system_simulate(sys, 0:50)
##' matplot(t(y[1, , ]), type = "l", col = "#0000ff55", lty = 1,
##'         xlab = "Time", ylab = "Value")
dust_example <- function(name) {
  examples <- list(sir = sir,
                   sirode = sirode,
                   walk = walk)
  examples[[match_value(name, names(examples))]]
}
