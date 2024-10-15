# dust2 <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->
[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![R-CMD-check](https://github.com/mrc-ide/dust2/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/mrc-ide/dust2/actions/workflows/R-CMD-check.yaml)
[![codecov.io](https://codecov.io/github/mrc-ide/dust2/coverage.svg?branch=main)](https://codecov.io/github/mrc-ide/dust2?branch=main)
<!-- badges: end -->

The `dust2` package provides an engine for running dynamical systems in discrete or continuous time and where the processes are stochastic or deterministic.  We focus on Markov models where the problem reduces to describing how model state changes as a function of its current state (and possibly time) but without reference to where it has come from.  Superficially, the problem is not very hard (see `vignette("design")`), but `dust2` takes care of many practical and bookkeeping details such as:

* Running systems in parallel on multi-core machines, even those involving random numbers
* Providing useful verbs for efficiently working with groups of simulations (different parameters, starting conditions or stochastic realisations)
* Comparing simulations to time-series of data, and implementing sequential Monte Carlo methods such as a bootstrap particle filter

## Get started

* `vignette("design")` describes the problems `dust` tries to solve ([read on package website](https://mrc-ide.github.io/dust/articles/design.html))
* `vignette("dust2")` describes `dust` by example, showing two simple systems and the methods that can drive them ([read on package website](https://mrc-ide.github.io/dust/articles/dust.html))
* [`odin2`](https://mrc-ide.github.io/odin2) is the way most `dust2` systems are written
* The [odin & monty book](https://mrc-ide.github.io/odin-monty) shows how the package works in context
* `vignette("writing")` shows how to write a `dust2` system by hand in C++
* `vignette("packaging")` describes how to package a system for easy reuse

## Roadmap

This package is a ground-up rewrite of [`dust`](https://mrc-ide.github.io/dust) and will eventually become version 2.0.0 of `dust`, which we will then release to CRAN.  It exists separately for now to facilitate development and use alongside the original `dust`, and is being developed in parallel with [`odin2`](https://mrc-ide.github.io/odin2) and [`monty`](https://mrc-ide.github.io/monty) (previously [`mcstate`](https://mrc-ide.github.io/mcstate)).  Some of the functionality here was originally found in `mcstate` and some of the previous version of `dust` can now be found in `monty` (e.g., the random number library).

## Installation

To install `dust2`:

```r
remotes::install_github("mrc-ide/dust2", upgrade = FALSE)
```

## License

MIT © Imperial College of Science, Technology and Medicine
