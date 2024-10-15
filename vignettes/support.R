knitr::opts_chunk$set(
  error = FALSE,
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5)

lang_output <- function(x, lang) {
  cat(c(sprintf("```%s", lang), x, "```"), sep = "\n")
}
cc_output <- function(x) lang_output(x, "cc")
r_output <- function(x) lang_output(x, "r")
plain_output <- function(x) lang_output(x, "plain")
