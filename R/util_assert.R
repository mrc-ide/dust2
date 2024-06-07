assert_scalar <- function(x, name = deparse(substitute(x)), arg = name,
                          call = NULL)  {
  if (length(x) != 1) {
    cli::cli_abort(
      c("'{name}' must be a scalar",
        i = "{name} has length {length(x)}"),
      call = call, arg = arg)
  }
  invisible(x)
}


assert_character <- function(x, name = deparse(substitute(x)),
                             arg = name, call = NULL) {
  if (!is.character(x)) {
    cli::cli_abort("Expected '{name}' to be character", call = call, arg = arg)
  }
  invisible(x)
}


assert_numeric <- function(x, name = deparse(substitute(x)),
                             arg = name, call = NULL) {
  if (!is.numeric(x)) {
    cli::cli_abort("Expected '{name}' to be numeric", call = call, arg = arg)
  }
  invisible(x)
}


assert_integer <- function(x, name = deparse(substitute(x)),
                           arg = name, call = NULL) {
  ## See the help in rlang; we might be being too strict here.
  if (!rlang::is_integerish(x)) {
    cli::cli_abort("Expected '{name}' to be integer", call = call, arg = arg)
  }
  if (!is.integer(x)) {
    x <- as.integer(round(x))
  }
  invisible(x)
}


assert_nonmissing <- function(x, name = deparse(substitute(x)),
                              arg = name, call = NULL) {
  if (anyNA(x)) {
    cli::cli_abort("Expected '{name}' to be non-NA", arg = arg, call = call)
  }
  invisible(x)
}


assert_scalar_character <- function(x, name = deparse(substitute(x)),
                                    arg = name, call = NULL) {
  assert_scalar(x, name, arg = arg, call = call)
  assert_character(x, name, arg = arg, call = call)
  assert_nonmissing(x, name, arg = arg, call = call)
}


assert_scalar_numeric <- function(x, name = deparse(substitute(x)),
                                    arg = name, call = NULL) {
  assert_scalar(x, name, arg = arg, call = call)
  assert_numeric(x, name, arg = arg, call = call)
  assert_nonmissing(x, name, arg = arg, call = call)
}


assert_scalar_integer <- function(x, name = deparse(substitute(x)),
                                  arg = name, call = call) {
  assert_scalar(x, name, arg = arg, call = call)
  assert_nonmissing(x, name, arg = arg, call = call)
  assert_integer(x, name, arg = arg, call = call)
}


assert_scalar_size <- function(x, allow_zero = TRUE,
                               name = deparse(substitute(x)),
                               arg = name, call = call) {
  assert_scalar_integer(x, name = name, arg = arg, call = call)
  min <- if (allow_zero) 0 else 1
  if (x < min) {
    cli::cli_abort("'{name}' must be at least {min}")
  }
}


assert_length <- function(x, len, name = deparse(substitute(x)), arg = name,
                          call = parent.frame()) {
  if (length(x) != len) {
    cli::cli_abort(
      "Expected '{name}' to have length {len}, but was length {length(x)}",
      arg = arg, call = call)
  }
  invisible(x)
}


assert_raw_vector <- function(x, len, name = deparse(substitute(x)),
                              arg = name, call = parent.frame()) {
  if (!is.raw(x)) {
    cli::cli_abort("'{name}' must be a raw vector", arg = arg, call = call)
  }
  assert_length(x, len)
  invisible(x)
}


assert_list <- function(x, name = deparse(substitute(x)), arg = name,
                        call = NULL) {
  if (!is.list(x)) {
    cli::cli_abort("Expected '{name}' to be a list",
                   arg = arg, call = call)
  }
  invisible(x)
}
