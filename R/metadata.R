## Each system will be marked up with decor-style comments, this allows
## us to describe a system without needing to try and parse anything,
## which is fraught at best with C++.
##
## The main things we need to know are:
##
## * what is the name of the class we are wrapping?
## * what should the name be of the system exported to R (if different)
## * what parameters does it accept (optional?)
## * does the system support comparison with data?
## * what data entries does it accept (optional?)
##
## Later, we'll want to do the same thing with GPU support if that
## looks different, and possibly with MPI support, though with those
## it might depend a bit on if we actually need much special there -
## it might be worth looking at the GPU stuff again fairly soon
## actually.
parse_metadata <- function(filename, call = NULL) {
  if (!file.exists(filename)) {
    cli::cli_abort("File '{filename}' does not exist", call = call)
  }
  data <- decor::cpp_decorations(files = filename)

  class <- parse_metadata_class(data, call)
  list(class = class,
       name = parse_metadata_name(data, call) %||% class,
       time_type = parse_metadata_time_type(data, call),
       has_compare = parse_metadata_has_compare(data, call),
       parameters = parse_metadata_parameters(data, call))
}


## All of the errors here could benefit from line number and context
## information, but it's not that important as users should never see
## these errors - systems will mostly be written by odin
parse_metadata_class <- function(data, call = NULL) {
  data <- find_attribute_value_single(data, "dust2::class", required = TRUE,
                                      call = call)
  if (length(data) != 1 || nzchar(names(data))) {
    cli::cli_abort(
      "Expected a single unnamed argument to '[[dust2::class()]]'",
      call = call)
  }
  if (!is.symbol(data[[1]])) {
    cli::cli_abort(
      "Expected an unquoted string argument to '[[dust2::class()]]'",
      call = call)
  }
  deparse(data[[1]])
}


parse_metadata_name <- function(data, call = NULL) {
  data <- find_attribute_value_single(data, "dust2::name", required = FALSE,
                                      call = call)
  if (is.null(data)) {
    return(NULL)
  }
  if (length(data) != 1 || nzchar(names(data))) {
    cli::cli_abort(
      "Expected a single unnamed argument to '[[dust2::name()]]'",
      call = call)
  }
  if (!is.symbol(data[[1]])) {
    cli::cli_abort(
      "Expected an unquoted string argument to '[[dust2::name()]]'",
      call = call)
  }
  deparse(data[[1]])
}


parse_metadata_time_type <- function(data, call = NULL) {
  data <- find_attribute_value_single(data, "dust2::time_type",
                                      required = TRUE, call = call)
  if (length(data) != 1 || nzchar(names(data))) {
    cli::cli_abort(
      "Expected a single unnamed argument to '[[dust2::time_type()]]'",
      call = call)
  }
  if (!is.symbol(data[[1]])) {
    cli::cli_abort(
      "Expected an unquoted string argument to '[[dust2::time_type()]]'",
      call = call)
  }
  value <- deparse(data[[1]])
  if (!(value %in% c("discrete", "continuous"))) {
    cli::cli_abort(
      paste("Expected argument to '[[dust2::time_type()]]' to be one of",
            "'discrete' or 'continuous'"),
      call = call)
  }
  value
}


parse_metadata_has_compare <- function(data, call = NULL) {
  data <- find_attribute_value_single(data, "dust2::has_compare",
                                      required = FALSE, call = call)
  if (is.null(data)) {
    return(FALSE)
  }
  if (length(data) != 0) {
    cli::cli_abort(
      "Expected no arguments to '[[dust2::has_compare()]]'",
      call = call)
  }
  TRUE
}


parse_metadata_parameters <- function(data, call = NULL) {
  res <- data$params[data$decoration == "dust2::parameter"]
  ok <- vlapply(res, function(x) {
   length(x) == 1 && !nzchar(names(x)[[1]]) && is.symbol(x[[1]])
  })
  if (!all(ok)) {
    cli::cli_abort(
      paste("Expected an unnamed unquoted string argument to",
            "'[[dust2::parameter()]]'"),
      call = call)
  }
  data_frame(name = vcapply(res, function(x) deparse(x[[1]])))
}


find_attribute_value_single <- function(data, name, required, call = NULL) {
  i <- data$decoration == name
  if (!any(i)) {
    if (required) {
      cli::cli_abort(
        "Attribute '[[{name}()]]' is required, but was not found",
        call = call)
    }
    return(NULL)
  }

  if (sum(i) > 1) {
    cli::cli_abort(
      "More than one '[[{name}()]]' attribute found",
      call = call)
  }

  data$params[[which(i)]]
}
