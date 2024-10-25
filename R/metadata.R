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
  time_type <- parse_metadata_time_type(data, call)
  list(class = class,
       name = parse_metadata_name(data, call) %||% class,
       time_type = time_type,
       default_dt = parse_metadata_default_dt(data, time_type, call),
       has_compare = parse_metadata_has_compare(data, call),
       has_adjoint = parse_metadata_has_adjoint(data, call),
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
  valid <- c("discrete", "continuous", "mixed")
  if (!(value %in% valid)) {
    cli::cli_abort(
      paste("Expected argument to '[[dust2::time_type()]]' to be one of",
            "{squote(valid)}"),
      call = call)
  }
  value
}


parse_metadata_default_dt <- function(data, time_type, call = NULL) {
  data <- find_attribute_value_single(data, "dust2::default_dt",
                                      required = FALSE, call = call)
  if (is.null(data)) {
    return(if (time_type == "discrete") 1 else NULL)
  }
  if (time_type == "continuous") {
    cli::cli_abort(
      "Can't use '[[dust::default_dt()]]' with continuous-time systems",
      call = call)
  }
  if (length(data) != 1 || nzchar(rlang::names2(data))) {
    cli::cli_abort(
      "Expected a single unnamed argument to '[[dust2::default_dt()]]'",
      call = call)
  }
  if (!is.numeric(data[[1]])) {
    cli::cli_abort(
      "Expected a numerical argument to '[[dust2::default_dt()]]'",
      call = call)
  }
  value <- data[[1]]
  check_dt(value, name = "[[dust2::default_dt()]]", call = call)
  value
}


parse_metadata_has_compare <- function(data, call = NULL) {
  parse_metadata_has_feature("compare", data, call)
}


parse_metadata_has_adjoint <- function(data, call = NULL) {
  parse_metadata_has_feature("adjoint", data, call)
}


parse_metadata_has_feature <- function(name, data, call = NULL) {
  attribute <- sprintf("dust2::has_%s", name)
  data <- find_attribute_value_single(data, attribute,
                                      required = FALSE, call = call)
  if (is.null(data)) {
    return(FALSE)
  }
  if (length(data) != 0) {
    cli::cli_abort(
      "Expected no arguments to '[[{attribute}()]]'",
      call = call)
  }
  TRUE
}


parse_metadata_parameters <- function(data, call = NULL) {
  data <- data[data$decoration == "dust2::parameter", ]
  res <- lapply(seq_len(nrow(data)), function(i) {
    tryCatch(
      parse_metadata_parameter_single(data$params[[i]]),
      error = function(e) {
        line <- data$line[[i]]
        cli::cli_abort(
          "Invalid arguments to [[dust2::parameter()]] on line {line}",
          call = call, parent = e)
      })
  })
  ret <- data_frame(
    name = vcapply(res, "[[", "name"),
    type = vcapply(res, "[[", "type"),
    constant = vlapply(res, "[[", "constant"),
    required = vlapply(res, "[[", "required"),
    rank = viapply(res, "[[", "rank"),
    min = vnapply(res, "[[", "min"),
    max = vnapply(res, "[[", "max"))


  if (anyDuplicated(ret$name)) {
    dups <- unique(ret$name[duplicated(ret$name)])
    cli::cli_abort(
      paste("Duplicate parameter names across '[[dust2::paramneter()]]'",
            "entries: {squote(dups)}"),
      call = call)
  }

  for (v in c("constant", "required", "rank")) {
    i <- is.na(ret[[v]])
    if (any(i) && !all(i)) {
      cli::cli_abort(
        c("Only some '[[dust2::parameter()]]' entries have arguments for '{v}'",
          i = paste("If any of your parameters include this argument, then",
                    "all parameters must have it, otherwise it's confusing")),
        call = call)
    }
  }

  ret[!vlapply(ret, function(x) all(is.na(x)))]
}


parse_metadata_parameter_single <- function(args) {
  fn <- function(name, type = "real_type", required = NA, constant = NA,
                 min = NA_real_, max = NA_real_,
                 rank = NA_integer_) {
  }

  call <- rlang::call2("fn", !!!args)
  d <- tryCatch(
    as.list(rlang::call_match(call, fn, defaults = TRUE)[-1]),
    error = function(e) {
      cli::cli_abort(
        "Invalid arguments to attribute [[dust2::parameter()]]",
        call = call, parent = e)
    })
  name_arg <- if (is.null(args$name)) "first argument" else "'name'"
  if (rlang::is_missing(d$name) || !is.name(d$name)) {
    cli::cli_abort("Expected {name_arg} to be an unquoted string",
                   call = call)
  }
  d$name <- as.character(d$name)
  match_value(d$type, c("real_type", "int", "bool"), "type")
  if (d$type == "bool") {
    for (v in c("min", "max")) {
      if (!is.na(d[[v]])) {
        cli::cli_abort(
          "Argument '{v}' does not make sense with 'type = \"bool\"'")
      }
    }
  } else {
    if (!is.na(d$min)) {
      assert_scalar_numeric(d$min, "min")
    }
    if (!is.na(d$max)) {
      assert_scalar_numeric(d$max, "max")
    }
    if (isTRUE(d$min > d$max)) {
      cli::cli_abort("'min' is greater than 'max'")
    }
  }
  if (!is.na(d$rank)) {
    assert_scalar_size(d$rank, name = "rank")
  }
  if (!is.na(d$required)) {
    assert_scalar_logical(d$required, name = "required")
  }
  if (!is.na(d$constant)) {
    assert_scalar_logical(d$constant, name = "constant")
  }
  d$rank <- as.integer(d$rank)
  d
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
