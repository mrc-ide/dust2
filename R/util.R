`%||%` <- function(x, y) { # nolint
  if (is.null(x)) y else x
}


set_names <- function(x, nms) {
  if (length(nms) == 1 && length(x) != 1) {
    if (is.null(x)) {
      return(x)
    }
    nms <- rep_len(nms, length(x))
  }
  names(x) <- nms
  x
}


vlapply <- function(...) {
  vapply(..., FUN.VALUE = TRUE)
}


viapply <- function(...) {
  vapply(..., FUN.VALUE = 1L)
}


vnapply <- function(...) {
  vapply(..., FUN.VALUE = 1)
}


vcapply <- function(...) {
  vapply(..., FUN.VALUE = "")
}


data_frame <- function(...) {
  data.frame(..., stringsAsFactors = FALSE, check.names = FALSE)
}


dust2_file <- function(path) {
  system.file(path, mustWork = TRUE, package = "dust2")
}


glue_whisker <- function(template, data) {
  glue::glue_data(data, template, .open = "{{", .close = "}}", .trim = FALSE)
}


read_lines <- function(path) {
  paste(readLines(path), collapse = "\n")
}


dir_create <- function(path) {
  for (p in path) {
    dir.create(p, showWarnings = FALSE, recursive = TRUE)
  }
}


is_directory <- function(path) {
  file.info(path, extra_cols = FALSE)$isdir
}


squote <- function(x) {
  sprintf("'%s'", x)
}


writelines_if_changed <- function(text, workdir, path, quiet) {
  path_full <- file.path(workdir, path)
  skip <- file.exists(path_full) && same_content(path_full, text)
  if (skip) {
    if (!quiet) {
      cli::cli_alert_info("'{path}' is up to date")
    }
  } else {
    writeLines(text, path_full)
    if (!quiet) {
      cli::cli_alert_success("Wrote '{path}'")
    }
  }
}


same_content <- function(path, text) {
  identical(read_lines(path), paste(as.character(text), collapse = "\n"))
}


protect <- function(fn, on_error) {
  force(fn)
  force(on_error)
  function(...) {
    tryCatch(fn(...), error = function(e) on_error)
  }
}


tail_errors <- function(x, n = 5) {
  if (length(x) > n) {
    x <- c(
      x[seq_len(n - 1)],
      sprintf("...and %d other errors", length(x) - (n - 1)))
  }
  set_names(x, "x")
}

fmod <- function(n, m) {
  (floor(n / m) * m) - n
}


envvar_is_truthy <- function(name) {
  tolower(Sys.getenv(name, "false")) %in% c("t", "true", "1")
}


last <- function(x) {
  x[[length(x)]]
}


drop_last <- function(x, n = 1) {
  len <- length(x)
  x[if (len < n) integer(0) else seq_len(len - n)]
}


dim2 <- function(x) {
  dim(x) %||% length(x)
}


rank_description <- function(rank) {
  if (rank == 0) {
    "scalar"
  } else if (rank == 1) {
    "vector"
  } else if (rank == 2) {
    "matrix"
  } else {
    sprintf("%d-dimensional array", rank)
  }
}


deparse_df <- function(df, indent) {
  if (is.null(df)) {
    return("NULL")
  }
  values <- vcapply(df, function(x) {
    str <- paste(vcapply(x, deparse), collapse = ", ")
    if (length(x) == 1) str else sprintf("c(%s)", str)
  })
  indent <- strrep(" ", indent)
  ret <- c("data.frame(",
           sprintf("%s%s = %s", indent, names(df), values))
  suffix <- rep(c("\n", ",\n", ")"), c(1, ncol(df) - 1, 1))
  paste0(ret, suffix, collapse = "")
}


disable_write <- function(cls, name, call = parent.frame()) {
  if (is.character(name)) {
    detail <- c(x = "Was trying to modify {squote(name)}")
  } else {
    detail <- NULL
  }
  cli::cli_abort(
    c("Cannot write to '{cls}' objects, they are read-only",
      detail),
    call = NULL)
}
