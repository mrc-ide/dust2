`%||%` <- function(x, y) { # nolint
  if (is.null(x)) y else x
}


set_names <- function(x, nms) {
  names(x) <- nms
  x
}


vlapply <- function(...) {
  vapply(..., FUN.VALUE = TRUE)
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
