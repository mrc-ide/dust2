##' Creates or updates the generated code for a set of dust systems in
##' a package. The user-provided code is assumed to be in `inst/dust`
##' as a series of C++ files; a file `inst/dust/x.cpp` will be
##' transformed into a file `src/x.cpp`.
##'
##' Classes used within a package must be distinct; typically these
##' will match the filenames.
##'
##' We add "cpp11 attributes" to the created functions, and will run
##' [cpp11::cpp_register()] on them once the generated code
##' has been created.
##'
##' Your package needs a `src/Makevars` file to enable OpenMP (if your
##' system supports it). If it is not present then a suitable Makevars
##' will be written, containing
##'
##' ```
##' PKG_CXXFLAGS=$(SHLIB_OPENMP_CXXFLAGS)
##' PKG_LIBS=$(SHLIB_OPENMP_CXXFLAGS)
##' ```
##'
##' following "Writing R Extensions" (see section "OpenMP support").
##' If your package does contain a `src/Makevars` file we do not
##' attempt to edit it but will error if it looks like it does not
##' contain these lines or similar.
##'
##' You also need to make sure that your package loads the dynamic
##' library; if you are using roxygen, then you might create a file
##' (say, `R/zzz.R`) containing
##'
##' ```
##' #' @useDynLib packagename, .registration = TRUE
##' NULL
##' ```
##'
##' substituting `packagename` for your package name as
##' appropriate. This will create an entry in `NAMESPACE`.
##'
##' @param path Path to the package
##'
##' @param quiet Passed to `cpp11::cpp_register`, if `TRUE` suppresses
##'   informational notices about updates to the cpp11 files.  If
##'   `NULL`, uses the value of the environment variable `DUST_QUIET`
##'   if set or `FALSE` otherwise.
##'
##' @title Create dust system in package
##'
##' @return Nothing, this function is called for its side effects
##'
##' @export
dust_package <- function(path, quiet = NULL) {
  call <- environment()
  pkg <- package_validate_root(path, call)
  path_dust <- file.path(path, "inst/dust")

  quiet <- dust_quiet(quiet)
  if (!quiet) {
    cli::cli_alert_info("Working in package '{pkg}' at '{path}'")
  }

  files <- dir(path_dust, pattern = "\\.cpp$")
  if (length(files) == 0L) {
    cli::cli_abort("No dust files found in 'inst/dust'")
  }
  if (!quiet) {
    cli::cli_alert_info("Found {length(files)} system{?s}")
  }

  package_validate_destination(path, files, call)

  data <- lapply(file.path(path_dust, files), package_generate)

  dir_create(file.path(path, c("src", "R")))
  for (i in seq_along(files)) {
    writelines_if_changed(
      c(dust_header("//"), data[[i]]$cpp),
      path,
      file.path("src", basename(files[[i]])),
      quiet)
  }

  code_r <- c(dust_header("##"), vcapply(data, "[[", "r"))
  writelines_if_changed(code_r, path, "R/dust.R", quiet)

  if (file.exists(file.path(path, "src/Makevars"))) {
    makevars <- read_lines(file.path(path, "src/Makevars"))
    package_validate_makevars_openmp(makevars)
  } else {
    makevars <- read_lines(dust2_file("template/Makevars.pkg"))
    writelines_if_changed(makevars, path, "src/Makevars", quiet)
  }

  cpp11::cpp_register(path, quiet = quiet)

  invisible()
}


package_validate_root <- function(path, call) {
  paths <- c("DESCRIPTION", "NAMESPACE")
  for (p in paths) {
    if (!file.exists(file.path(path, p))) {
      cli::cli_abort("Expected a file '{p}' at path '{path}'",
                     call = call)
    }
  }

  desc <- pkgload::pkg_desc(path)
  pkg <- desc$get_field("Package")
  if (pkg != "dust2") {
    package_validate_has_dep(desc, "dust2", "Imports")
    package_validate_has_dep(desc, "dust2", "LinkingTo")
  }
  package_validate_has_dep(desc, "cpp11", "LinkingTo")
  package_validate_has_dep(desc, "monty", "LinkingTo")

  name <- desc$get_field("Package")
  if (grepl("[._]+", name)) {
    cli::cli_abort(
      "Package name must not contain '.' or '_' (found '{name}')",
      call = call)
  }

  package_validate_namespace(file.path(path, "NAMESPACE"), name, call)

  pkg
}


package_validate_has_dep <- function(desc, package, type, call = NULL) {
  if (!desc$has_dep(package, type)) {
    cli::cli_abort("Expected package '{package}' as '{type}' in DESCRIPTION",
                   call = call)
  }
}


package_validate_destination <- function(path, files, call) {
  check <- c(file.path(path, "src", files),
             file.path(path, "R", "dust.R"))
  for (f in check[file.exists(check)]) {
    if (!isTRUE(grepl("^(//|##) Generated by dust", readLines(f, 1)))) {
      cli::cli_abort(
        "File '{f}' does not look like it was created by dust - stopping",
        call = call)
    }
  }
}

package_validate_namespace <- function(path, name, call) {
  ## There are two seemingly reasonable ways of sniffing use of
  ## roxygen, both are probably reasonable.  We can look in
  ## DESCRIPTION for RoxygenNote, or we can look in NAMESPACE for the
  ## "generated by roxygen2" header.  We go for the latter here
  ## because it localises things into this function and because
  ## roxygen will not update the NAMESPACE if the header is not there.
  ##
  ## any() here copes with the case of an entirely empty NAMESPACE
  ## file, converting logical() to FALSE
  using_roxygen <- any(grepl("roxygen2", readLines(path, n = 1)))
  exprs <- as.list(parse(file = path))
  package_validate_namespace_usedynlib(exprs, name, using_roxygen, call)
}


package_validate_namespace_usedynlib <- function(exprs, name, using_roxygen,
                                                 call) {
  unexpected <- character()
  for (e in exprs) {
    if (rlang::is_call(e, "useDynLib")) {
      lib <- e[[2]]
      if (is.name(lib)) {
        lib <- deparse(lib)
      }
      if (identical(lib, name)) {
        return()
      }
      unexpected <- c(unexpected, lib)
    }
  }
  if (using_roxygen) {
    hint <- paste(
      "It looks like you are using Roxygen, so you should add",
      "{.code #' @useDynLib {name}, .registration = TRUE}",
      "into one of your source files and then run",
      "{.run devtools::document()} to regenerate your",
      "NAMESPACE file")
  } else {
    hint <- paste(
      "It looks like you are managing your NAMESPACE file without Roxygen,",
      "so you should add {.code useDynLib({name}, .registration = TRUE)}",
      "into your NAMESPACE file")
  }
  if (length(unexpected) > 0) {
    warn <-
      c(x = "Found unexpected 'useDynLib()' call for {squote(unexpected)}")
  } else {
    warn <- NULL
  }

  cli::cli_abort(
    c("Did not find a 'useDynLib()' call for '{name}' in NAMESPACE",
      warn,
      i = hint),
    call = call)
}


package_validate_makevars_openmp <- function(text, call) {
  ok <- grepl("PKG_CXXFLAGS", text) &&
    grepl("PKG_LIBS", text) &&
    grepl("SHLIB_OPENMP_CXXFLAGS", text)
  if (!ok) {
    cli::cli_abort(
      "Package has a 'src/Makevars' but no openmp flags support",
      call = call)
  }
}


package_generate <- function(filename, call) {
  config <- parse_metadata(filename, call = call)
  system <- read_lines(filename)
  data <- dust_template_data(config$name, config$class, config$time_type,
                             config$default_dt)
  list(r = substitute_dust_template(data, "package.R"),
       cpp = dust_generate_cpp(system, config, data))
}


dust_header <- function(comment) {
  sprintf("%s Generated by dust2 (version %s) - do not edit",
          comment, utils::packageVersion("dust2"))
}
