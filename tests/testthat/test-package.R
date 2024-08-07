test_that("can generate basic test package", {
  path <- withr::local_tempfile()
  create_test_package(path)
  expect_no_message(dust_package(path, quiet = TRUE))
  files <- dir(path, recursive = TRUE, include.dirs = FALSE,
               no.. = TRUE, all.files = TRUE)
  expect_setequal(
    files,
    c("DESCRIPTION", "NAMESPACE",
      "inst/dust/sir.cpp",
      "inst/dust/walk.cpp",
      "R/cpp11.R",
      "R/dust.R",
      "src/Makevars",
      "src/cpp11.cpp",
      "src/sir.cpp",
      "src/walk.cpp"))
})


test_that("can be chatty when creating package", {
  path <- withr::local_tempfile()
  create_test_package(path)
  res <- evaluate_promise(dust_package(path))
  expect_match(res$messages, "Working in package 'pkg' at", all = FALSE)
  expect_match(res$messages, "Found 2 systems", all = FALSE)
  expect_match(res$messages, "Wrote 'R/dust.R'", all = FALSE)
})


test_that("Fail to run if DESCRIPTION missing", {
  path <- withr::local_tempfile()
  create_test_package(path)
  unlink(file.path(path, "DESCRIPTION"))
  expect_error(
    package_validate_root(path),
    "Expected a file 'DESCRIPTION' at path '.+'")
})


test_that("Fail to run if NAMESPACE missing", {
  path <- withr::local_tempfile()
  create_test_package(path)
  unlink(file.path(path, "NAMESPACE"))
  expect_error(
    package_validate_root(path),
    "Expected a file 'NAMESPACE' at path '.+'")
})


test_that("Fail to run if DESCRIPTION lacks dependencies", {
  path <- withr::local_tempfile()
  create_test_package(path)
  desc <- readLines(file.path(path, "DESCRIPTION"))
  writeLines(desc[!grepl("^Imports:", desc)], file.path(path, "DESCRIPTION"))
  expect_error(
    package_validate_root(path),
    "Expected package 'dust2' as 'Imports' in DESCRIPTION")
})


test_that("Fail to run if no dust files found", {
  path <- withr::local_tempfile()
  create_test_package(path)
  unlink(dir(file.path(path, "inst", "dust"), full.names = TRUE))
  expect_error(
    dust_package(path, quiet = TRUE),
    "No dust files found in 'inst/dust'")
})


test_that("validate destination notices existing C++ code", {
  path <- withr::local_tempfile()
  create_test_package(path)

  path_cpp <- file.path(path, "src", "walk.cpp")
  msg <- "File '.+\\.cpp' does not look like it was created by dust - stopping"

  file.create(path_cpp)
  expect_error(
    package_validate_destination(path, c("sir.cpp", "walk.cpp")),
    msg)

  writeLines("// some actual content", path_cpp)
  expect_error(
    package_validate_destination(path, c("sir.cpp", "walk.cpp")),
    msg)

  writeLines("// Generated by dust2", path_cpp)
  expect_no_error(
    package_validate_destination(path, c("sir.cpp", "walk.cpp")))
})


test_that("validate destination notices existing R code", {
  path <- withr::local_tempfile()
  create_test_package(path)

  path_r <- file.path(path, "R", "dust.R")
  msg <- "File '.+\\.R' does not look like it was created by dust - stopping"

  file.create(path_r)
  expect_error(
    package_validate_destination(path, character()),
    msg)

  writeLines("## some actual content", path_r)
  expect_error(
    package_validate_destination(path, character()),
    msg)

  writeLines("## Generated by dust", path_r)
  expect_silent(
    package_validate_destination(path, character()))
})


test_that("Validate NAMESPACE has correct useDynLib call", {
  path <- withr::local_tempfile()
  create_test_package(path)
  path_ns <- file.path(path, "NAMESPACE")
  txt <- readLines(path_ns)

  expect_no_error(
    package_validate_namespace(path_ns, "pkg"))

  err <- expect_error(
    package_validate_namespace(path_ns, "other"),
    "Did not find a 'useDynLib()' call for 'other' in NAMESPACE",
    fixed = TRUE)
  expect_match(
    conditionMessage(err),
    "Found unexpected 'useDynLib()' call for 'pkg'",
    fixed = TRUE)
  expect_match(err$body, "without Roxygen", all = FALSE)

  ## Nice error if using roxygen
  writeLines(c("# Generated by roxygen2: do not edit by hand", txt), path_ns)
  err <- expect_error(
    package_validate_namespace(path_ns, "other"),
    "Did not find a 'useDynLib()' call for 'other' in NAMESPACE",
    fixed = TRUE)
  expect_match(err$body, "It looks like you are using Roxygen", all = FALSE)

  writeLines(gsub('"', "", txt), path_ns)
  expect_no_error(
    package_validate_namespace(path_ns, "pkg"))
  expect_error(
    package_validate_namespace(path_ns, "other"),
    "Did not find a 'useDynLib()' call for 'other' in NAMESPACE",
    fixed = TRUE)

  file.create(path_ns)
  err <- expect_error(
    package_validate_namespace(path_ns, "other"),
    "Did not find a 'useDynLib()' call for 'other' in NAMESPACE",
    fixed = TRUE)
  expect_match(err$body, "without Roxygen", all = FALSE)
})


test_that("Validate openmp support", {
  text <- read_lines(dust2_file("template/Makevars.pkg"))
  expect_no_error(package_validate_makevars_openmp(text))
  msg <- "Package has a 'src/Makevars' but no openmp flags support"
  expect_error(
    package_validate_makevars_openmp(""),
    msg)
  expect_error(
    package_validate_makevars_openmp("PKG_CXXFLAGS=$(SHLIB_OPENMP_CXXFLAGS)"),
    msg)
  expect_error(
    package_validate_makevars_openmp("PKG_LIBS=$(SHLIB_OPENMP_CXXFLAGS)"),
    msg)
})


test_that("Validate openmp support in real package", {
  path <- withr::local_tempfile()
  create_test_package(path)
  file.create(file.path(path, "src/Makevars"))
  expect_error(
    dust_package(path, quiet = TRUE),
    "Package has a 'src/Makevars' but no openmp flags support")
})


test_that("guide user to sensible package name", {
  path <- withr::local_tempfile()
  create_test_package(path)

  path_descr <- file.path(path, "DESCRIPTION")
  descr <- sub("Package: pkg", "Package: my.pkg", readLines(path_descr),
               fixed = TRUE)
  writeLines(descr, path_descr)

  path_namespace <- file.path(path, "NAMESPACE")
  namespace <- sub("pkg", "my.pkg", readLines(path_namespace),
               fixed = TRUE)
  writeLines(namespace, path_namespace)

  expect_error(
    dust_package(path, quiet = TRUE),
    "Package name must not contain '.' or '_' (found 'my.pkg')",
    fixed = TRUE)
})
