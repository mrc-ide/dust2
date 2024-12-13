test_that("can construct template data", {
  expect_equal(
    dust_template_data("foo", "foo", "discrete", FALSE, FALSE, NULL, 1),
    list(name = "foo", class = "foo",
         time_type_property = "discrete", time_type = "discrete",
         has_compare = "FALSE", has_adjoint = "FALSE", parameters = "NULL",
         default_dt = "1", package = "foo",
         linking_to = "cpp11, dust2, monty", cpp_std = NULL,
         compiler_options = ""))
})

test_that("can mangle a package name", {
  res <- dust_template_data("foo", "bar", "discrete", FALSE, FALSE, NULL, 1,
                            mangle = "abc")
  expect_equal(res$name, "foo")
  expect_equal(res$class, "bar")
  expect_equal(res$package, "fooabc")
})

test_that("can add linking to", {
  res <- dust_template_data("foo", "foo", "discrete", FALSE, FALSE, NULL, 1,
                            linking_to = "baz")
  expect_equal(res$linking_to, "cpp11, dust2, monty, baz")
})

test_that("can link to multiple packages and avoid relinking dust", {
  res <- dust_template_data("foo", "foo", "discrete", FALSE, FALSE, NULL, 1,
                            linking_to = c("x", "dust2", "y"))
  expect_equal(res$linking_to, "cpp11, dust2, monty, x, y")
})

test_that("can pass compiler options", {
  res <- dust_template_data("foo", "foo", "discrete", FALSE, FALSE, NULL, 1,
                            optimisation_level = "none",
                            compiler_options = "-Xf")
  expect_equal(res$compiler_options, "-Xf -O0")
})

test_that("can control the C++ standard", {
  res <- dust_template_data("foo", "foo", "discrete", FALSE, FALSE, NULL, 1,
                            cpp_std = "c++14")
  expect_equal(res$cpp_std, "c++14")
})

test_that("can set a default dt", {
  res <- dust_template_data("foo", "foo", "discrete", FALSE, FALSE, NULL, 0.25)
  expect_equal(res$default_dt, "0.25")
  res <- dust_template_data("foo", "foo", "continuous", FALSE, FALSE, NULL,
                            NULL)
  expect_equal(res$default_dt, "NULL")
  res <- dust_template_data("foo", "foo", "mixed", FALSE, FALSE, NULL, NULL)
  expect_equal(res$default_dt, "NULL")
})

test_that("can add parameter information", {
  df <- data.frame(name = c("a", "b"))
  res <- dust_template_data("foo", "foo", "discrete", FALSE, FALSE, df, 0.25)
  expect_equal(res$parameters, deparse_df(df, 4))
})


test_that("can cope with class name containing underscore", {
  res <- dust_template_data(
    "my_thing", "my_thing", "discrete", FALSE, FALSE, NULL, NULL)
  expect_equal(res$class, "my_thing")
  expect_equal(res$name, "my_thing")
  expect_equal(res$package, "my.thing")
})


test_that("can select optimisation level", {
  expect_equal(validate_compiler_options(NULL, NULL), "")
  expect_equal(validate_compiler_options("-Xf", NULL), "-Xf")
  expect_equal(validate_compiler_options(NULL, "none"), "-O0")
  expect_equal(validate_compiler_options("-Xf", "none"), "-Xf -O0")
  expect_equal(validate_compiler_options(NULL, "standard"), "-O2")
  expect_equal(validate_compiler_options(NULL, "max"), "-O3 -ffast-math")
  expect_error(validate_compiler_options(NULL, "min"),
               "Unknown 'optimisation_level' value 'min'")
})


test_that("can validate C++ standard", {
  expect_null(validate_cpp_std(NULL))
  expect_equal(validate_cpp_std("c++14"), "c++14")
  expect_equal(validate_cpp_std("C++20"), "C++20")
  expect_error(
    validate_cpp_std("verynew"),
    "'cpp_std' does not look like a valid C++ standard name (e.g., C++14)",
    fixed = TRUE)
})


test_that("validate that working directory is suitable", {
  hash <- "abc1234"
  path <- withr::local_tempfile()
  expect_equal(dust_workdir(path, hash), path)

  dir_create(file.path(path, c("src", "R")))
  expect_equal(dust_workdir(path, hash), path)
  file.create(file.path(path, c("DESCRIPTION", "NAMESPACE",
                                "src/Makevars", "src/cpp11.cpp", "src/dust.cpp",
                                "R/cpp11.R", "R/dust.R")))
  expect_equal(dust_workdir(path, hash), path)
  file.create(file.path(path, c("src/dust.o", "src/other.o", "src/dust.so")))
  expect_equal(dust_workdir(path, hash), path)

  file.create(file.path(path, "other"))
  expect_error(dust_workdir(path, hash),
               "Path '.+' does not look like a dust directory")
})


test_that("validate that working directory is suitable", {
  path <- withr::local_tempfile()
  hash <- "abc1234"
  file.create(path)
  expect_error(dust_workdir(path, hash),
               "already exists but is not a directory")
})


test_that("can work in a stable temporary directory", {
  hash <- "abc123456789"
  withr::with_envvar(
    c(DUST_WORKDIR_ROOT = "my/path"),
    expect_equal(dust_workdir(NULL, hash), "my/path/dust_abc1234"))
})


test_that("quiet responds to envvar", {
  withr::with_envvar(c(DUST_QUIET = NA_character_), {
    expect_false(dust_quiet(NULL))
    expect_false(dust_quiet(FALSE))
    expect_true(dust_quiet(TRUE))
    expect_error(dust_quiet(1), "Expected 'quiet' to be logical")
  })

  withr::with_envvar(c(DUST_QUIET = "true"), {
    expect_true(dust_quiet(NULL))
    expect_false(dust_quiet(FALSE))
  })
})


test_that("debug responds to envvar", {
  withr::with_envvar(c(DUST_DEBUG = NA_character_), {
    expect_false(dust_debug(NULL))
    expect_false(dust_debug(FALSE))
    expect_true(dust_debug(TRUE))
    expect_error(dust_debug(1), "Expected 'debug' to be logical")
  })

  withr::with_envvar(c(DUST_DEBUG = "true"), {
    expect_true(dust_debug(NULL))
    expect_false(dust_debug(FALSE))
  })
})
