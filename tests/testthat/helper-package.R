create_test_package <- function(path, name = "pkg",
                                examples = c("walk.cpp", "sir.cpp")) {
  dir.create(path, FALSE, TRUE)
  dir.create(file.path(path, "inst/dust"), FALSE, TRUE)
  dir.create(file.path(path, "R"), FALSE, TRUE)
  dir.create(file.path(path, "src"), FALSE, TRUE)

  ns <- sprintf('useDynLib("%s", .registration = TRUE)', name)
  desc <- c(
    sprintf("Package: %s", name),
    "Imports: dust2",
    "LinkingTo: cpp11, dust2, mcstate2",
    "Version: 0.0.1",
    "Authors@R: c(person('A', 'Person', role = c('aut', 'cre'),",
    "                    email = 'person@example.com'))",
    "Title: An Example Using Dust",
    "Description: Provides an example of dust within a packge. Not intended to",
    "    be used for anything really.",
    "License: CC0")

  writeLines(desc, file.path(path, "DESCRIPTION"))
  writeLines(ns, file.path(path, "NAMESPACE"))
  file.copy(dust2_file(file.path("examples", examples)),
            file.path(path, "inst/dust"))

  invisible(path)
}
