test_history <- function(time, state, order = NULL, index_group = NULL,
                         index_particle = NULL, reorder = FALSE) {
  ## TODO: move index_particle ahead of index_group everywhere, once
  ## we're stable.
  test_history_(time, state, order, index_group, index_particle, reorder)
}
