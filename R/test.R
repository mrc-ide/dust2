test_trajectories <- function(time, state, order = NULL,
                         index_state = NULL, index_group = NULL,
                         select_particle = NULL, reorder = FALSE) {
  test_trajectories_(time, state, order,
                     index_state, index_group, select_particle, reorder)
}


test_snapshots <- function(time, save_snapshots, state, order = NULL,
                           index_group = NULL, select_particle = NULL,
                           reorder = FALSE) {
  test_snapshots_(time, save_snapshots, state, order,
                  index_group, select_particle, reorder)
}
