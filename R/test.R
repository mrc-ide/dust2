test_trajectories <- function(time, state, order = NULL,
                              index_state = NULL, index_group = NULL,
                              select_particle = NULL, times_snapshot = NULL,
                              save_state = TRUE, reorder = FALSE) {
  test_trajectories_(time, state, order,
                     index_state, index_group, select_particle,
                     times_snapshot, save_state, reorder)
}
