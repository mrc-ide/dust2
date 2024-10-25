#include <dust2/common.hpp>
#include <numeric>

// [[dust2::class(logistic)]]
// [[dust2::time_type(continuous)]]
// [[dust2::parameter(n, type = "int", constant = TRUE, required = FALSE, rank = 0)]]
// [[dust2::parameter(r, required = TRUE, constant = FALSE, rank = 1)]]
// [[dust2::parameter(K, required = TRUE, constant = FALSE, rank = 1)]]
class logistic {
public:
  logistic() = delete;

  using mixed_time = std::false_type;
  using real_type = double;

  struct shared_state {
    size_t n;
    std::vector<real_type> r;
    std::vector<real_type> K;
  };

  struct internal_state {};
  using rng_state_type = monty::random::generator<real_type>;

  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"x", {shared.n}}, {"total", {}}};
  }

  // Based on packing_state; the *last* 'n' elements here are output.
  // We might change this later but this is quite convenient from a
  // data layout view and has the advantage that this can be computed
  // without using any reference to shared.
  static size_t size_output() {
    return 1;
  }

  static void initial(real_type time,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state) {
    for (size_t i = 0; i < shared.n; ++i) {
      state[i] = 1;
    }
  }

  static void rhs(real_type time,
                  const real_type * state,
                  const shared_state& shared,
                  internal_state& internal,
                  real_type * state_deriv) {
    for (size_t i = 0; i < shared.n; ++i) {
      state_deriv[i] = shared.r[i] * state[i] * (1 - state[i] / shared.K[i]);
    }
  }

  static void output(real_type time,
                     real_type * state,
                     const shared_state& shared,
                     internal_state& internal) {
    state[shared.n] = std::accumulate(state, state + shared.n,
                                      static_cast<real_type>(0));
  }

  static shared_state build_shared(cpp11::list pars) {
    const auto n = dust2::r::read_size(pars, "n");
    std::vector<real_type> r(n);
    std::vector<real_type> K(n);
    const auto dim = dust2::array::dimensions<1>{n};
    dust2::r::read_real_array(pars, dim, r.data(), "r", true);
    dust2::r::read_real_array(pars, dim, K.data(), "K", true);
    return shared_state{n, r, K};
  }

  static void update_shared(cpp11::list pars, shared_state& shared) {
    const auto dim = dust2::array::dimensions<1>{shared.n};
    dust2::r::read_real_array(pars, dim, shared.r.data(), "r", false);
    dust2::r::read_real_array(pars, dim, shared.K.data(), "K", false);
  }
};
