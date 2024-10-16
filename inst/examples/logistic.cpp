#include <dust2/common.hpp>
#include <numeric>

// [[dust2::class(logistic)]]
// [[dust2::time_type(continuous)]]
// [[dust2::parameter(n)]]
// [[dust2::parameter(r)]]
// [[dust2::parameter(K)]]
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
    return dust2::packing{{"x", {shared.n}}};
  }

  static size_t size_output(const shared_state& shared) {
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
                     const real_type * state,
                     const shared_state& shared,
                     internal_state& internal,
                     real_type * output) {
    // We will change this to use a delay (e.g., growth over last
    // period) to give the history a good workout right away.
    output[0] = std::accumulate(state, state + shared.n, 0);
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
