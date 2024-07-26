#pragma once

// #include <algorithm>
// #include <numeric>
#include <sstream>
#include <vector>

namespace dust2 {
namespace utils {

class errors {
public:
  errors() : errors(0) {
  }
  errors(size_t len) :
    count(0), err(len), seen(len) {
  }

  void reset() {
    count = 0;
    std::fill(seen.begin(), seen.end(), false);
    std::fill(err.begin(), err.end(), "");
  }

  bool unresolved() const {
    return count > 0;
  }

  template <typename T>
  void capture(const T& e, size_t i) {
    err[i] = e.what();
    seen[i] = true;
  }

  void report(bool clear = false, const char *title = "particles",
              size_t n_max = 4) {
    count = std::accumulate(std::begin(seen), std::end(seen), 0);
    if (count == 0) {
      return;
    }

    std::stringstream msg;
    msg << count << " " << title << " reported errors.";

    const size_t n_report = std::min(n_max, count);
    for (size_t i = 0, j = 0; i < seen.size() && j < n_report; ++i) {
      if (seen[i]) {
        msg << std::endl << "  - " << i + 1 << ": " << err[i];
        ++j;
      }
    }
    if (n_report < count) {
      msg << std::endl << "  - (and " << (count - n_report) << " more)";
    }

    if (clear) {
      reset();
    }

    throw std::runtime_error(msg.str());
  }

private:
  size_t count;
  std::vector<std::string> err;
  std::vector<bool> seen;
};

}
}
