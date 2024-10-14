#pragma once

#include <cstddef>
#include <utility>
#include <vector>

namespace dust2 {

template <typename real_type>
using zero_every_type = std::vector<std::pair<real_type, std::vector<size_t>>>;

}
