#pragma once

#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <dust2/tools.hpp>

namespace dust2 {

class packing {
  using mapping_type = std::pair<std::string, std::vector<size_t>>;
public:
  packing(std::initializer_list<mapping_type> data)
    : data_(data) {
    len_.reserve(data_.size());
    for (auto& el : data_) {
      len_.push_back(tools::prod(el.second));
    }
    size_ = std::accumulate(len_.begin(), len_.end(), 0);
  }

  auto size() const {
    return size_;
  }

  auto& len() const {
    return len_;
  }

  auto& data() const {
    return data_;
  }

  bool operator!=(const packing& other) const {
    return data_ != other.data();
  }

private:
  std::vector<mapping_type> data_;
  std::vector<size_t> len_;
  size_t size_;
};

}
