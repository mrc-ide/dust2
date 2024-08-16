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
    bool require_nonscalar = false;
    for (auto& el : data_) {
      len_.push_back(tools::prod(el.second));
    }
    size_ = std::accumulate(len_.begin(), len_.end(), 0);
  }

  // We can support an incremental interface easily enough like this
  // if the initializer list version proves too annoying:
  //
  // void add(mapping_type x) {
  //   data_.push_back(x);
  //   len_.push_back(tools::prod(data_.back().second));
  //   size_ += len_.back();
  // }

  auto size() const {
    return size_;
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
