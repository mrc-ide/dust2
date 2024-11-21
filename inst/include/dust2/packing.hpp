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
    : data_(data), size_(0) {
    len_.reserve(data_.size());
    offset_.reserve(data_.size());
    for (auto& el : data_) {
      const auto len_i = tools::prod(el.second);
      len_.push_back(len_i);
      offset_.push_back(size_);
      size_ += len_i;
    }
  }

  packing() : size_(0) {
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

  template <typename Iter>
  void copy_offset(Iter it) {
    std::copy(offset_.begin(), offset_.end(), it);
  }

  bool operator!=(const packing& other) const {
    return data_ != other.data();
  }

private:
  std::vector<mapping_type> data_;
  std::vector<size_t> len_;
  std::vector<size_t> offset_;
  size_t size_;
};

}
