#pragma once

#include <numeric>
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

  void add(mapping_type x) {
    data_.push_back(x);
    len_.push_back(tools::prod(data_.back().second));
    size_ += len_.back();
  }

  auto names() const {
    std::vector<std::string> ret;
    ret.reserve(data_.size());
    for (const auto& el : data_) {
      ret.push_back(el.first);
    }
    return ret;
  }

  auto size() const {
    return size_;
  }

private:
  std::vector<mapping_type> data_;
  std::vector<size_t> len_;
  size_t size_;
};

}
