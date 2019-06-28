#include "assignments/ev/euclidean_vector.h"

#include <algorithm>  // Look at these - they are helpful https://en.cppreference.com/w/cpp/algorithm

EuclideanVector::EuclideanVector(int i) {
  EuclideanVector(i, 0.0);
}

EuclideanVector::EuclideanVector(int i, double d)
    : len_{(i == 0) ? 1 : i} {
  magnitudes_ = std::make_unique<double[]>(len_);

  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] = d;
  }
}

EuclideanVector::EuclideanVector(std::vector<double>::const_iterator begin,
                                 std::vector<double>::const_iterator end) {
  len_ = end - begin;

  magnitudes_ = std::make_unique<double[]>(len_);

  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] = *begin++;
  }
}

EuclideanVector::EuclideanVector(const EuclideanVector& orig)
    : len_{orig.len_} {
  //magnitudes_  orig.magnitudes_;
}

EuclideanVector::EuclideanVector(EuclideanVector&& orig) noexcept
    : len_{orig.len_} {
  orig.len_ = 0;
}
