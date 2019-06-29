#include "assignments/ev/euclidean_vector.h"

#include <algorithm>  // Look at these - they are helpful https://en.cppreference.com/w/cpp/algorithm
#include <cmath>
#include <iostream>

EuclideanVector::EuclideanVector(int i) : EuclideanVector(i, 0.0) {}

EuclideanVector::EuclideanVector(int i, double d) : len_{(i == 0) ? 1 : i} {
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
    : len_{orig.len_},
      magnitudes_{std::make_unique<double[]>(orig.len_)} {
  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] = orig.magnitudes_[j];
  }
}

EuclideanVector::EuclideanVector(EuclideanVector&& orig) noexcept
    : len_{orig.len_}, magnitudes_{std::move(orig.magnitudes_)} {
  orig.len_ = 0;
  orig.magnitudes_.reset();
}





double EuclideanVector::at(int x) {
  if (x < 0 || x >= len_) {
    throw EuclideanVectorError("Index X is not valid for this EuclideanVector object");
  }

  return magnitudes_[x];
}

int EuclideanVector::GetNumDimensions() {
  return len_;
}

double EuclideanVector::GetEuclideanNorm() {
  if (this->GetNumDimensions() == 0) {
    throw EuclideanVectorError("EuclideanVector with no dimensions does not have a norm");
  }

  double sum = 0;
  for (int j = 0; j < len_; ++j) {
    sum += std::pow(magnitudes_[j], 2);
  }

  return std::sqrt(sum);
}
/*
EuclideanVector EuclideanVector::CreateUnitVector() {
  if (this->GetNumDimensions() == 0) {
    throw EuclideanVectorError("EuclideanVector with no dimensions does not have a unit vector");
  }

  EuclideanVector ev(len_);

  double norm = this->GetEuclideanNorm();
  for (int j = 0; j < len_; ++j) {
    ev.magnitudes_[j] = magnitudes_[j]/norm;
  }

  return ev;
}*/

std::ostream& operator<<(std::ostream& os, const EuclideanVector& v) {
  int len = v.len_;

  os << "[";
  for (int j = 0; j < len - 1; ++j) {
    os << v.magnitudes_[j] << " ";
  }
  os << v.magnitudes_[len - 1];
  os << "]";

  return os;
}

