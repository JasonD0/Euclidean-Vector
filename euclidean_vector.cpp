#include "assignments/ev/euclidean_vector.h"

#include <algorithm>  // Look at these - they are helpful https://en.cppreference.com/w/cpp/algorithm
#include <cmath>
#include <iostream>
#include <cassert>

EuclideanVector::EuclideanVector(int i) : EuclideanVector(i, 0.0) {}

EuclideanVector::EuclideanVector(int i, double d)
    : len_{(i == 0) ? 1 : i},
      magnitudes_{std::make_unique<double[]>(len_)} {
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
    : len_{std::move(orig.len_)},
      magnitudes_{std::move(orig.magnitudes_)} {
  orig.len_ = 0;
  orig.magnitudes_ = nullptr;
}


EuclideanVector& EuclideanVector::operator=(const EuclideanVector& v) {
  len_ = v.len_;
  magnitudes_ = std::make_unique<double[]>(len_);

  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] = v.magnitudes_[j];
  }

  return *this;
}

EuclideanVector& EuclideanVector::operator=(EuclideanVector&& v) noexcept {
  len_ = v.len_;
  magnitudes_ = std::move(v.magnitudes_);

  return *this;
}

double& EuclideanVector::operator[](int x) {
  assert(x >= 0 && x < len_);
  return this->magnitudes_[x];
}

double EuclideanVector::operator[](int x) const {
  assert(x >= 0 && x < len_);
  return this->magnitudes_[x];
}

EuclideanVector& EuclideanVector::operator+=(const EuclideanVector& v) {
  if (len_ != v.len_) {
    throw EuclideanVectorError("Dimensions of LHS(X) and RHS(Y) do not match");
  }

  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] += v.magnitudes_[j];
  }

  return *this;
}

EuclideanVector& EuclideanVector::operator-=(const EuclideanVector& v) {
  if (len_ != v.len_) {
    throw EuclideanVectorError("Dimensions of LHS(X) and RHS(Y) do not match");
  }

  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] -= v.magnitudes_[j];
  }

  return *this;
}

EuclideanVector& EuclideanVector::operator*=(int x) {
  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] *= x;
  }

  return *this;
}

EuclideanVector& EuclideanVector::operator/=(int x) {
  if (x == 0) {
    throw EuclideanVectorError("Invalid vector division by 0");
  }

  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] /= x;
  }

  return *this;
}

EuclideanVector::operator std::vector<double>() {
  std::vector<double> v;

  for (int j = 0; j < len_; ++j) {
    v.push_back(magnitudes_[j]);
  }

  return v;
}

EuclideanVector::operator std::list<double>() {
  std::list<double> l;

  for (int j = 0; j < len_; ++j) {
    l.push_back(magnitudes_[j]);
  }

  return l;
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

EuclideanVector EuclideanVector::CreateUnitVector() {
  if (this->GetNumDimensions() == 0) {
    throw EuclideanVectorError("EuclideanVector with no dimensions does not have a unit vector");
  }

  if (this->GetEuclideanNorm() == 0) {
    throw EuclideanVectorError("EuclideanVector with euclidean normal of 0 does not have a unit vector");
  }

  std::vector<double> v;

  double norm = this->GetEuclideanNorm();
  for (int j = 0; j < len_; ++j) {
    v.push_back(magnitudes_[j]/norm);
  }

  return EuclideanVector{v.begin(), v.end()};
}


bool operator==(const EuclideanVector& v1, const EuclideanVector& v2) {
  if (v1.len_ != v2.len_) {
    return false;
  }

  for (int j = 0; j < v1.len_; ++j) {
    if (v1.magnitudes_[j] != v2.magnitudes_[j]) {
      return false;
    }
  }

  return true;
}

bool operator!=(const EuclideanVector& v1, const EuclideanVector& v2) {
  return !(v1 == v2);
}

EuclideanVector operator+(const EuclideanVector& v1, const EuclideanVector& v2) {
  if (v1.len_ != v2.len_) {
    throw EuclideanVectorError("Dimensions of LHS(X) and RHS(Y) do not match");
  }

  std::vector<double> v;

  for (int j = 0; j < v1.len_; ++j) {
    v.push_back(v1.magnitudes_[j] + v2.magnitudes_[j]);
  }

  return EuclideanVector{v.begin(), v.end()};
}

EuclideanVector operator-(const EuclideanVector& v1, const EuclideanVector& v2) {
  if (v1.len_ != v2.len_) {
    throw EuclideanVectorError("Dimensions of LHS(X) and RHS(Y) do not match");
  }

  std::vector<double> v;

  for (int j = 0; j < v1.len_; ++j) {
    v.push_back(v1.magnitudes_[j] - v2.magnitudes_[j]);
  }

  return EuclideanVector{v.begin(), v.end()};
}

double operator*(const EuclideanVector& v1, const EuclideanVector& v2) {
  if (v1.len_ != v2.len_) {
    throw EuclideanVectorError("Dimensions of LHS(X) and RHS(Y) do not match");
  }

  double sum {0.0};

  for (int j = 0; j < v1.len_; ++j) {
    sum += v1.magnitudes_[j] * v2.magnitudes_[j];
  }

  return sum;
}

EuclideanVector operator*(const EuclideanVector& ev, const int& i) {
  std::vector<double> v;

  for (int j = 0; j < ev.len_; ++j) {
    v.push_back(ev.magnitudes_[j] * i);
  }

  return EuclideanVector{v.begin(), v.end()};
}

EuclideanVector operator*(const int& i, const EuclideanVector& ev) {
  /*for (int j = 0; j < len_; ++j) {
    this->magnitudes_[j] = v1.magnitudes_[j] * i;
  }
*/
  return ev*i;
}

EuclideanVector operator/(const EuclideanVector& ev, const int& i) {
  if (i == 0) {
    throw EuclideanVectorError("Invalid vector division by 0");
  }

  std::vector<double> v;

  for (int j = 0; j < ev.len_; ++j) {
    v.push_back(ev.magnitudes_[j] / i);
  }

  return EuclideanVector{v.begin(), v.end()};
}

std::ostream& operator<<(std::ostream& os, const EuclideanVector& v) {
  int len = v.len_;

  if (len == 0) {
    os << "[]";
    return os;
  }

  os << "[";
  for (int j = 0; j < len - 1; ++j) {
    os << v.magnitudes_[j] << " ";
  }
  os << v.magnitudes_[len - 1];
  os << "]";

  return os;
}
