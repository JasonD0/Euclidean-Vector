#include "assignments/ev/euclidean_vector.h"

#include <algorithm>  // Look at these - they are helpful https://en.cppreference.com/w/cpp/algorithm
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>

/**
 * Default Constructor
 *  this constructor is a delegating constructor
 *  constructs a Euclidean Vector
 * @param i   number of dimensions
 */
EuclideanVector::EuclideanVector(int i) : EuclideanVector(i, 0.0) {}

/**
 * Constructor
 *  takes in the number of dimensions and the magnitude of each dimension
 *  and constructs a Euclidean Vector
 * @param i   number of dimensions
 * @param d   magnitude in each dimension
 */
EuclideanVector::EuclideanVector(int i, double d)
  : len_{(i == 0) ? 1 : i}, magnitudes_{std::make_unique<double[]>(len_)} {
  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] = d;
  }
}

/**
 * Constructor
 *  constructs a Euclidean Vector from the endpoints of a vector
 * @param begin   start of a vector
 * @param end     end of a vector
 */
EuclideanVector::EuclideanVector(std::vector<double>::const_iterator begin,
                                 std::vector<double>::const_iterator end) {
  len_ = end - begin;
  magnitudes_ = std::make_unique<double[]>(len_);

  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] = *begin++;
  }
}

/**
 * Copy Constructor
 *  copies each component of the given Euclidean Vector and constructs a new
 *  Euclidean Vector
 * @param orig    copy-from Euclidean Vector
 */
EuclideanVector::EuclideanVector(const EuclideanVector& orig)
  : len_{orig.len_}, magnitudes_{std::make_unique<double[]>(orig.len_)} {
  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] = orig.magnitudes_[j];
  }
}

/**
 * Move Constructor
 *  moves each component of the given Euclidean Vector and constructs a new
 *  Euclidean Vector
 * @param orig   moved-from Euclidean Vector
 */
EuclideanVector::EuclideanVector(EuclideanVector&& orig) noexcept
  : len_{std::move(orig.len_)}, magnitudes_{std::move(orig.magnitudes_)} {
  orig.len_ = 0;
}

/**
 * Copy Assignment
 *  copies each component of the given Euclidean Vector to the current
 *  Euclidean Vector object
 * @param v   copy-from Euclidean Vector
 * @return    the current Euclidean Vector object
 */
EuclideanVector& EuclideanVector::operator=(const EuclideanVector& v) {
  len_ = v.len_;
  magnitudes_ = std::make_unique<double[]>(len_);

  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] = v.magnitudes_[j];
  }

  return *this;
}

/**
 * Move Assignment
 *  moves each component of the given EuclideanVector to the current
 *  Euclidean Vector object
 * @param v   move-from Euclidean Vector
 * @return    the current Euclidean Vector object
 */
EuclideanVector& EuclideanVector::operator=(EuclideanVector&& v) noexcept {
  len_ = std::move(v.len_);
  magnitudes_ = std::move(v.magnitudes_);
  v.len_ = 0;

  return *this;
}

/**
 * Subscript Operator
 *  Gets the reference to the value in a given dimension of the Euclidean Vector
 * @param x   dimension of the Euclidean Vector
 * @return    reference to the value in the given dimension
 */
double& EuclideanVector::operator[](const int& x) {
  assert(x >= 0 && x < len_);
  return this->magnitudes_[x];
}

/**
 * Subscript Operator
 *  Gets the value in a given dimension of the Euclidean Vector
 * @param x   dimension of the Euclidean Vector
 * @return    value at the given dimension
 */
double EuclideanVector::operator[](const int& x) const {
  assert(x >= 0 && x < len_);
  return this->magnitudes_[x];
}

/**
 * Addition of two Euclidean Vectors
 *  point-wise addition onto the current Euclidean Vector object
 * @param v   a Euclidean Vector
 * @return    the current Euclidean Vector object containing each point-wise
 *            addition value
 */
EuclideanVector& EuclideanVector::operator+=(const EuclideanVector& v) {
  if (len_ != v.len_) {
    throw EuclideanVectorError("Dimensions of LHS(X) and RHS(Y) do not match");
  }

  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] += v.magnitudes_[j];
  }

  return *this;
}

/**
 * Subtraction of two Euclidean Vectors
 *  point-wise subtraction onto the current Euclidean Vector object
 * @param v   a Euclidean Vector
 * @return    the current Euclidean Vector object containing each point-wise
 *            subtraction value
 */
EuclideanVector& EuclideanVector::operator-=(const EuclideanVector& v) {
  if (len_ != v.len_) {
    throw EuclideanVectorError("Dimensions of LHS(X) and RHS(Y) do not match");
  }

  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] -= v.magnitudes_[j];
  }

  return *this;
}

/**
 * Scalar multiplication of the current Euclidean Vector object
 *  scales each element of the current Euclidean Vector object by the scalar
 * @param x   the scalar
 * @return    the current Euclidean Vector object containing each scaled value
 */
EuclideanVector& EuclideanVector::operator*=(const int& x) {
  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] *= x;
  }

  return *this;
}

/**
 * Scalar division of the current Euclidean Vector object
 *  divides each element of the current Euclidean Vector object by the scalar
 * @param x   the scalar
 * @return    the current EuclideanVector object containg each scaled value
 */
EuclideanVector& EuclideanVector::operator/=(const int& x) {
  if (x == 0) {
    throw EuclideanVectorError("Invalid vector division by 0");
  }

  for (int j = 0; j < len_; ++j) {
    magnitudes_[j] /= x;
  }

  return *this;
}

/**
 * Vector Type Conversion
 *  convert the current Euclidean Vector object into a vector
 * @return    vector with the same elements as the current Euclidean Vector object
 */
EuclideanVector::operator std::vector<double>() const {
  std::vector<double> v;

  for (int j = 0; j < len_; ++j) {
    v.push_back(magnitudes_[j]);
  }

  return v;
}

/**
 * List Type Conversion
 *  convert the current Euclidean Vector object into a list
 * @return    list with the same elements as the current Euclidean Vector object
 */
EuclideanVector::operator std::list<double>() const {
  std::list<double> l;

  for (int j = 0; j < len_; ++j) {
    l.push_back(magnitudes_[j]);
  }

  return l;
}

/**
 * Gets the reference of the magnitude in the given dimension
 * @param x   dimension
 * @return    reference to the value in the given dimension
 */
double& EuclideanVector::at(const int& x) {
  if (x < 0 || x >= len_) {
    throw EuclideanVectorError("Index X is not valid for this EuclideanVector object");
  }

  return magnitudes_[x];
}

/**
 * Gets the value of the magnitude in the given dimension
 * @param x   dimension
 * @return    value of the magnitude in the given dimension
 */
double EuclideanVector::at(const int& x) const {
  if (x < 0 || x >= len_) {
    throw EuclideanVectorError("Index X is not valid for this EuclideanVector object");
  }

  return magnitudes_[x];
}

/**
 * Gets the number of dimensions of the current Euclidean Vector object
 * @return    number of dimensions
 */
int EuclideanVector::GetNumDimensions() const {
  return len_;
}

/**
 * Gets the Euclidean norm of the current Euclidean Vector object
 *  calculates the square of the magnitudes in each dimension
 * @return    the Euclidean norm
 */
double EuclideanVector::GetEuclideanNorm() const {
  if (this->GetNumDimensions() == 0) {
    throw EuclideanVectorError("EuclideanVector with no dimensions does not have a norm");
  }

  double sum = 0;
  for (int j = 0; j < len_; ++j) {
    sum += std::pow(magnitudes_[j], 2);
  }

  return std::sqrt(sum);
}

/**
 * Creates a Euclidean Vector that is the unit vector of the current Euclidean
 * Vector object
 *  the magnitude for each dimension of the unit vector is the original vector's
 *  magnitude divided by the Euclidean Norm
 * @return    the unit vector
 */
EuclideanVector EuclideanVector::CreateUnitVector() const {
  if (this->GetNumDimensions() == 0) {
    throw EuclideanVectorError("EuclideanVector with no dimensions does not have a unit vector");
  }

  if (this->GetEuclideanNorm() == 0) {
    throw EuclideanVectorError(
        "EuclideanVector with euclidean normal of 0 does not have a unit vector");
  }

  std::vector<double> v;

  double norm = this->GetEuclideanNorm();
  for (int j = 0; j < len_; ++j) {
    v.push_back(magnitudes_[j] / norm);
  }

  return EuclideanVector{v.begin(), v.end()};
}

/**
 * Equals Operator (==)
 *  compares two Euclidean Vectors
 * @param v1    a Euclidean Vector
 * @param v2    a Euclidean Vector
 * @return      if equal then true
 *              otherwise false
 */
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

/**
 * Not Equal Operator (!=)
 *  compares two Euclidean Vectors
 * @param v1    a Euclidean Vector
 * @param v2    a Euclidean Vector
 * @return      if not equal then true
 *              otherwise false
 */
bool operator!=(const EuclideanVector& v1, const EuclideanVector& v2) {
  return !(v1 == v2);
}

/**
 * Addition of two Euclidean Vectors
 *  point-wise addition of two Euclidean Vectors
 * @param v1    a Euclidean Vector
 * @param v2    a Euclidean Vector
 * @return      a Euclidean Vector containing each point-wise addition value
 */
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

/**
 * Subtraction of two Euclidean Vectors
 *  point-wise subtraction of two Euclidean Vectors
 * @param v1    a Euclidean Vector
 * @param v2    a Euclidean Vector
 * @return      a Euclidean Vector containing each point-wise subtraction value
 */
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

/**
 * Dot Product of two Euclidean Vectors
 *  calculates the summation of the point-wise multiplication of the Euclidean
 *  Vectors
 * @param v1    a Euclidean Vector
 * @param v2    a Euclidean Vector
 * @return      the dot product
 */
double operator*(const EuclideanVector& v1, const EuclideanVector& v2) {
  if (v1.len_ != v2.len_) {
    throw EuclideanVectorError("Dimensions of LHS(X) and RHS(Y) do not match");
  }

  double sum{0.0};

  for (int j = 0; j < v1.len_; ++j) {
    sum += v1.magnitudes_[j] * v2.magnitudes_[j];
  }

  return sum;
}

/**
 * Scalar Multiplication
 *  scales each element of the Euclidean Vector by the scalar
 * @param ev    a Euclidean Vector
 * @param i     the scalar
 * @return      a Euclidean Vector containing the scaled values
 */
EuclideanVector operator*(const EuclideanVector& ev, const int& i) {
  std::vector<double> v;

  for (int j = 0; j < ev.len_; ++j) {
    v.push_back(ev.magnitudes_[j] * i);
  }

  return EuclideanVector{v.begin(), v.end()};
}

/**
 * Scalar Multiplication
 *  scales each element of the Euclidean Vector by the scalar
 * @param i     the scalar
 * @param ev    a Euclidean Vector
 * @return      a Euclidean Vector containing the scaled values
 */
EuclideanVector operator*(const int& i, const EuclideanVector& ev) {
  return ev * i;
}

/**
 * Scalar Division
 *  divides each element of the Euclidean Vector by the scalar
 * @param ev    a Euclidean Vector
 * @param i     the scalar
 * @return      a Euclidean Vector containing the scaled values
 */
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

/**
 * Output Stream
 *  Prints the magnitude in each dimension of the Euclidean Vector
 * @param os    the output stream
 * @param v     a Euclidean Vector
 * @return      the output stream with the Euclidean Vector
 */
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
