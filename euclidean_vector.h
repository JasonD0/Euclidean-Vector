#ifndef ASSIGNMENTS_EV_EUCLIDEAN_VECTOR_H_
#define ASSIGNMENTS_EV_EUCLIDEAN_VECTOR_H_

#include <exception>
#include <string>
#include <memory>
#include <vector>
#include <list>

class EuclideanVectorError : public std::exception {
 public:
  explicit EuclideanVectorError(const std::string& what) : what_(what) {}
  const char* what() const noexcept{ return what_.c_str(); }
 private:
  std::string what_;
};

class EuclideanVector {
 public:
  explicit EuclideanVector(int i);
  explicit EuclideanVector(int i, double d);
  explicit EuclideanVector(std::vector<double>::const_iterator begin,
                           std::vector<double>::const_iterator end);
  EuclideanVector(const EuclideanVector& orig);
  EuclideanVector(EuclideanVector&& orig) noexcept;
  ~EuclideanVector() = default;

  EuclideanVector& operator=(const EuclideanVector& v);
  EuclideanVector& operator=(EuclideanVector&& v) noexcept;
  double& operator[](int x);
  double operator[](int x) const;
  EuclideanVector& operator+=(const EuclideanVector& v);
  EuclideanVector& operator-=(const EuclideanVector& v);
  EuclideanVector& operator*=(int x);
  EuclideanVector& operator/=(int x);
  explicit operator std::vector<double>();
  explicit operator std::list<double>();

  double at(int x);
  int GetNumDimensions();
  double GetEuclideanNorm();
  EuclideanVector CreateUnitVector();

  friend bool operator==(const EuclideanVector& v1, const EuclideanVector& v2);
  friend bool operator!=(const EuclideanVector& v1, const EuclideanVector& v2);
  friend EuclideanVector operator+(const EuclideanVector& v1, const EuclideanVector& v2);
  friend EuclideanVector operator-(const EuclideanVector& v1, const EuclideanVector& v2);
  friend double operator*(const EuclideanVector& v1, const EuclideanVector& v2);
  friend EuclideanVector operator*(const EuclideanVector& ev, const int& i);
  friend EuclideanVector operator*(const int& i, const EuclideanVector& ev);
  friend EuclideanVector operator/(const EuclideanVector& ev, const int& i);
  friend std::ostream& operator<<(std::ostream& os, const EuclideanVector& v);

 private:
  int len_;
  std::unique_ptr<double[]> magnitudes_;
};

#endif
