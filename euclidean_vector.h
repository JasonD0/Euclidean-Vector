// TODO(you): Include header guards

#include <exception>
#include <string>
#include <memory>
#include <vector>

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
  explicit EuclideanVector(const EuclideanVector& orig);
  explicit EuclideanVector(EuclideanVector&& orig) noexcept;
  ~EuclideanVector() = default;

  EuclideanVector& operator=(const EuclideanVector& v);
  EuclideanVector& operator=(EuclideanVector&& v) noexcept;

  double at(int x);
  int GetNumDimensions();
  double GetEuclideanNorm();
  //EuclideanVector CreateUnitVector();

  friend std::ostream& operator<<(std::ostream& os, const EuclideanVector& v);

 private:
  int len_;
  std::unique_ptr<double[]> magnitudes_;
};
