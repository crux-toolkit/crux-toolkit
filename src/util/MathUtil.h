#ifndef MATHUTIL_H
#define MATHUTIL_H

#include <cmath>
#include <string>
#include <sstream>
#include <vector>

class MathUtil {
 public:
  static double Round(double x, int decimals = 0);
  static bool AlmostEqual(double x, double y, int precision);

  template<typename T>
  static double Sum(const T& values) {
    double sum = 0;
    for (typename T::const_iterator i = values.begin(); i != values.end(); i++) {
      sum += *i;
    }
    return sum;
  }

  template<typename T>
  static double Mean(const T& values) {
    return Sum(values) / values.size();
  }

  template<typename T>
  static double Variance(const T& values, bool population = true) {
    double mean = Mean(values);
    double result = 0.0;
    for (typename T::const_iterator i = values.begin(); i != values.end(); i++) {
      result += pow(*i - mean, 2);
    }
    return population ? result / values.size() : result / (values.size() - 1);
  }

  template<typename T>
  static double StdDev(const T& values, bool population = true) {
    return sqrt(Variance(values, population));
  }

  class Combination {
   public:
    Combination(size_t n, size_t k);
    size_t N() const;
    size_t Size() const;
    const std::vector<size_t>& Values() const;
    bool HasNext() const;
    void Next();
   private:
    size_t n_;
    std::vector<size_t> values_;
  };

 private:
  MathUtil();
  ~MathUtil();
};

#endif

