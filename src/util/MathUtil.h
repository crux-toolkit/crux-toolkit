#ifndef MATHUTIL_H
#define MATHUTIL_H

#include <cmath>
#include <string>
#include <sstream>
#include <vector>

#include "boost/tuple/tuple.hpp"

class MathUtil {
 public:
  static double Round(double x, int decimals = 0);
  static bool AlmostEqual(double x, double y, int precision = 6);

  // added by Yang
  static int factorial(int n);
  static double LogNChooseK(int n, int k);
  static double LogSumExp(std::vector<double>* log_values); // en.wikipedia.org/wiki/LogSumExp

  static double MaxInArr(double* arr_values, int size);
  static double NormalizedDotProduct(double* src_values, double* tgt_values, int size, bool take_sqrt = true);
  static std::vector<double> linspace(double start, double end, int num);

  // find the closest match index from the data array to the query in O(log(n))
  static int binarySearch(const double* data_arr, int data_size, double query);
  static int binarySearch(const std::vector<double>* data_vec, double query);
  // find the closest match index from the data array to the query in O(n), as a sanity check
  static int linearSearch(const double* data_arr, int data_size, double query);
  static int linearSearch(const std::vector<double>* data_vec, double query);
  // fit linear regression
  static boost::tuple<double, double> fitLinearRegression(std::vector<double>* x_values, std::vector<double>* y_values);

  static double gammaln(double xx);

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

