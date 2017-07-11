#ifndef MATHUTIL_H
#define MATHUTIL_H

#include <string>
#include <sstream>
#include <vector>

class MathUtil {
 public:
  static double Round(double x, int decimals = 0);
  static bool AlmostEqual(double x, double y, int precision);

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

