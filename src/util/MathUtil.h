#ifndef MATHUTIL_H
#define MATHUTIL_H

class MathUtil {
 public:
  static double Round(double x, int decimals = 0);
  static bool AlmostEqual(double x, double y, int precision);

 private:
  MathUtil();
  ~MathUtil();
};

#endif

