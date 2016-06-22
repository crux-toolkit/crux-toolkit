#include "MathUtil.h"

#include <cmath>

using namespace std;

double MathUtil::Round(double x, int decimals) {
  double shift = pow(10, decimals);
  return (x >= 0)
    ? floor(x * shift + 0.5) / shift
    : ceil(x * shift - 0.5) / shift;
}

bool MathUtil::AlmostEqual(double x, double y, int precision) {
  return abs(x - y) < 5*pow(10, -(precision + 1));
}

MathUtil::MathUtil() {}
MathUtil::~MathUtil() {}

