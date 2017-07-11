#include "MathUtil.h"

#include <cmath>
#include <stdexcept>

using namespace std;

double MathUtil::Round(double x, int decimals) {
  double shift = pow(10.0, decimals);
  return (x >= 0)
    ? floor(x * shift + 0.5) / shift
    : ceil(x * shift - 0.5) / shift;
}

bool MathUtil::AlmostEqual(double x, double y, int precision) {
  return abs(x - y) < 5*pow(10.0, -(precision + 1));
}

MathUtil::Combination::Combination(size_t n, size_t k)
  : n_(n) {
  if (n < 0) {
    throw runtime_error("Combination constructed with negative parameter");
  }
  for (size_t i = 0; i < k; i++) {
    values_.push_back(i);
  }
}

size_t MathUtil::Combination::N() const {
  return n_;
}

size_t MathUtil::Combination::Size() const {
  return values_.size();
}

const vector<size_t>& MathUtil::Combination::Values() const {
  return values_;
}

bool MathUtil::Combination::HasNext() const {
  return !values_.empty() && values_.front() != n_ - values_.size();
}

void MathUtil::Combination::Next() {
  if (!HasNext()) {
    throw runtime_error("Combination does not have next");
  }
  int k = (int)Size();
  int i;
  for (i = k - 1; i > 0 && values_[i] == n_ - k + i; --i);
  ++values_[i];
  for (int j = i; j < k - 1; ++j) {
    values_[j + 1] = values_[j] + 1;
  }
}

MathUtil::MathUtil() {}
MathUtil::~MathUtil() {}

