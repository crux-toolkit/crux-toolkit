#include "MathUtil.h"

#include <cmath>
#include <stdexcept>
#include <algorithm>

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

// added by Yang
int MathUtil::factorial(int n) {
  int product = 1;
  for (int i = 1; i <= n; i++) { product *= i; }
  return product;
}

double MathUtil::LogNChooseK(int n, int k) {
	double result = 0;
	for (int i = 1; i <= k; i++) { result += log(n-i+1) - log(i); }
	return result;
}

double MathUtil::LogSumExp(std::vector<double>* log_values) {
	if (log_values->size() <= 0) { return 0.0; }
	else {
		double max_value = *max_element(log_values->begin(), log_values->end());
		double sum = 0;
		for (int i = 0; i < log_values->size() ; i++){ sum += exp(log_values->at(i) - max_value); }
		return log(sum) + max_value;
	}
}

double MathUtil::MaxInArr(double* arr_values, int size) {
	if (arr_values == NULL || size <= 0) { return 0; }

	double max_value = arr_values[0];
	for (int i=0; i<size; ++i) {
		if (arr_values[i] > max_value) { max_value = arr_values[i]; }
	}
	return max_value;
}

double MathUtil::NormalizedDotProduct(double* src_values, double* tgt_values, int size, bool take_sqrt) {
	if (src_values == NULL || tgt_values == NULL || size <= 0) { return 0; }

	double prod_sum = 0, src_sum = 0, tgt_sum = 0;
	for (int i=0; i<size; ++i) {
		if (take_sqrt) {
			prod_sum += sqrt(src_values[i]) * sqrt(tgt_values[i]);
			src_sum += src_values[i];
			tgt_sum += tgt_values[i];
		} else {
			prod_sum += src_values[i] * tgt_values[i];
			src_sum += src_values[i] * src_values[i];
			tgt_sum += tgt_values[i] * tgt_values[i];
		}
	}

	if (AlmostEqual(prod_sum, 0)) { return 0; }
	else { return prod_sum / sqrt(src_sum * tgt_sum); }
}

std::vector<double> MathUtil::linspace(double start, double end, int num) {
	std::vector<double> linspaced;

	if (num < 2) { throw runtime_error("linspace should contain at least two points"); }
	double delta = (end - start) / (num - 1);
	for(int i=0; i < num-1; ++i) { linspaced.push_back(start + delta * i); }
	linspaced.push_back(end);

	return linspaced;
}

int MathUtil::binarySearch(const double* data_arr, int data_size, double query) {
	if (data_arr == NULL || data_size <= 0) { return 0; }

	if (query < data_arr[0]) { return data_arr[0]; }
	if (query > data_arr[data_size-1]) { return data_arr[data_size-1]; }

	int lo = 0, hi = data_size-1;
	while (lo <= hi) {
		int mid = int((hi + lo)/2);
		if (query < data_arr[mid]) { hi = mid - 1; }
		else if (query > data_arr[mid]) { lo = mid + 1; }
		else { return data_arr[mid]; }
	}

	return fabs(data_arr[lo] - query) < fabs(query - data_arr[hi]) ? lo : hi;
}

int MathUtil::linearSearch(const double* data_arr, int data_size, double query) {
	if (data_arr == NULL || data_size <= 0) { return 0; }

	int min_idx = 0;
	double min_diff = fabs(data_arr[0] - query);

	for (int idx=0; idx<data_size; ++idx) {
		double curr_diff = fabs(data_arr[idx] - query);
		if (curr_diff <= min_diff) {
			min_idx = idx;
			min_diff = curr_diff;
		}
	}
	return min_idx;
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

