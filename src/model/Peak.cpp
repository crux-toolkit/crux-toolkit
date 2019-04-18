/**
 * AUTHOR: Kha Nguyen
 * CREATE DATE: 4/21/2011
 * 
 * \file Peak.cpp
 */

#include "Peak.h"
#include "util/utils.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "io/carp.h"
#include <algorithm>

Peak::Peak():
  intensity_(0.0), intensity_rank_(0.0), location_(0.0) {
}

/**
 * Construct a Peak object
 * @param intensity: intensity for the new peak
 * @param location: location for the new peak
 */
Peak::Peak(FLOAT_T intensity, FLOAT_T location):
  intensity_(intensity), intensity_rank_(0.0), location_(location) {
}

Peak::Peak(const Peak& other):
  intensity_(other.intensity_), intensity_rank_(other.intensity_rank_), location_(other.location_) {
}

void swap(Peak& x, Peak& y) {
  using std::swap;
  swap(x.intensity_, y.intensity_);
  swap(x.intensity_rank_, y.intensity_rank_);
  swap(x.location_, y.location_);
}

Peak& Peak::operator=(Peak rhs) {
  swap(*this, rhs);
  return *this;
}

/**
 * Return the intensity of this Peak
 */
FLOAT_T Peak::getIntensity() const {
    return this->intensity_;
}

 /**
  * Return the intensity rank of this Peak
  */
FLOAT_T Peak::getIntensityRank() const {
    return this->intensity_rank_;
}

/**
 * Return the location of this Peak
 */
FLOAT_T Peak::getLocation() const {
    return this->location_;
}

/**
 * Set the intensity of this Peak
 */
void Peak::setIntensity(FLOAT_T intensity) {
    this->intensity_ = intensity;
}

/**
 * Set the intensity rank of this peak
 */
void Peak::setIntensityRank(FLOAT_T intensity_rank) {
    this->intensity_rank_ = intensity_rank;
}

/**
 * Set the location of this Peak
 */
void Peak::setLocation(FLOAT_T location) {
    this->location_ = location;
}

/**
 * Print the intensity and location of this peak to stdout
 */
void Peak::print() {
    printf("%.1f %.1f\n",
            this->location_,
            this->intensity_);
}

/**
 * Compare the intensity of this Peak and another Peak
 * Return true if this Peak is greater, false otherwise
 */
bool Peak::compareByIntensity(const Peak& x, const Peak& y) {
    return x.intensity_ > y.intensity_;
}

/**
 * Compare the mz(location) of this Peak and another Peak
 * Return true if the other Peak is greater, false otherwise
 */
bool Peak::compareByMZ(const Peak& x, const Peak& y) {
    return x.location_ < y.location_;
}


std::vector<Peak*> allocate_peak_vector(unsigned int num_peaks) {
  std::vector<Peak*> ans;

  for (unsigned int idx = 0; idx < num_peaks; idx++) {
    ans.push_back(new Peak(0,0));
  }

  return ans;
}

void free_peak_vector(std::vector<Peak*> &peaks) {
  for (unsigned int idx = 0; idx < peaks.size(); idx++) {
    delete peaks[idx];
  }
  peaks.clear();
}

bool compare_peaks_by_intensity(Peak* peak_one, Peak* peak_two) {
  return (peak_one->getIntensity() > peak_two->getIntensity());
}

bool compare_peaks_by_mz(Peak* peak_one, Peak* peak_two) {
  return (peak_one->getLocation() < peak_two->getLocation());
}

void sort_peaks(std::vector<Peak*> &peak_array, PEAK_SORT_TYPE_T sort_type) {
  if (sort_type == _PEAK_INTENSITY) {
    sort(peak_array.begin(), peak_array.end(), compare_peaks_by_intensity);
  } else if (sort_type == _PEAK_LOCATION) {
    sort(peak_array.begin(), peak_array.end(), compare_peaks_by_mz);
  } else {
    carp(CARP_ERROR, "no matching peak sort type");
  }
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
