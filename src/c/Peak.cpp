/**
 * AUTHOR: Kha Nguyen
 * CREATE DATE: 4/21/2011
 * 
 * \file Peak.cpp
 */

#include "Peak.h"
#include "utils.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "carp.h"
#include <algorithm>

/**
 * Construct a Peak object
 * @param intensity: intensity for the new peak
 * @param location: location for the new peak
 */
Peak::Peak(FLOAT_T intensity, FLOAT_T location) {
    this->intensity_ = intensity;
    this->intensity_rank_ = 0.0;
    this->location_ = location;
}

/**
 * Return the intensity of this Peak
 */
FLOAT_T Peak::getIntensity() {
    return this->intensity_;
}

 /**
  * Return the intensity rank of this Peak
  */
FLOAT_T Peak::getIntensityRank() {
    return this->intensity_rank_;
}

/**
 * Return the location of this Peak
 */
FLOAT_T Peak::getLocation() {
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
bool Peak::compareByIntensity(Peak other) {
    return getIntensity() > other.getIntensity();
}

/**
 * Compare the mz(location) of this Peak and another Peak
 * Return true if the other Peak is greater, false otherwise
 */
bool Peak::compareByMZ(Peak other) {
    return getLocation() < other.getLocation();
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
