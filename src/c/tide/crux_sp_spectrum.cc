// Tahmina Baker
//
// This file contains implementations for the classes defined in 
// crux_sp_spectrum.h. Please see the header file for details.

#include <math.h>
#include "crux_sp_spectrum.h"

SpSpectrum::SpSpectrum(const Spectrum& spectrum, int charge, double max_mz) 
  : beta_(0.075), max_intensity_(0.0), last_idx_(0) {
  max_mz_ = MassConstants::mass2bin(max_mz);  
  intensity_array_ = new double[IntensityArraySize()];
  memset(intensity_array_, 0, sizeof(double)*IntensityArraySize());

  PreprocessSpectrum(spectrum, charge);
}

SpSpectrum::~SpSpectrum() {
  delete [] intensity_array_;
}

void SpSpectrum::PreprocessSpectrum(const Spectrum& spectrum, int charge) {
  double precursor_mz = spectrum.PrecursorMZ();
  double experimental_mass_cut_off = precursor_mz*charge + 50;
  double max_intensity = 0;

  for (int i = 0; i < spectrum.Size(); ++i) {
    double peak_location = spectrum.M_Z(i);

    // skip all peaks larger than experimental mass
    if(peak_location > experimental_mass_cut_off)
      continue;

    // skip all peaks within precursor ion mz +/- 15
    if((peak_location < precursor_mz + 15) &&  
       (peak_location > precursor_mz - 15))
      continue;
    
    // map peak location to bin
    int bin = MassConstants::mass2bin(peak_location);
    double intensity = sqrt(spectrum.Intensity(i));

    // set intensity in array with correct mz, only if max peak in the bin
    if(intensity_array_[bin] < intensity) {
      intensity_array_[bin] = intensity;
      
      // check if max_intensity
      if(intensity > max_intensity)
        max_intensity = intensity;
    }
    
    // set last idx to the largest added peak mz value
    if(last_idx_ < bin)
      last_idx_ = bin;
  }

  // set max_intensity
  max_intensity_ = max_intensity;

  NormalizeIntensityArray(100);

  SmoothPeaks();

  ZeroPeaks();

  // Sequest28 modifications.
  // Determine number of top peaks to select based on the experimental mass
  // In Sequest27, the top peaks were always selected as 200.
  // Keep top ions of square-root(16*experimental mass) ranking, but not 
  // exceeding 200 ions.
  int top_bins = 200;
  if(experimental_mass_cut_off-50 < 3200){
    // top bins are sqrt of 16* experimental mass
    top_bins = (int)(sqrt((experimental_mass_cut_off-50)*16) + 0.5); 
    // however cannot exceed 200
    if(top_bins > 200){
      top_bins = 200;
    }
  }
  else{
    top_bins = (int)((experimental_mass_cut_off-50)/14.00);
  }

  // extract the top ions
  ExtractPeaks(top_bins);

  // equalize peaks
  EqualizePeaks();
}

void SpSpectrum::NormalizeIntensityArray(double threshold) {
    // return if max_intensity is 0
    if(max_intensity_ < 0.00001)
      return;
 
    // normalize all peaks
    for(int mz_idx = 0; mz_idx < last_idx_ + 1; ++mz_idx) {
      intensity_array_[mz_idx] = intensity_array_[mz_idx] * 
                                    threshold / max_intensity_;
    }
}

void SpSpectrum::SmoothPeaks() {
  // create a new array, which will replace the original intensity array
  double* new_array = new double[IntensityArraySize()];
  memset(new_array, 0, sizeof(double)*IntensityArraySize());

  // iterate over all peaks
  for(int idx = 2; idx < IntensityArraySize()-2; ++idx) {
    // smooooth
    new_array[idx] = (intensity_array_[idx-2] + 4*intensity_array_[idx-1] + 
                      6*intensity_array_[idx] + 4*intensity_array_[idx+1] +
                      intensity_array_[idx+2])/16;

    // set last idx in the array
    if(last_idx_ < idx && new_array[idx] == 0) {
      last_idx_ = idx -1;
      break;
    }
  }

  delete [] intensity_array_;
  intensity_array_ = new_array;
}

void SpSpectrum::ZeroPeaks() {
  // create a new array, which will replace the original intensity array
  double* new_array = new double[IntensityArraySize()];
  memset(new_array, 0, sizeof(double)*IntensityArraySize());
  
  // step 1,
  ZeroPeakMeanStdev(1, new_array);

  // step 2,
  ZeroPeakMeanStdev(2, new_array);

  delete[] intensity_array_;
  intensity_array_ = new_array;
}

void SpSpectrum::ExtractPeaks(int top_rank) {
  // copy all peaks to temp_array
  double* temp_array = new double[IntensityArraySize()];
  memset(temp_array, 0, sizeof(double)*IntensityArraySize());
  int temp_idx = 0;
  for(int idx = 0; idx < IntensityArraySize(); ++idx){
    if(intensity_array_[idx] > 0){
      temp_array[temp_idx] = intensity_array_[idx];
      ++temp_idx;
    }
  }
  
  // if there's over top_rank peaks, keep only top_rank peaks
  // quick sort
  QuickSort(temp_array, temp_idx);
  
  // set max and cut_off
  double max_intensity = temp_array[0];
  double cut_off = temp_array[top_rank-1];
  
  // remove peaks bellow cut_off 
  // also, normalize peaks to max_intensity to 100
  for(int idx = 0; idx < IntensityArraySize(); ++idx) {
    if(intensity_array_[idx] > 0) {
      if(intensity_array_[idx] < cut_off) {
        // If it's below cutoff, remove peak
        intensity_array_[idx] = 0;
      } else {
        // nomalize peak to max 100
        intensity_array_[idx] = intensity_array_[idx] / max_intensity * 100;
      }
    }
  }
  
  delete [] temp_array;
}

void SpSpectrum::EqualizePeaks()
{
  // equalize peaks to it's greatest intensity
  // should use array size, but sequest seems to have a bug
  // last idx is thus, modification to fit sequest
  // consequences are we will not equalize the very last peak.
  double max_intensity = 0;
  int end_idx = 0;
  int last_idx = IntensityArraySize();
  for(int idx = 0; idx < last_idx; ++idx) {
    // are we inside a continous block?
    if(intensity_array_[idx] > 0){
      max_intensity = intensity_array_[idx];
      end_idx = idx + 1;
      
      // loop to find the largest peak in the continuous block
      while(end_idx < last_idx && intensity_array_[end_idx] > 0) {
        // reset max intensity
        if(intensity_array_[end_idx] > max_intensity)
           max_intensity = intensity_array_[end_idx];
        
        ++end_idx;
      }
      
      // set all peaks in block to max_intensity
      for(; idx < end_idx; ++idx)
        intensity_array_[idx] = max_intensity;
    }
  }
}

double SpSpectrum::GetIntensityArrayMean(int peak_idx, int* peak_count) {
  // set upper bound
  int end_idx = peak_idx + 50;
  if(peak_idx + 50 >= IntensityArraySize())
    end_idx = IntensityArraySize()-1;
  
  // set start index
  int start_idx = peak_idx - 50;
  if(peak_idx - 50 <= 0)
    start_idx = 0;
  
  // sum up the intensities
  double total_intensity = 0;
  for(; start_idx <= end_idx; ++start_idx) {
    ++*peak_count;
    total_intensity += intensity_array_[start_idx];
  }

  // BUG! it should divide by 101 but Sequest uses 100
  return (total_intensity / (*peak_count-1));
}

double SpSpectrum::GetIntensityArrayStdev(int peak_idx, double mean, 
                                        int peak_count) {
  // set upper bound
  int end_idx = peak_idx + 50;
  if(peak_idx + 50 >= IntensityArraySize())
    end_idx = IntensityArraySize()-1;

  // set start index
  int start_idx = peak_idx - 50;
  if(peak_idx - 50 <= 0)
    start_idx = 0;
  
  // sum up the intensities
  double variance = 0;
  double dev = 0;
  for(; start_idx <= end_idx; ++start_idx){
    // sum up all deviations squared
    dev = intensity_array_[start_idx] - mean;
    variance += (dev*dev);
  }
  
  // return the stdev
  return sqrt(variance/peak_count);
}

void SpSpectrum::ZeroPeakMeanStdev(int step, double* new_array) {
  // iterate over all peaks
  double mean = 0;
  double stdev = 0;
  int peak_count = 0;
  for(int idx = 0; idx < IntensityArraySize(); ++idx){
    peak_count = 0;
  
    // get mean
    mean = GetIntensityArrayMean(idx, &peak_count);

    // get stdev
    stdev = GetIntensityArrayStdev(idx, mean, peak_count);
    
    // iterate over all positions and extract peaks
    if(intensity_array_[idx] > (mean + step*stdev)) {
      new_array[idx] = intensity_array_[idx] - (mean - stdev);

      // reset the last idx
      if(last_idx_ < idx)
        last_idx_ = idx;

      // for only step 1, zero out original peak
      if(step == 1)
        intensity_array_[idx] = 0;
    }
  }
}
