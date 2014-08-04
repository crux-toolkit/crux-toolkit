// Tahmina Baker
//
// The code in this file is pretty much a direct copy  of the 2009 version of
// crux's "create_intensity_array_sp" and related functions in scorer.c.
//
// The code preprocesses a spectrum for sp calculations as follows:
// 1. Create an intensity array for the spectrum to be indexed via m/z bin
// 2. Normalize the array so that maximum peak equals 100
// 3. Smooth the peaks in the array
// 4. Extract peaks that are larger than mean + #step*stdev into new array,
//    zero out the peaks that have been extracted, repeat twice, then replace 
//    old array with extracted peak array
// 5. Extract peaks up to top rank peaks remove other peaks, do second 
//    normalization on the top peaks back to max 100 intensity, and replace old
//    array with normalized top peak array.
// 6. Equalize all peaks in a continous region to the largest peak within the
//    continous bins from left to right.
//
// The following major modifications were made when the code was ported over:
// 1. Functions were renamed to follow Tide conventions. 
//    E.g. "PreprocessSpectrum" performs the same function as crux's 
//    "create_intensity_array_sp". 
// 2. Crux uses the data structure SCORER_T to calculate different types of
//    scores represented by SCORER_TYPE_T, whereas the code here is specific
//    to sp calculations. The SpSpectrum class is equivalent to a crux SCORER_T
//    struct of SCORER_TYPE_T=SP, with only the following subset of data: 
//    sp_beta, sp_max_mz, intensity_array, max_intensity, last_idx.
// 3. SpSpectrum is a private class that preprocesses a given spectrum on 
//    creation using private methods, then ONLY exposes the pieces of 
//    information needed by the creator. In crux, the static function 
//    "create_intensity_array_sp" processes a given spectrum with other static
//    helper functions, then populates an SCORER_T public struct.
// 4. Tide stores spectrum peak information slightly differently than crux; 
//    changes were made to the code to traverse peaks and access peak data.
// 5. SpSpectrum::SmoothPeaks() is specific to SP calculations while crux's
//    "smooth_peaks" function can handle various types of scores represented 
//    by SCORER_TYPE_T.
// 6. The "double" data type was used in place of crux's "FLOAT_T".


#ifndef CRUX_SP_SPECTRUM_H
#define CRUX_SP_SPECTRUM_H

#include "mass_constants.h"
#include "spectrum_collection.h"
#include "spectrum_preprocess.h"
#include "utils.h"


class SpSpectrum {
 public:
  SpSpectrum(const Spectrum& spectrum, int charge, double max_mz);
  ~SpSpectrum();

  double Intensity(int index) const { return intensity_array_[index]; }
  double Beta() const { return beta_; }
  double TotalIonIntensity() {
    double total_ion_intensity = 0.0;
    for (int i = 0; i < IntensityArraySize(); i++)
      total_ion_intensity += intensity_array_[i];
    return total_ion_intensity;
  }

 private:

  void PreprocessSpectrum(const Spectrum& spectrum, int charge);

  // Normalizes array so that maximum peak equals threshold
  void NormalizeIntensityArray(double threshold);

  // Smooth all peaks in intensity array and replace the original array with 
  // the newly smooothed array.
  void SmoothPeaks();

  // Zero and extract peaks: 
  // Extract peaks that are larger than mean + #step*stdev into new array
  // Zero out the peaks that have been extracted
  // Repeat twice, then replace old array with extracted peak array
  void ZeroPeaks();

  // Keep only the peaks up to top rank peaks remove other peaks.
  // Do second normalization on the top peaks back to max 100 intensity.
  // Replace old array with normalized top peak array.
  void ExtractPeaks(int top_rank);

  // Equalize all peaks in a continous region to the largest peak within the
  // continous bins from left to right
  void EqualizePeaks();
  
  // Get the mean of intensity in array within +/- 50 mz of the working peak.
  double GetIntensityArrayMean(int peak_idx, int* peak_count);

  // Get the stdev of intensity in array within +/- 50 mz of the working peak
  double GetIntensityArrayStdev(int peak_idx, double mean, int peak_count);

  // Helper to "zero and extract peaks". The fact that a peak has removed will 
  // affect the following peaks.
  void ZeroPeakMeanStdev(int step, double* new_array);

  int IntensityArraySize() {
    return (int)max_mz_ + 1;
  }

  void SwapQuick(double* a, int idx, int jdx) {
    double temp = 0;
    temp = a[idx];
    a[idx] = a[jdx];
    a[jdx] = temp;
  }
   
  int Random(int i, int j) {
    return i + myrandom_limit(j-i+1);
  }

  void QuickSort(double a[], int left, int right) {
    int last = left, i;

    if (left >= right) return;
    
    SwapQuick(a,left,Random(left,right));
    for (i = left + 1; i <= right; i++)
      if (a[i] > a[left]) /// CHECK THIS!!
        SwapQuick(a,++last,i);
    SwapQuick(a,left,last);
    QuickSort(a,left,last-1);
    QuickSort(a,last+1,right);
  }

  void QuickSort(double a[], int array_size) {
    QuickSort(a, 0, array_size-1);
  }
  
  // The beta variable = 0.075
  double beta_; 
  
  // The max mz for the intensity array
  double max_mz_; 

  // Intensity array that can be indexed using the m/z bin
  double* intensity_array_; 

  // the max intensity in the intensity array
  double max_intensity_; 

  // the last index in the array, the data size of the array
  int last_idx_; 
};

#endif // CRUX_SP_SPECTRUM_H
