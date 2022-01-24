// Benjamin Diament and Attila Kertesz-Farkas

// Originally implemented by Benjamin Diament. It was reimplemented by Attila 
// Kertesz-Farkas in December 2021.
// Benjamin's implementetion used a structure to keep track of the types of the
// theoretical and the experimental peaks. There were 10 types of peaks such as:
// PeakMain = 0,
// LossPeak = 1,
// FlankingPeak = 2,
// PrimaryPeak = 3,
// PeakCombinedB1 = 4, // Represents charge 1 B ion, its flanks and neutral losses.
// PeakCombinedY1 = 5, // Represents charge 1 Y ion, its flanks and neutral losses.
// // Each of the charge 2 peak types below includes an 'a' version and a 'b'
// // version. In the 'a' version the loss of ammonia is 9 Da. per charge.
// // In the 'b' version the loss of ammonia is 8 Da. per charge.
// PeakCombinedB2 = 6, // Represents charge 2 B ion, its flanks and neutral losses.
// PeakCombinedY2 = 7, // Represents charge 2 Y ion, its flanks and neutral losses.
// //  PeakCombinedB2b = 8, // Represents charge 2 B ion, its flanks and neutral losses.
// //  PeakCombinedY2b = 9, // Represents charge 2 Y ion, its flanks and neutral losses.
// NUM_PEAK_TYPES = 10 // Total number of peak types.
//
// In the current implementation (by AKF) these types for the theoretical and 
// experimental peaks have been removed. Note that, these types are not need for the
// original XCorr implementation; however, these types would be hand for other scorings
// such as the hyperScore of X!Tandem.
// 
// The current implementation has also removed the other classes for generating theoretical 
// peaks, those classes were not used. 

#ifndef THEORETICAL_PEAK_SET_H
#define THEORETICAL_PEAK_SET_H

#include <iostream>
#include <algorithm>
#include "mass_constants.h"
#include "peptides.pb.h"
#include "max_mz.h"
#include "math.h"
#include <functional>
#include "fixed_cap_array.h"

using namespace std;
typedef google::protobuf::RepeatedField<int>::const_iterator FieldIter;

typedef FixedCapacityArray<int> TheoreticalPeakArr;
#define NUM_PEAK_TYPES 2

// A class that generates a theoretical peak for B and Y ion at charge states
// of 1+ and 2+. The Original XCorr requires that the thoeretical peaks are 
// unique and this, indeed, boosts the XCorr's performance. That is, a 
// theoretical peak, say, in bin 459 for B-ion cannot score again for in 
// the same bin 459 for a Y-ion of any charges. In order to keep track of
// the unique peaks, this class employes a mask verctor to check if a peak
// (of any type (b or y) at any charge states (1+ or 2+) has already been 
// implemented). 
//
// The peaks are stored sparesly as a list in a fixed capacity vector, 
// respectively for each charge states. Note that, theoretical ions of 
// charge 1+ must be generated first, before generating the ions of carge 2+.

class TheoreticalPeakSetBYSparse {
 public:
  explicit TheoreticalPeakSetBYSparse(int capacity) {
    peaks_[0].Init(capacity);
    peaks_[1].Init(capacity);
    peak_mask_end = MaxBin::Global().CacheBinEnd();
    cache_end = MaxBin::Global().CacheBinEnd()*NUM_PEAK_TYPES;
    peak_mask = new int[peak_mask_end];
    memset(peak_mask, 0,  sizeof(int)*peak_mask_end);
  }

  virtual ~TheoreticalPeakSetBYSparse() {}

  void Clear() {
    TheoreticalPeakArr::iterator itr; 
    int idx;
    for (itr = peaks_[0].begin(); itr != peaks_[0].end(); ++itr) {
        idx = int((*itr)/2);
        peak_mask[idx]   = 0;
        peak_mask[idx+1] = 0;
    }
    for (itr = peaks_[1].begin(); itr != peaks_[1].end(); ++itr) {
        idx = int((*itr)/2);
        peak_mask[idx]   = 0;
        peak_mask[idx+1] = 0;
    }
    peaks_[0].clear();
    peaks_[1].clear();
  }

  void AddYIon(double mass, int charge) {
    assert(charge <= 2);
    int index_y = MassConstants::mass2bin(mass + MassConstants::Y + MASS_PROTON, charge);
    if (index_y >= 0 && index_y < peak_mask_end && !peak_mask[index_y]) {
      peak_mask[index_y] = 1;
      index_y += index_y;
      if (charge == 2) {
        ++index_y;
      }  
      if (index_y < cache_end)
        peaks_[charge-1].push_back(index_y);
    }
  }

  void AddBIon(double mass, int charge) {
    assert(charge <= 2);
    int index_b = MassConstants::mass2bin(mass + MassConstants::B + MASS_PROTON, charge);
    if (index_b >= 0 && index_b < peak_mask_end && !peak_mask[index_b]) {
      peak_mask[index_b] = 1;
      index_b += index_b;
      if (charge == 2) {
        ++index_b;
      } 
      if (index_b < cache_end)
        peaks_[charge-1].push_back(index_b);
    }
  }
  // Faster interface needing no copying at all.
  const TheoreticalPeakArr* GetPeaks() const { return peaks_; }

 private:
    int* peak_mask = 0; //MaxBin::Global().BackgroundBinEnd()
    int peak_mask_end = 0;
    int cache_end = 0;
    TheoreticalPeakArr peaks_[2];
};

// This class is used to store theoretical b ions only
// for use in exact p-value calculations.
class TheoreticalPeakSetBIons {
 public:
  TheoreticalPeakSetBIons() {}
  explicit TheoreticalPeakSetBIons(int capacity) {
    unordered_peak_list_.reserve(capacity);
  }
  virtual ~TheoreticalPeakSetBIons() {}

  void Clear() { unordered_peak_list_.clear(); }
  void AddBIon(double mass) {
    int index = MassConstants::mass2bin(mass, 1);
    if (index >= 0) {
      unordered_peak_list_.push_back((unsigned int)index);
    }
  }
  vector<unsigned int> unordered_peak_list_;
  double binWidth_;   // TODO: This should be removed from here.
  double binOffset_;  // TODO: This should be removed from here.
};


#endif // THEORETICAL_PEAK_SET_H

