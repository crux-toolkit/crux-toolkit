// Benjamin Diament

// A TheoreticalPeakPair represents the ordered pair (bucket, type) where
// bucket is an integer refering to a bucketed m/z position in the theoretical
// spectrum, and type identifies the type of the peak. This type may be
// LossPeak, FlankingPeak, PrimaryPeak which represents an intensity, or it could be one of
// the combined peak types (see TheoreticalPeakType below) which represent a 
// linear combination of several peaks.
#ifndef THEORETICAL_PEAK_PAIR_H
#define THEORETICAL_PEAK_PAIR_H

#include <utility>
#include "fixed_cap_array.h"

using namespace std;

enum TheoreticalPeakType {
  PeakMain = 0,
  LossPeak = 1,
  FlankingPeak = 2,
  PrimaryPeak = 3,
  PeakCombinedB1 = 4, // Represents charge 1 B ion, its flanks and neutral losses.
  PeakCombinedY1 = 5, // Represents charge 1 Y ion, its flanks and neutral losses.
  // Each of the charge 2 peak types below includes an 'a' version and a 'b'
  // version. In the 'a' version the loss of ammonia is 9 Da. per charge.
  // In the 'b' version the loss of ammonia is 8 Da. per charge.
  PeakCombinedB2 = 6, // Represents charge 2 B ion, its flanks and neutral losses.
  PeakCombinedY2 = 7, // Represents charge 2 Y ion, its flanks and neutral losses.
//  PeakCombinedB2b = 8, // Represents charge 2 B ion, its flanks and neutral losses.
//  PeakCombinedY2b = 9, // Represents charge 2 Y ion, its flanks and neutral losses.

  NUM_PEAK_TYPES = 10 // Total number of peak types.
};

class TheoreticalPeakPair {
 public:
  TheoreticalPeakPair() {} // no intialization; beware.

  TheoreticalPeakPair(int code) : code_(code) {}

  TheoreticalPeakPair(int bin, TheoreticalPeakType peak_type)
    : code_(bin * NUM_PEAK_TYPES + peak_type) {
  }

  // Bin() and Type() give the components of the ordered pair as described
  // above. Code() should be used for serialization to and from a file, and
  // for indexing into the precomputed cache of transformed observed peaks.
  // Note that Bin() and Type() are not fast; Code() is fast.
  int Bin() const { return code_ / NUM_PEAK_TYPES; }
  int Type() const { return code_ % NUM_PEAK_TYPES; }
  int Code() const { return code_; }

  bool operator<(const TheoreticalPeakPair& other) const {
    return code_ < other.code_;
  }

  bool operator>(const TheoreticalPeakPair& other) const {
    return code_ > other.code_;
  }

 private:
   // code_ represents the bin position and the peak type together. See
   // constructors and accessors.
  int code_;
};

typedef FixedCapacityArray<TheoreticalPeakPair> TheoreticalPeakArr;

#endif // THEORETICAL_PEAK_PAIR_H
