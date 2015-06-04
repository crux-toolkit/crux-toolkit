// Benjamin Diament
// For handling the arcana and complexity of limiting m/z.
// Most of this comes up in an effort to mimic Crux's results prefectly, and 
// generally concerns at which m/z value we should stop loops. If there is any 
// theme to draw from these arcana it seems that 0's are implicitly tacked on to the 
// high-end of a spectrum after a cut-off is programatically selected.
//
// Many of the methods of the class MaxMZ concern converting back and forth 
// between m/z values (in daltons/charge) to buckets.
//
// Once a MaxMZ is Init()'d the following three functions return
// appropriate calculated values for the spectrum under consideration:
// MaxBin() : highest bin value that will be returned by MaxMZ::Bin().
// BackgroundBinEnd() : one bin beyond where background subtraction should end.
// CacheBinEnd() : one bin beyond where the cache lookups can extend.


// This new class called MZX_BIN is to replace the old MAX_MZ class.
//written by Attila Kertesz-Farkas

#ifndef MAX_BIN_H
#define MAX_BIN_H

#include <assert.h>
#include <gflags/gflags.h>
#include "mass_constants.h"

DECLARE_double(max_bin);

// Maximum offset in the denominator of the XCorr function.
const int MAX_XCORR_OFFSET = 75;

class MaxBin {
 public:
  MaxBin() : max_bin_(0), background_bin_end_(0), cache_bin_end_(0) {}

  //MassConstants has to be initialized before this call.
  void InitBin(int highest_mz) {
    max_bin_ = MassConstants::mass2bin(highest_mz, 1);
    background_bin_end_ = MassConstants::mass2bin(highest_mz + MAX_XCORR_OFFSET + 1, 1);
    // room for diff between main ion and losses, A ion being most
    // distant, and 29 bins being as far away as an A-ion can be
    cache_bin_end_ = MassConstants::mass2bin(highest_mz + MAX_XCORR_OFFSET + 30, 1);
  }

  int MaxBinEnd() const { return max_bin_; }
  int BackgroundBinEnd() const { return background_bin_end_; }
  int CacheBinEnd() const { return cache_bin_end_; }

  static const MaxBin& Global() {
    // during indexing (as opposed to searching), Global().MaxBin() may be 0,
    // which indicates we should index all theoretical peaks, regardless of
    // m/z
    return global;
  }

  static void SetGlobalMax(double highest_mz) {
    global.InitBin(highest_mz);
    FLAGS_max_bin = global.MaxBinEnd();
  }

  static void SetGlobalMaxFromFlag() {
    if (FLAGS_max_bin > 0)
    SetGlobalMax(FLAGS_max_bin);
  }

 private:
  int max_bin_;
  int background_bin_end_;
  int cache_bin_end_;

  static MaxBin global;
};

#endif // MAX_BIN_H
