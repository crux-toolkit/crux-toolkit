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

#ifndef MAX_MZ_H
#define MAX_MZ_H

#include <assert.h>
#include <gflags/gflags.h>

DECLARE_double(max_mz);

// Bin width for integerizing m/z axis.  This value was taken directly
// from SEQUEST.
const double bin_width_mono = 1.0005079;

// Maximum offset in the denominator of the XCorr function.
const int MAX_XCORR_OFFSET = 75;

class MaxMZ {
 public:
  MaxMZ() : max_bin_(0), background_bin_end_(0), cache_bin_end_(0) {}

  static int Bin(double mz) {
    return int(mz/bin_width_mono + 0.5);
  }

  static double BinInvert(int bin) {
    // returns highest m/z that would be binned (by Bin() above) into bin.
    // result should err on the side of being too large, hence extra 0.005
    return (double(bin) + 0.505) * bin_width_mono;
  }

  void Init(double highest_mz) {
    InitBin(Bin(highest_mz));
  }

  void InitBin(int highest_bin) {
    max_bin_ = highest_bin;
    background_bin_end_ = max_bin_ + MAX_XCORR_OFFSET + 1;
    // room for diff between main ion and losses, A ion being most
    // distant, and 29 bins being as far away as an A-ion can be
    cache_bin_end_ = background_bin_end_ + 29;
  }

  int MaxBin() const { return max_bin_; }
  int BackgroundBinEnd() const { return background_bin_end_; }
  int CacheBinEnd() const { return cache_bin_end_; }

  static const MaxMZ& Global() {
    // during indexing (as opposed to searching), Global().MaxBin() may be 0,
    // which indicates we should index all theoretical peaks, regardless of
    // m/z
    return global;
  }

  static void SetGlobalMax(double highest_mz) {
    assert(global.max_bin_ == 0); // this should be done only once
    FLAGS_max_mz = highest_mz;
    global.Init(highest_mz);
  }

  static void SetGlobalMaxFromFlag() {
    if (FLAGS_max_mz > 0);
    SetGlobalMax(FLAGS_max_mz);
  }

 private:
  int max_bin_;
  int background_bin_end_;
  int cache_bin_end_;

  static MaxMZ global;
};

#endif // MAX_MZ_H
