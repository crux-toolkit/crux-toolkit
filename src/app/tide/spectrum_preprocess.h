// Benjamin Diament

// ObservedPeakSet is a workspace for preprocessing an observed spectrum,
// passed in as a protocol buffer to PreprocessSpectrum().
//
// Memory allocated upon construction may be reused for multiple spectra,
// provided that the spectrum adheres to the limits of MaxMZ::Global().
// (This guarantee is provided at search time by code that finds the maximimum
// values as a side effect of sorting. In search.cc the code is this:
//    spectra.Sort();
//    MaxMZ::SetGlobalMax(spectra.HighestMZ());
// If such a maximumum weren't known in advance the client would construct
// a fresh ObervedPeakSet.)
//
// PreprocessSpectrum() first performs a normalization procedure, then
// performs some scaling operations and linear combinations and caches the
// results.  The normalization procedure is as described in the XCORR
// literature: the sqrt of each peak intensity is taken, then peaks are
// divided into 10 regions and normalized to have a maximum intensity of 50
// within each region.
//
// Notionally what we wish to do is to take a dot product of the normalized peak
// vector with a potentially large collection of different theoretical peak
// vectors, one dot product for each candidate peptide. In the interest of
// efficiency, and because the structure of the theoretical peak vector has some
// regularities, we can gain efficiency by caching certain results.  The cache,
// accessed by GetCache(), consists of various transformations of the normalized
// peak vector.
//
// These efficiencies are availlable because certain patterns are seen in these
// theoretical spectra over and over. For instance, whenever we see a Y ion at a
// particular m/z bin in the theoretical spectrum it will have intensity of 50,
// and the two bins flanking it will have an intensity of 25. If this Y ion were
// of charge 1 then a neutral loss of ammonia will almost always appear 17 bins
// to the left with intensity 10.  If this Y ion were of charge 2 then the
// theoretical spectrum will have a peak of intensity 10 at either 8 or 9 bins
// to the left -- both versions occur commonly in different theoretical spectra.
// If this ion were a B ion of charge 1 instead, we would additionally expect a
// corresponding A ion 28 bins to the left in the theoretical spectrum and a
// neutral loss of water 18 bins to the left.  Because these patterns are so
// common, it pays very well to compute these transformations and cache the
// results.
//
// Let's call the normalized peak vector u. The cache itself contains several
// vectors, corresponding to each of the TheoreticalPeakType's (see
// theoretical_peak_pair.h). The first few are straightforward:
//
// PeakMain --> u
// LossPeak --> 10 * u
// FlankingPeak --> 25 * u
// PrimaryPeak --> 50 * u
//
// Actually, instead of computing 10 * u, 25 * u, and 50 * u, we compute 2 * u,
// 5 * u and 10 * u. This results in dot products that are 5 times too small,
// but the adjustments can be made at the last moment e.g. when results are
// displayed. These smaller multiplications allow us to use addition operations
// instead of multiplications.
//
// The remaining cache vectors are linear combinations. The idea is that
// multiplying the ith entry of the appropriate vector below by the ith entry in
// a theoretical vector of only B and Y ions will capture the contribution of
// the flanks and the neutral lossses and A ion all in one operation.
//
// PeakCombinedB1 represents a charge 1 B ion. u[i-1] and u[i+1] are the two
// flanking bins. u[i-17] represents the loss of NH3 at charge 1, u[i-18]
// represents the loss of H2O  at charge 1, and u[i-28] represents the A ion at
// charge 1. It's useful to cache this vector because, most often, a
// TheoreticalPeakSet will need to represent all six of these peaks at these
// intensities and relative displacements whenever a B ion of charge 1 is
// present. Moreover, in the rare case where one or two of these displacements
// are off by one, say, with respect to the actual TheoreticalPeakSet, a
// difference vector reprsented by TheoreticalPeakSetDiff (q.v.) accommodates
// the discrepancy.
// The ith entry of this cache vector is:
//   50*u[i] + 25*u[i-1] + 25*u[i+1] + 10*u[i-17] + 10*u[i-18] + 10*u[i-28]
//
// Cache vetors for other common theoretical peak conformations are defined
// similarly:
//
// PeakCombinedY1 represents a charge 1 Y ion, its flanks and neutral losses.
// The ith entry of this cache vector is:
//   50*u[i] + 25*u[i-1] + 25*u[i+1] + 10*u[i-17]
//
// Each of the charge 2 peak types below includes an 'a' version and a 'b'
// version. In the 'a' version the loss of ammonia is 9 Da. per charge.  In the
// 'b' version the loss of ammonia is 8 Da. per charge. This is because both of
// these versions occur commonly.
//
// PeakCombinedB2a represents a charge 2 B ion, its flanks and neutral losses.
// The ith entry of this cache vector is:
//   50*u[i] + 25*u[i-1] + 25*u[i+1] + 10*u[i-9] + 10*u[i-14]
//
// PeakCombinedY2a represents a charge 2 Y ion, its flanks and neutral losses.
// The ith entry of this cache vector is:
//   50*u[i] + 25*u[i-1] + 25*u[i+1] + 10*u[i-9]
//
// PeakCombinedB2b represents a charge 2 B ion, its flanks and neutral losses.
// The ith entry of this cache vector is:
//   50*u[i] + 25*u[i-1] + 25*u[i+1] + 10*u[i-8] + 10*u[i-9] + 10*u[i-14]
//
// PeakCombinedY2b represents a charge 2 Y ion, its flanks and neutral losses.
// The ith entry of this cache vector is:
//   50*u[i] + 25*u[i-1] + 25*u[i+1] + 10*u[i-8]
#ifndef SPECTRUM_PREPROCESS_H
#define SPECTRUM_PREPROCESS_H

#include <iostream>
#include <vector>
#include "theoretical_peak_set.h"
#include "max_mz.h"
#include "mass_constants.h"

using namespace std;

class Spectrum;

class ObservedPeakSet {
 public:

  // This constructor is to be removed when marging TideSearchLite into the main branch.
/*  ObservedPeakSet(double bin_width = MassConstants::bin_width_,
     double bin_offset = MassConstants::bin_offset_,
     bool NL = false, bool FP = false)
    : peaks_(NULL),
    cache_(NULL) {

    bin_width_  = bin_width;
    bin_offset_ = bin_offset;
    NL_ = NL; //NL means neutral loss
    FP_ = FP; //FP means flanking peaks
    background_bin_end_ = 0;
    cache_end_ = 0;
  }

  // Keep this constructor when TideSearchLite is accepted.
*/  ObservedPeakSet( bool NL = false, bool FP = false) {
    peaks_ = NULL;
    cache_ = NULL;

    bin_width_  = MassConstants::bin_width_;
    bin_offset_ = MassConstants::bin_offset_;
    NL_ = NL; //NL means neutral loss
    FP_ = FP; //FP means flanking peaks
    background_bin_end_ = 0;
    cache_end_ = 0;
  }

//  ~ObservedPeakSet() { delete[] peaks_; delete[] cache_; }
  ~ObservedPeakSet() { 
    if (peaks_ != NULL)
      delete peaks_;
    if (cache_ != NULL)
    delete[] cache_; 
  }

  const int* GetCache() const { return cache_; } //TODO 261: access restriction?

  void PreprocessSpectrum(const Spectrum& spectrum, int charge) {
    long int dummy1, dummy2, dummy3, dummy4;
    PreprocessSpectrum(spectrum, charge, &dummy1, &dummy2, &dummy3, &dummy4);
  }

  void PreprocessSpectrum(const Spectrum& spectrum, int charge,
                          long int* num_range_skipped,
                          long int* num_precursors_skipped,
                          long int* num_isotopes_skipped,
                          long int* num_retained,
                          bool dia_mode = false);

  // created by Andy Lin 2/11/2016
  // Method for creating residue evidence matrix from Spectrum
  void CreateResidueEvidenceMatrix(const Spectrum& spectrum,
                                   int charge,
                                   int maxPrecurMassBin,
                                   double precursorMass,
                                   int nAA,
                                   const vector<double>& aaMass,
                                   double fragTol, int granularityScale,
                                   double nTermMass, double cTermMass,
                                   long int* num_range_skipped,
                                   long int* num_precursors_skipped,
                                   long int* num_isotopes_skipped,
                                   long int* num_retained,
                                   vector< vector<double> >& residueEvidenceMatrix);
// created by Andy Lin in Feb 2018
// help method for CreateResidueEvidenceMatrix
void addEvidToResEvMatrix(vector<double>& ionMass,
                    vector<int>& ionMassBin,
                    vector<double>& ionMasses,
                    vector<double>& ionIntens,
                    vector<double>& ionIntensitiesSort,
                    int numSpecPeaks,
                    int maxPrecurMassBin,
                    int nAA,
                    const vector<double>& aaMass,
                    const vector<int>& aaMassBin,
                    const double residueToleranceMass,
                    vector< vector<double> >& residueEvidenceMatrix);


  // added by Yang
  int LargestMzbin() const { return largest_mzbin_; };
  int SmallestMzbin() const { return smallest_mzbin_; };
  vector<pair<int, double>>& FilteredPeakTuples() { return dyn_filtered_peak_tuples_; }
  int getBackgroundBinEnd() {return background_bin_end_; }
  int getCacheEnd() {return cache_end_; }

  void SpectrumTopN(const Spectrum& spectrum, size_t n, int charge,
                          long int* num_range_skipped,
                          long int* num_precursors_skipped,
                          long int* num_isotopes_skipped,
                          long int* num_retained,
                          bool dia_mode = false);
  std::vector<std::pair<unsigned int, double>> top_N_peaks_; 

 private:
  void PreparePeaks(const Spectrum& spectrum, int charge);
  void FillCache();
  double FillPeaks(const Spectrum& spectrum, int charge,
                          long int* num_range_skipped,
                          long int* num_precursors_skipped,
                          long int* num_isotopes_skipped,
                          long int* num_retained);
  
  void KeepTopNPeaks(size_t n);
  void PreprocessSpectrum(const Spectrum& spectrum, double* intensArrayObs,
                          int* intensRegion, int maxPrecurMass, int charge);

  double* peaks_;
  int* cache_;

  bool NL_;
  bool FP_;
  double bin_width_;
  double bin_offset_;
  int background_bin_end_;
  int cache_end_;

  // added by Yang
  vector<pair<int, double>> dyn_filtered_peak_tuples_;
  int largest_mzbin_, smallest_mzbin_;

  friend class ObservedPeakTester;
};

#endif

