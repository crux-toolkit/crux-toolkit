// Benjamin Diament
//
// Mostly an implementation of ObservedPeakSet::PreprocessSpectrum(). See .h
// file.
//
// Portions of this code are taken verbatim from Crux, specifically
// NormalizeEachRegion(). The region-defining code in PreprocessSpectrum() is
// also from Crux and inherits the attendant arcana. Some of this code is
// certainly at odds with the published intention of XCORR, but in order to
// preserve the results from Crux exactly we retain this legacy for now.
// TODO 251: revisit this.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <algorithm>
#include <math.h>
#include <gflags/gflags.h>
#include "spectrum_collection.h"
#include "spectrum_preprocess.h"
#include "mass_constants.h" //added by Andy Lin
#include "max_mz.h"
#include "util/mass.h"
#include "util/Params.h"
#include <cmath>

using namespace std;

#ifdef DEBUG
DEFINE_int32(debug_spectrum_id, -1, "Spectrum number to debug");
DEFINE_int32(debug_charge, 0, "Charge to debug. 0 for all");
#endif

// This computes that part of the XCORR function where an average value of the
// peaks within a window surrounding each peak is subtracted from that peak.
// This version is a linear-time implementation of the subtraction. Linearity is
// accomplished by computing an array of partial sums.
static void SubtractBackground(double* observed, int end) {
  // operation is as follows: new_observed = observed -
  // average_within_window but average is computed as if the array
  // extended infinitely: denominator is same throughout array, even
  // near edges (where fewer elements have been summed)
  static const double multiplier = 1.0 / (MAX_XCORR_OFFSET * 2);

  double total = 0;
  vector<double> partial_sums(end+1);
  for (int i = 0; i < end; ++i)
    partial_sums[i] = (total += observed[i]);
  partial_sums[end] = total;

  for (int i = 0; i < end; ++i) {
    int right_index = min(end, i + MAX_XCORR_OFFSET);
    int left_index = max(0, i - MAX_XCORR_OFFSET - 1);
    observed[i] -= multiplier * (partial_sums[right_index] - partial_sums[left_index] - observed[i]);
  }
}

void ObservedPeakSet::PreprocessSpectrum(const Spectrum& spectrum, int charge,
                                         long int* num_range_skipped,
                                         long int* num_precursors_skipped,
                                         long int* num_isotopes_skipped,
                                         long int* num_retained) {
#ifdef DEBUG
  bool debug = (FLAGS_debug_spectrum_id == spectrum.SpectrumNumber()
                && (FLAGS_debug_charge == 0 || FLAGS_debug_charge == charge));
  if (debug)
    debug = true; // allows a breakpoint
#endif
  double precursor_mz = spectrum.PrecursorMZ();
  double experimental_mass_cut_off = (precursor_mz-MASS_PROTON)*charge+MASS_PROTON + 50;
  double max_peak_mz = spectrum.M_Z(spectrum.Size()-1);

  assert(MaxBin::Global().MaxBinEnd() > 0);

  max_mz_.InitBin(min(experimental_mass_cut_off, max_peak_mz));
  cache_end_ = MaxBin::Global().CacheBinEnd() * NUM_PEAK_TYPES;

  memset(peaks_, 0, sizeof(double) * MaxBin::Global().BackgroundBinEnd());

  if (Params::GetBool("skip-preprocessing")) {
    for (int i = 0; i < spectrum.Size(); ++i) {
      double peak_location = spectrum.M_Z(i);
      if (peak_location >= experimental_mass_cut_off) {
        (*num_range_skipped)++;
        continue;
      }
      int mz = MassConstants::mass2bin(peak_location);
      double intensity = spectrum.Intensity(i);
      if (intensity > peaks_[mz]) {
        peaks_[mz] = intensity;
      }
    }
  } else {
    bool remove_precursor = Params::GetBool("remove-precursor-peak");
    double precursor_tolerance = Params::GetDouble("remove-precursor-tolerance");
    double deisotope_threshold = Params::GetDouble("deisotope");
    int max_charge = spectrum.MaxCharge();

    // Fill peaks
    int largest_mz = 0;
    double highest_intensity = 0;
    for (int i = spectrum.Size() - 1; i >= 0; --i) {
      double peak_location = spectrum.M_Z(i);

      // Get rid of peaks beyond the possible range, given charge and precursor.
      if (peak_location >= experimental_mass_cut_off) {
        (*num_range_skipped)++;
        continue;
      }

      // Remove precursor peaks.
      if (remove_precursor &&
          fabs(peak_location - precursor_mz) <= precursor_tolerance ) {
        (*num_precursors_skipped)++;
        continue;
      }

      if (deisotope_threshold != 0.0 && spectrum.Deisotope(i, deisotope_threshold)) {
        (*num_isotopes_skipped)++;
        continue;
      }

      (*num_retained)++;

      int mz = MassConstants::mass2bin(peak_location);
      double intensity = spectrum.Intensity(i);
      if ((mz > largest_mz) && (intensity > 0)) {
        largest_mz = mz;
      }

      intensity = sqrt(intensity);
      if (intensity > highest_intensity) {
        highest_intensity = intensity;
      }
      if (intensity > peaks_[mz]) {
        peaks_[mz] = intensity;
      }
    }

    double intensity_cutoff = highest_intensity * 0.05;

    double normalizer = 0.0;
    int region_size = largest_mz / NUM_SPECTRUM_REGIONS + 1;
    for (int i = 0; i < NUM_SPECTRUM_REGIONS; ++i) {
      highest_intensity = 0;
      int high_index = i;
      for (int j = 0; j < region_size; ++j) {
        int index = i * region_size + j;
        if (peaks_[index] <= intensity_cutoff) {
          peaks_[index] = 0;
        }
        if (peaks_[index] > highest_intensity) {
          highest_intensity = peaks_[index];
          high_index = index;
        }
      }
      if (highest_intensity == 0) {
        continue;
      }
      normalizer = 50.0 / highest_intensity;
      for (int j = 0; j < region_size; ++j) {
        int index = i * region_size + j;
        if (peaks_[index] != 0) {
          peaks_[index] *= normalizer;
        }
      }
    }

#ifdef DEBUG
    if (debug) {
      cout << "GLOBAL MAX MZ: " << MaxMZ::Global().MaxBin() << ", " << MaxMZ::Global().BackgroundBinEnd()
           << ", " << MaxMZ::Global().CacheBinEnd() << endl;
      cout << "MAX MZ: " << max_mz_.MaxBin() << ", " << max_mz_.BackgroundBinEnd()
           << ", " << max_mz_.CacheBinEnd() << endl;
      ShowPeaks();
      cout << "====== SUBTRACTING BACKGROUND ======" << endl;
    }
#endif
  }
  SubtractBackground(peaks_, max_mz_.BackgroundBinEnd());

#ifdef DEBUG
  if (debug)
    ShowPeaks();
#endif
  MakeInteger();
  ComputeCache();
#ifdef DEBUG
  if (debug)
    ShowCache();
#endif
}

inline int round_to_int(double x) {
  if (x >= 0)
    return int(x + 0.5);
  return int(x - 0.5);
}

void ObservedPeakSet::MakeInteger() {
  // essentially cheap fixed-point arithmetic for peak intensities
  for(int i = 0; i < max_mz_.BackgroundBinEnd(); i++)
    Peak(PeakMain, i) = round_to_int(peaks_[i]*50000);
}

// See .h file. Computes and stores all transformations of the observed peak
// set.
void ObservedPeakSet::ComputeCache() {
  for (int i = 0; i < max_mz_.BackgroundBinEnd(); ++i) {
    // Instead of computing 10 * x, 25 * x, and 50 * x, we compute 2 *
    // x, 5 * x and 10 * x. This results in dot products that are 5
    // times too small, but the adjustments can be made at the last
    // moment e.g. when results are displayed. These smaller
    // multiplications allow us to use addition operations instead of
    // multiplications.
    int x = Peak(PeakMain, i);
    int y = x+x;
    Peak(LossPeak, i) = y;
    int z = y+y+x;
    Peak(FlankingPeak, i) = z;
    Peak(PrimaryPeak, i) = z+z;
  }

  for (int i = max_mz_.BackgroundBinEnd() * NUM_PEAK_TYPES; i < cache_end_; ++i) {
    cache_[i] = 0;
  }

  for (int i = 0; i < max_mz_.CacheBinEnd(); ++i) {
    int flanks = Peak(PrimaryPeak, i);
    if ( FP_ == true) {
        if (i > 0) {
          flanks += Peak(FlankingPeak, i-1);
        }
        if (i < max_mz_.CacheBinEnd() - 1) {
          flanks += Peak(FlankingPeak, i+1);
        }
    }
    int Y1 = flanks;
    if ( NL_ == true) {
        if (i > MassConstants::BIN_NH3) {
          Y1 += Peak(LossPeak, i-MassConstants::BIN_NH3);
        }
        if (i > MassConstants::BIN_H2O) {
          Y1 += Peak(LossPeak, i-MassConstants::BIN_H2O);
        }
    }
    Peak(PeakCombinedY1, i) = Y1;
    int B1 = Y1;
    Peak(PeakCombinedB1, i) = B1;
    Peak(PeakCombinedY2, i) = flanks;
    Peak(PeakCombinedB2, i) = flanks;
  }
}

// This dot product is replaced by calls to on-the-fly compiled code.
int ObservedPeakSet::DotProd(const TheoreticalPeakArr& theoretical) {
  int total = 0;
  TheoreticalPeakArr::const_iterator i = theoretical.begin();
  for (; i != theoretical.end(); ++i) {
    //if (i->Code() >= cache_end_)
    //  break;
    total += cache_[i->Code()];
  }
  return total;
}

#if 0
int ObservedPeakSet::DebugDotProd(const TheoreticalPeakArr& theoretical) {
  cout << "cache_end_=" << cache_end_ << endl;
  int total = 0;
  TheoreticalPeakArr::const_iterator i = theoretical.begin();
  for (; i != theoretical.end(); ++i) {
    cout << "DotProd Lookup(" << i->Bin() << "," << i->Type() << "):";
    if (i->Code() >= cache_end_) {
      cout << "code=" << i->Code() << "past cache_end_=" << cache_end_ << "; ignoring" << endl;
      continue;
    }
    total += cache_[i->Code()];
    cout << cache_[i->Code()] << "; total=" << total << endl;
  }
  return total;
}
#endif

// Written by Andy Lin in Feb 2018
// Helper function for CreateResidueEvidenceMatrix
void ObservedPeakSet::addEvidToResEvMatrix(
  vector<double>& ionMass,
  vector<int>& ionMassBin,
  vector<double>& ionMasses,
  vector<double>& ionIntens,
  vector<double>& ionIntensitiesSort,
  int numSpecPeaks,
  int nAA,
  const vector<double>& aaMass,
  const vector<int>& aaMassBin,
  const double residueToleranceMass,
  vector<vector<double> >& residueEvidenceMatrix
  ) {
  double bIonMass; int bIonMassBin;
  for (int ion = 0; ion < ionMass.size(); ion++) {
    // bIonMass is named correctly because we assume that all 
    // y ions have been converted to their b ion equivalent
    // before ionMass is passed in
    bIonMass = ionMass[ion];
    bIonMassBin = ionMassBin[ion];
  
    for (int curAaMass = 0; curAaMass < nAA; curAaMass++) {
      int newResMassBin = bIonMassBin + aaMassBin[curAaMass];
      
      // Find all ion mass bins that match newResMassBin
      int index = find(ionMassBin.begin(), ionMassBin.end(),newResMassBin) - ionMassBin.begin();
      double score = 0.0;
      for (int i = index; i < ionMasses.size(); i++) {
        if (newResMassBin != ionMassBin[i]) {
          continue;
        } 
       
        double ionMassDiff = ionMass[i] - bIonMass;
        //double aaTolScore = residueToleranceMass - std::abs(ionMassDiff - aaMass[curAaMass]);
        double aaTolScore = 1.0 - (std::abs(ionMassDiff - aaMass[curAaMass]) / residueToleranceMass);

        if (aaTolScore > 0.0) {
          // Find where intensities are in sorted vector
          // Based upon 10-bin normalized intensities
          // If multiple peaks have the same intensity, all 
          // peaks have the same rank
          int loc1 = std::find(ionIntensitiesSort.begin(), ionIntensitiesSort.end(), ionIntens[ion]) - ionIntensitiesSort.begin();
          int loc2 = std::find(ionIntensitiesSort.begin(), ionIntensitiesSort.end(), ionIntens[i]) - ionIntensitiesSort.begin();
  
          //determine rank intensity
          double rank1 = (double)loc1 / numSpecPeaks;
          double rank2 = (double)loc2 / numSpecPeaks;
          double tmpScore = aaTolScore * (rank1 + rank2);

          if (tmpScore > score) {
            score = tmpScore;
          }
        }
      }

      // Add evidence to matrix
      // Use -1 since all mass bins are index 1 instead of index 0
      residueEvidenceMatrix[curAaMass][newResMassBin-1] += score;
    }
  }
}

// Written by Andy Lin in Feb 2016
// From an observed spectrum, calculate and populate the residue evidence matrix
void ObservedPeakSet::CreateResidueEvidenceMatrix(
  const Spectrum& spectrum,
  int charge,
  int maxPrecurMassBin,
  double precursorMass, // neutral mass
  int nAA, //TODO different than one used in CreateEvidenceVector
  const vector<double>& aaMass, //TODO different aaMass then one used in CreateEvidenceVector
  double fragTol,
  int granularityScale,
  double nTermMass,
  double cTermMass,
  long int* num_range_skipped,
  long int* num_precursors_skipped,
  long int* num_isotopes_skipped,
  long int* num_retained,
  vector<vector<double> >& residueEvidenceMatrix
  ) {

  assert(MaxBin::Global().MaxBinEnd() > 0);

  //TODO move to constants file?
  const double massHMono = MassConstants::mono_h;  // mass of hydrogen (monoisotopic)

  //TODO used in Xcorr function..not used in residue-evidence
  //bool flanking_peak = Params::GetBool("use-flanking-peaks");
  //bool neutral_loss_peak = Params::GetBool("use-neutral-loss-peaks");

  int ma;
  int pc;

  double precurMz = spectrum.PrecursorMZ();
  int nIon = spectrum.Size();
  int precurCharge = charge;
  double experimentalMassCutoff = precursorMass + 50.0;
  double residueToleranceMass = fragTol;
  const double maxIntensPerRegion = 50.0;

  // Determining max ion mass and max ion intensity
  bool skipPreprocess = Params::GetBool("skip-preprocessing");
  bool remove_precursor = !skipPreprocess && Params::GetBool("remove-precursor-peak");
  double precursorMZExclude = Params::GetDouble("remove-precursor-tolerance");
  double deisotope_threshold = Params::GetDouble("deisotope");
  double maxIonIntens = 0.0;
  double maxIonMass = 0.0;
  set<int> peakSkip;
  for (int ion = 0; ion < nIon; ion++) {
    double ionMass = spectrum.M_Z(ion);
    double ionIntens = sqrt(spectrum.Intensity(ion));
    if (ionMass >= experimentalMassCutoff) {
      peakSkip.insert(ion);
      if (num_range_skipped) {
        (*num_range_skipped)++;
      }
      continue;
    } else if (remove_precursor && ionMass > precurMz - precursorMZExclude && 
               ionMass < precurMz + precursorMZExclude) {
      peakSkip.insert(ion);
      if (num_precursors_skipped) {
        (*num_precursors_skipped)++;
      }
      continue;
    } else if (deisotope_threshold != 0.0 && spectrum.Deisotope(ion, deisotope_threshold)) {
      peakSkip.insert(ion);
      if (num_isotopes_skipped) {
        (*num_isotopes_skipped)++;
      }
      continue;
    }

    if (num_retained) {
      (*num_retained)++;
    }
    if (maxIonIntens < ionIntens) {
      maxIonIntens = ionIntens;
    }
    if (maxIonMass < ionMass) {
      maxIonMass = ionMass;
    }
  }

  // peaks that pass preprocessing filters
  vector<double> ionMasses;
  vector<double> ionIntensities;
  vector<double> ionIntensitiesSort;
  vector<int> intensRegion;
  double numSpecPeaks;

  // grass filtering
  int regionSelector = (int)floor(MassConstants::mass2bin(maxIonMass) / (double)NUM_SPECTRUM_REGIONS);
  for(int ion = 0; ion < nIon; ion++) {
    double ionMass = spectrum.M_Z(ion);
    double ionIntens = sqrt(spectrum.Intensity(ion));

    if (peakSkip.find(ion) != peakSkip.end()) {
      continue;
    }

    if (ionIntens < 0.05 * maxIonIntens) {
      continue;
    }
    ionMasses.push_back(ionMass);
    ionIntensities.push_back(ionIntens);
    ionIntensitiesSort.push_back(ionIntens);

    int ionBin = MassConstants::mass2bin(ionMass);
    int region = (int)floor((double)(ionBin) / (double)regionSelector);
    if (region >= NUM_SPECTRUM_REGIONS) {
      region = NUM_SPECTRUM_REGIONS - 1;
    }
    intensRegion.push_back(region);
  }

  // +2 because we artificially add to the spectrum
  // a peak corresponding to the NTerm and CTerm mass
  numSpecPeaks = ionIntensities.size() + 2;

  // 10-bin intensity normalization
  vector<double> maxRegion(NUM_SPECTRUM_REGIONS, 0);
  for (int i = 0; i < ionMasses.size(); i++) {
    int curRegion = intensRegion[i];
    if (maxRegion[curRegion] < ionIntensities[i]) {
      maxRegion[curRegion] = ionIntensities[i];
    }
    if (curRegion < 0) {
      carp(CARP_FATAL, "ion is in a negative region, should never happen");
    }
  }

  for (int i = 0; i < ionMasses.size(); i++) {
    int curRegion = intensRegion[i];
    if (curRegion >= 0 && maxRegion[curRegion] > 0.0) {
      ionIntensities[i] *= (maxIntensPerRegion / maxRegion[curRegion]);
      ionIntensitiesSort[i] *= (maxIntensPerRegion / maxRegion[curRegion]);
    }
  }

  // sort ion intensities
  sort(ionIntensitiesSort.begin(), ionIntensitiesSort.end());

  // Determine which bin each amino acid mass is in
  vector<int> aaMassBin;
  for (int i = 0; i < nAA; i++) {
    int binMass = (int)floor(MassConstants::mass2bin(aaMass[i]));
    aaMassBin.push_back(binMass);
  }

  // Need to add lines for matlab line 56?
  // Basically bounds the ion masses and ion mass bins so that
  // 1: ionMassBinB >=1
  // 2: ionMassBinB <= maxMassBin - max(aaMassBin)
  // Also need to figure out maxMassBin (currently hard coded to 8500)
  // need to do 4 times (+1/+2/B ion/ yion)

  // Populate ionMass and ionMassBin for b ions in 1+ charge state
  vector<double> ionMass;
  vector<int> ionMassBin;
  vector<double> ionIntens;

  ionMass.push_back(nTermMass);
  ionMassBin.push_back(MassConstants::mass2bin(nTermMass));
  ionIntens.push_back(0.0);

  for (int ion = 0; ion < ionMasses.size(); ion++) {
    double tmpIonMass = ionMasses[ion];
    int binTmpIonMass = (int)floor(MassConstants::mass2bin(tmpIonMass));

    ionMass.push_back(tmpIonMass);
    ionMassBin.push_back(binTmpIonMass);
    ionIntens.push_back(ionIntensities[ion]);
  }
  ionMass.push_back(precursorMass - cTermMass);
  ionMassBin.push_back(MassConstants::mass2bin(precursorMass - cTermMass));
  ionIntens.push_back(0.0);

  addEvidToResEvMatrix(ionMass, ionMassBin, ionMasses, ionIntens, ionIntensitiesSort,
                       numSpecPeaks, nAA, aaMass, aaMassBin, residueToleranceMass,
                       residueEvidenceMatrix);
  ionMass.clear();
  ionMassBin.clear();
  ionIntens.clear();

  int yIonMassBin;
  double yIonMass;
  double yIonIntens;

  // Find pairs of y ions in 1+ charge state
  ionMass.push_back(precursorMass - cTermMass);
  ionMassBin.push_back(MassConstants::mass2bin(precursorMass - cTermMass));
  ionIntens.push_back(0.0);

  for (int ion = 0; ion < ionMasses.size(); ion++) {
    // Convert to equivalent b ion masses for ease of processing
    double tmpIonMass = precursorMass - ionMasses[ion] + (2.0 * massHMono);
    // Determine which bin each ion mass is in
    int binTmpIonMass = (int)floor(MassConstants::mass2bin(tmpIonMass));

    if (tmpIonMass > 0) {
      ionMass.push_back(tmpIonMass);
      ionMassBin.push_back(binTmpIonMass);
      ionIntens.push_back(ionIntensities[ion]);
    }
  }
  ionMass.push_back(nTermMass);
  ionMassBin.push_back(MassConstants::mass2bin(nTermMass));
  ionIntens.push_back(0.0);

  reverse(ionMass.begin(), ionMass.end());
  reverse(ionMassBin.begin(), ionMassBin.end());
  reverse(ionIntens.begin(), ionIntens.end());

  addEvidToResEvMatrix(ionMass, ionMassBin, ionMasses, ionIntens, ionIntensitiesSort,
                       numSpecPeaks, nAA, aaMass, aaMassBin, residueToleranceMass,
                       residueEvidenceMatrix);
  ionMass.clear();
  ionMassBin.clear();
  ionIntens.clear();

  // Find pairs of b ions in 2+ charge state
  ionMass.push_back(nTermMass);
  ionMassBin.push_back(MassConstants::mass2bin(nTermMass));
  ionIntens.push_back(0.0);

  for (int ion = 0; ion < ionMasses.size(); ion++) {
    double tmpIonMass = 2.0 * ionMasses[ion] - massHMono;
    int binTmpIonMass = (int)floor(MassConstants::mass2bin(tmpIonMass));

    ionMass.push_back(tmpIonMass);
    ionMassBin.push_back(binTmpIonMass);
    ionIntens.push_back(ionIntensities[ion]);
  }
  ionMass.push_back(precursorMass - cTermMass);
  ionMassBin.push_back(MassConstants::mass2bin(precursorMass - cTermMass));
  ionIntens.push_back(0.0);

  addEvidToResEvMatrix(ionMass, ionMassBin, ionMasses, ionIntens, ionIntensitiesSort,
                       numSpecPeaks, nAA, aaMass, aaMassBin, residueToleranceMass,
                      residueEvidenceMatrix);
  ionMass.clear();
  ionMassBin.clear();
  ionIntens.clear();

  // Find pairs of y ions in 2+ charge state
  ionMass.push_back(precursorMass - cTermMass);
  ionMassBin.push_back(MassConstants::mass2bin(precursorMass - cTermMass));
  ionIntens.push_back(0.0);

  for (int ion = 0; ion < ionMasses.size(); ion++) {
    double tmpIonMass = precursorMass - (2.0 * ionMasses[ion] - massHMono) + (2.0 * massHMono);
    int binTmpIonMass = (int)floor(MassConstants::mass2bin(tmpIonMass));

    if (tmpIonMass > 0.0) {
      ionMass.push_back(tmpIonMass);
      ionMassBin.push_back(binTmpIonMass);
      ionIntens.push_back(ionIntensities[ion]);
    }
  }
  ionMass.push_back(nTermMass);
  ionMassBin.push_back(MassConstants::mass2bin(nTermMass));
  ionIntens.push_back(0.0);

  reverse(ionMass.begin(), ionMass.end());
  reverse(ionMassBin.begin(),ionMassBin.end());
  reverse(ionIntens.begin(), ionIntens.end());

  addEvidToResEvMatrix(ionMass, ionMassBin, ionMasses, ionIntens, ionIntensitiesSort,
                       numSpecPeaks, nAA, aaMass, aaMassBin, residueToleranceMass,
                       residueEvidenceMatrix);
  ionMass.clear();
  ionMassBin.clear();
  ionIntens.clear();

  // Get maxEvidence value
  double maxEvidence = -1.0;
  for (int i = 0; i < maxPrecurMassBin; i++) {
    for (int curAaMass = 0; curAaMass < nAA; curAaMass++) {
      if (residueEvidenceMatrix[curAaMass][i] > maxEvidence) {
        maxEvidence = residueEvidenceMatrix[curAaMass][i];
      }
    }
  }

  // Discretize residue evidence so largest value is residueEvidenceIntScale
  double residueEvidenceIntScale = (double)granularityScale;
  for (int i = 0; i < maxPrecurMassBin; i++) {
    for (int curAaMass = 0; curAaMass < nAA; curAaMass++) {
      if (residueEvidenceMatrix[curAaMass][i] > 0) {
        double residueEvidence = residueEvidenceMatrix[curAaMass][i];
        residueEvidenceMatrix[curAaMass][i] = round(residueEvidenceIntScale * residueEvidence / maxEvidence);
      }
    }
  }
}
