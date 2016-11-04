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
#include "util/Params.h"
#include "util/mass.h"
#include <cmath>

#include <iostream> //Remove this when done -- Andy Lin
#include <iomanip> //Remove this when done -- Andy Lin

using namespace std;

// Number of m/z regions in XCorr normalization.
#define NUM_SPECTRUM_REGIONS 10

// Integerization constant for the XCorr p-value calculation.
#define EVIDENCE_INT_SCALE 500.0

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

void ObservedPeakSet::PreprocessSpectrum(const Spectrum& spectrum,
                                       int charge ) {
#ifdef DEBUG
  bool debug = (FLAGS_debug_spectrum_id == spectrum.SpectrumNumber()
                && (FLAGS_debug_charge == 0 || FLAGS_debug_charge == charge));
  if (debug)
    debug = true; // allows a breakpoint
#endif
  double precursor_mz = spectrum.PrecursorMZ();
  double proton = MassConstants::proton;
  double experimental_mass_cut_off = (precursor_mz-proton)*charge+proton + 50;
  double max_peak_mz = spectrum.M_Z(spectrum.Size()-1);

  assert(MaxBin::Global().MaxBinEnd() > 0);

  max_mz_.InitBin(min(experimental_mass_cut_off, max_peak_mz));
  cache_end_ = MaxBin::Global().CacheBinEnd() * NUM_PEAK_TYPES;

  memset(peaks_, 0, sizeof(double) * MaxBin::Global().BackgroundBinEnd());

  if (Params::GetBool("skip-preprocessing")) {
    for (int i = 0; i < spectrum.Size(); ++i) {
      double peak_location = spectrum.M_Z(i);
      if (peak_location >= experimental_mass_cut_off) {
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

    // Fill peaks
    int largest_mz = 0;
    double highest_intensity = 0;
    for (int i = 0; i < spectrum.Size(); ++i) {
      double peak_location = spectrum.M_Z(i);
      if (peak_location >= experimental_mass_cut_off ||
          (remove_precursor && fabs(peak_location - precursor_mz) <= precursor_tolerance)) {
        continue;
      }

      double intensity = spectrum.Intensity(i);
      int mz = MassConstants::mass2bin(peak_location);
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

//TODO Code may be perform same thing as the above PreprocessSpectrum function
//Moved by Andy Lin on 2/12/2016
//Code was originally in CreateEvidenceVector
//Moved preprocessing code to own function
//Preprocesses spectum with following steps
// 1) Discretize mass bins of observed spectrum
// 2) Take square root of intensity
// 3) Grass removal
// 4) Observed spectrum is split into 10 segments and intensities
//    within segment are normalized to maxItensPerRegion (currently set to 50.0)
//Preseved code is below this function
void ObservedPeakSet::PreprocessSpectrum(
  const Spectrum& spectrum,
  double* intensArrayObs,
  int* intensRegion,
  int maxPrecurMass,
  int charge)
{
  // TODO need to review these constants, decide which can be moved to parameter file
  const int nRegion = NUM_SPECTRUM_REGIONS;
  const double maxIntensPerRegion = 50.0;
  const double precursorMZExclude = 15.0; //also used in CreateResidueEvidenceMatrix

  int ma;
  int ionBin;
  double precurMz = spectrum.PrecursorMZ();
  int nIon = spectrum.Size();
  int precurCharge = charge;
  double experimentalMassCutoff = precurMz * precurCharge + 50.0;

  //Determining max ion mass and max ion intensity
  double maxIonMass = 0.0;
  double maxIonIntens = 0.0;
  for (int ion = 0; ion < nIon; ion++) {
    double ionMass = spectrum.M_Z(ion);
    double ionIntens = spectrum.Intensity(ion);
    if (ionMass >= experimentalMassCutoff) {
      continue;
    }
    if (maxIonMass < ionMass) {
      maxIonMass = ionMass;
    }
    if (maxIonIntens < ionIntens) {
      maxIonIntens = ionIntens;
    }
  }
 
  //Mass axis is descretized into 1 dalton bins
  int regionSelector = (int)floor(MassConstants::mass2bin(maxIonMass) / (double)nRegion);
  for (int ion = 0; ion < nIon; ion++) {
    double ionMass = spectrum.M_Z(ion);
    double ionIntens = spectrum.Intensity(ion);
    if (ionMass >= experimentalMassCutoff) {
      continue;
    }
    if (ionMass > precurMz - precursorMZExclude && ionMass < precurMz + precursorMZExclude) {
      continue;
    }
    ionBin = MassConstants::mass2bin(ionMass);
    int region = (int)floor((double)(ionBin) / (double)regionSelector);
    if (region >= nRegion) {
      region = nRegion - 1;
    }
    intensRegion[ionBin] = region;
    if (intensArrayObs[ionBin] < ionIntens) {
      intensArrayObs[ionBin] = ionIntens;
    }
  }

  //Take square root of intensities as well as grass removal
  maxIonIntens = sqrt(maxIonIntens);
  for (ma = 0; ma < maxPrecurMass; ma++) {
    intensArrayObs[ma] = sqrt(intensArrayObs[ma]);
    if (intensArrayObs[ma] <= 0.05 * maxIonIntens) {
      intensArrayObs[ma] = 0.0;
    }
  }

  //Split ovserved spectrum int nRegions (currently 10)
  //Scale each region itensities to maxItensPerRegion (currently 50)
  double* maxRegion = new double[nRegion];
  for (int re = 0; re < nRegion; re++) {
    maxRegion[re] = 0.0;
  }

  for (ma = 0; ma < maxPrecurMass; ma++) {
    int reg = intensRegion[ma];
    if (reg >= 0 && maxRegion[reg] < intensArrayObs[ma]) {
      maxRegion[reg] = intensArrayObs[ma];
    }
  }
  for (ma = 0; ma < maxPrecurMass; ma++) {
    int reg = intensRegion[ma];
    if (reg >= 0 && maxRegion[reg] > 0.0) {
      intensArrayObs[ma] *= (maxIntensPerRegion / maxRegion[reg]);
    }
  }
  delete [] maxRegion;
}
//End -- Moved by Andy Lin

/*
TODO preserved from PreprocessSpectrum only for historical reference; get rid of eventually  
// Fill peaks
int largest_mz = 0;
double highest_intensity = 0;
for (int i = 0; i < spectrum.Size(); ++i) {
  double peak_location = spectrum.M_Z(i);
  if(peak_location >= experimental_mass_cut_off)
    continue;
  double intensity = spectrum.Intensity(i);
  int mz = (int)(peak_location / bin_width + 0.5);
 
  if ((mz > largest_mz) && (intensity > 0))
    largest_mz = mz;
 
  intensity = sqrt(intensity);
  if (intensity > highest_intensity)
    highest_intensity = intensity;
  if (intensity > peaks_[mz])
    peaks_[mz] = intensity;
}

double intensity_cutoff = highest_intensity * 0.05;
int region_size = largest_mz / NUM_SPECTRUM_REGIONS;
for (int i = 0; i < NUM_SPECTRUM_REGIONS; ++i) {
  highest_intensity = 0;
  for (int j = 0; j < region_size; ++j) {
    int index = i * region_size + j;
    if (peaks_[index] <= intensity_cutoff)
      peaks_[index] = 0;
    if (peaks_[index] > highest_intensity)
      highest_intensity = peaks_[index];
  }
  if (highest_intensity == 0)
    continue;
  double normalizer = 50.0/highest_intensity;
  for (int j = 0; j < region_size; ++j) {
    int index = i * region_size + j;
    if (peaks_[index] != 0)
      peaks_[index] *= normalizer;
  }    
}
//TODO end preserved
*/

/* Calculates vector of cleavage evidence for an observed spectrum, using XCorr
 * b/y/neutral peak sets and heights.
 *
 * Written by Jeff Howbert, May, 2013 (as function createEvidenceArrayObserved). 
 * Extended and modified by Jeff Howbert, October, 2013.
 * Ported to and integrated with Tide by Jeff Howbert, November, 2013.
 */
void ObservedPeakSet::CreateEvidenceVector(
  const Spectrum& spectrum,
  double binWidth,
  double binOffset,
  int charge,
  double pepMassMonoMean,
  int maxPrecurMass,
  int* evidenceInt
) {

  // TODO preserved from PreprocessSpectrum; may delete in future
#ifdef DEBUG
  bool debug = (FLAGS_debug_spectrum_id == spectrum.SpectrumNumber()
                && (FLAGS_debug_charge == 0 || FLAGS_debug_charge == charge));
  if (debug) {
    debug = true; // allows a breakpoint
  }
#endif
  // double experimental_mass_cut_off = (precursor_mz-proton)*charge+proton + 50;
  // double max_peak_mz = spectrum.M_Z(spectrum.Size()-1);
  assert(MaxBin::Global().MaxBinEnd() > 0);
  // max_mz_.Init(min(experimental_mass_cut_off, max_peak_mz));
  // TODO end preserved

  // TODO need to review these constants, decide which can be moved to parameter file
  const int nRegion = NUM_SPECTRUM_REGIONS;
  const double maxIntensPerRegion = 50.0;
  const double precursorMZExclude = 15.0;
  const double massNH3Mono = 17.02655;     // mass of NH3 (monoisotopic)
  const double massCOMono =  27.9949;      // mass of CO (monoisotopic)
  const double massH2OMono = 18.010564684; // mass of water (monoisotopic)
  const double massHMono = 1.0078246;      // mass of hydrogen (monoisotopic)
  const double BYHeight = 50.0;
  const double NH3LossHeight = 10.0;    
  const double COLossHeight = 10.0;    // for creating a ions on the fly from b ions
  const double H2OLossHeight = 10.0;
  const double FlankingHeight = BYHeight / 2;;
  // TODO end need to review
  bool flanking_peak = Params::GetBool("use-flanking-peaks");
  bool neutral_loss_peak = Params::GetBool("use-neutral-loss-peaks");

  int ma;
  int pc;
  int ionBin;
  double bIonMass;
  double yIonMass;
  double ionMZMultiCharge;
  double ionMassNH3Loss;
  double ionMassCOLoss;
  double ionMassH2OLoss;

  double* evidence = new double[maxPrecurMass];
  double* intensArrayObs = new double[maxPrecurMass];
  int* intensRegion = new int[maxPrecurMass];

  for (ma = 0; ma < maxPrecurMass; ma++) {
    evidence[ma] = 0.0;
    evidenceInt[ma] = 0;
    intensArrayObs[ma] = 0.0;
    intensRegion[ma] = -1;
  }

  double precurMz = spectrum.PrecursorMZ();
  int nIon = spectrum.Size();
  int precurCharge = charge;
  double experimentalMassCutoff = precurMz * precurCharge + 50.0;
  double proton = MassConstants::proton;

  //Added by Andy Lin on 2/12/2016
  //preprocesses spectrum 
  PreprocessSpectrum(spectrum, intensArrayObs, intensRegion,maxPrecurMass,charge);
  //End added by Andy Lin

  // ***** Adapted from tide/spectrum_preprocess2.cc.
  // TODO replace, if possible, with call to 
  // static void SubtractBackground(double* observed, int end).
  // Note numerous small changes from Tide code.
  double multiplier = 1.0 / (MAX_XCORR_OFFSET * 2.0 + 1.0);
  double total = 0.0;
  double* partial_sums = new double[maxPrecurMass];
  for (int i = 0; i < maxPrecurMass; ++i) {
    partial_sums[ i ] = ( total += intensArrayObs[i]);
  }
  for (int i = 0; i < maxPrecurMass; ++i) {
    int right_index = std::min(maxPrecurMass - 1, i + MAX_XCORR_OFFSET);
    int left_index = std::max(0, i - MAX_XCORR_OFFSET - 1);
    intensArrayObs[i] -= multiplier * (partial_sums[right_index] - partial_sums[left_index]);
  }
  delete [] partial_sums;
  
  int binFirst = MassConstants::mass2bin(30);
  int binLast = MassConstants::mass2bin(pepMassMonoMean - 47);

  for (ma = binFirst; ma <= binLast; ma++) {
    // b ion
    bIonMass = (ma-0.5+binOffset) * binWidth ;
    ionBin = MassConstants::mass2bin(bIonMass);
    evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * BYHeight;
    for (pc = 2; pc < precurCharge; pc++) {
      evidence[ma] = evidence[ma] + intensArrayObs[MassConstants::mass2bin(bIonMass, pc)] * BYHeight;
    }
    // y ion
    yIonMass = pepMassMonoMean + 2 * massHMono - bIonMass;
    ionBin = MassConstants::mass2bin(yIonMass);
    evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * BYHeight;
    for (pc = 2; pc < precurCharge; pc++) {
      evidence[ma] = evidence[ma] + intensArrayObs[MassConstants::mass2bin(yIonMass, pc)] * BYHeight;
    }
    if (flanking_peak == true){

      //flanking peaks for b ions
      ionBin = MassConstants::mass2bin(bIonMass, 1);
      evidence[ma] = evidence[ma] + intensArrayObs[ionBin + 1] * FlankingHeight;
      evidence[ma] = evidence[ma] + intensArrayObs[ionBin - 1] * FlankingHeight;
      for (pc = 2; pc < precurCharge; pc++) {
        evidence[ma] = evidence[ma] + intensArrayObs[MassConstants::mass2bin(bIonMass, pc) + 1] * FlankingHeight;
        evidence[ma] = evidence[ma] + intensArrayObs[MassConstants::mass2bin(bIonMass, pc) - 1] * FlankingHeight;
      }
      //flanking peaks for b ions
      ionBin = MassConstants::mass2bin(yIonMass, pc);
      evidence[ma] = evidence[ma] + intensArrayObs[ionBin + 1] * FlankingHeight;
      evidence[ma] = evidence[ma] + intensArrayObs[ionBin - 1] * FlankingHeight;
      for (pc = 2; pc < precurCharge; pc++) {
        evidence[ma] = evidence[ma] + intensArrayObs[MassConstants::mass2bin(yIonMass, pc) + 1] * FlankingHeight;
        evidence[ma] = evidence[ma] + intensArrayObs[MassConstants::mass2bin(yIonMass, pc) - 1] * FlankingHeight;
      }
    }
    if (neutral_loss_peak == true){
      // NH3 loss from b ion
      ionMassNH3Loss = bIonMass - massNH3Mono;
      ionBin = MassConstants::mass2bin(ionMassNH3Loss);
      evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * NH3LossHeight;
      for (pc = 2; pc < precurCharge; pc++) {
        evidence[ma] = evidence[ma] + intensArrayObs[MassConstants::mass2bin(ionMassNH3Loss, pc)] * NH3LossHeight;
      }
      // NH3 loss from y ion
      ionMassNH3Loss = yIonMass - massNH3Mono;
      ionBin = MassConstants::mass2bin(ionMassNH3Loss);//(int)floor(ionMassNH3Loss / binWidth + 1.0 - binOffset);
      evidence[ma] = evidence[ma] + intensArrayObs[ionBin] * NH3LossHeight;
      for (pc = 2; pc < precurCharge; pc++) {
        evidence[ma] = evidence[ma] + intensArrayObs[MassConstants::mass2bin(ionMassNH3Loss, pc)] * NH3LossHeight;
      }
      // CO and H2O loss from b ion
      ionMassCOLoss = bIonMass - massCOMono;
      ionMassH2OLoss = bIonMass - massH2OMono;
      evidence[ma] = evidence[ma] + intensArrayObs[MassConstants::mass2bin(ionMassCOLoss)] * COLossHeight;
      evidence[ma] = evidence[ma] + intensArrayObs[MassConstants::mass2bin(ionMassH2OLoss)] * H2OLossHeight;
      for (pc = 2; pc < precurCharge; pc++) {
        evidence[ma] = evidence[ma] + intensArrayObs[MassConstants::mass2bin(ionMassCOLoss, pc)] * COLossHeight;
        evidence[ma] = evidence[ma] + intensArrayObs[MassConstants::mass2bin(ionMassH2OLoss, pc)] * H2OLossHeight;
      }
      // H2O loss from y ion
      ionMassH2OLoss = yIonMass - massH2OMono;
      evidence[ma] = evidence[ma] + intensArrayObs[MassConstants::mass2bin(ionMassH2OLoss)] * H2OLossHeight;
      for (pc = 2; pc < precurCharge; pc++) {
        evidence[ma] = evidence[ma] + intensArrayObs[MassConstants::mass2bin(ionMassH2OLoss, pc)] * H2OLossHeight;
      }
    }
  }

  // discretize evidence array
  for (ma = 0; ma < maxPrecurMass; ma++) {
    evidenceInt[ma] = (int)floor(evidence[ma] / EVIDENCE_INT_SCALE + 0.5);
  }
  // clean up
  delete [] evidence;
  delete [] intensArrayObs;
  delete [] intensRegion;

  // TODO preserved from PreprocessSpectrum; may delete in future
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
#ifdef DEBUG
  if (debug)
    ShowPeaks();
#endif
  // TODO end preserved
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

//Written by Andy Lin in Feb 2016
//From an observed spectrum, calculate and populate the residue evidence matrix
void ObservedPeakSet::CreateResidueEvidenceMatrix(
  const Spectrum& spectrum,
  int charge,
  int maxPrecurMassBin,
  double precursorMass,
  int nAA, //TODO different than one used in CreateEvidenceVector
  const vector<double> aaMass, //TODO different aaMass then one used in CreateEvidenceVector
  int fragTol,
  int granularityScale,
  vector<vector<double> >& residueEvidenceMatrix
  ) {

  assert(MaxBin::Global().MaxBinEnd() > 0);

  const double massHMono = MassConstants::mono_h;  // mass of hydrogen (monoisotopic)

  //TODO used in Xcorr function..not used in residue-evidence?
  //bool flanking_peak = Params::GetBool("use-flanking-peaks");
  //bool neutral_loss_peak = Params::GetBool("use-neutral-loss-peaks");

  int ma;
  int pc;

  double precurMz = spectrum.PrecursorMZ();
  int nIon = spectrum.Size();
  int precurCharge = charge;
  double experimentalMassCutoff = precursorMass + 50.0;
//  double proton = MassConstants::proton;
  double residueToleranceMass = 100.0 * (double)fragTol * 0.000001;
  const double precursorMZExclude = 15.0; //also used in PreprocessSpectrum

  //Determining max sqrt(ion intensity)
  double maxIonIntens = 0.0;
  for (int ion = 0; ion < nIon; ion++) {
    double ionMass = spectrum.M_Z(ion);
    double ionIntens = sqrt(spectrum.Intensity(ion));
    if (ionMass >= experimentalMassCutoff) {
      continue;
    }
    if (maxIonIntens < ionIntens) {
      maxIonIntens = ionIntens;
    }
  }

  //MZ that pass preprocessing filters 
  vector<double> ionMasses;
  for(int ion=0 ; ion<nIon ; ion++) {
    double ionMass = spectrum.M_Z(ion);
    double ionIntens = sqrt(spectrum.Intensity(ion));

    //remove peaks larger than experimentalMassCutoff
    if (ionMass >= experimentalMassCutoff) {
      continue;
    }

    //remove peaks within precursorMZExclude of precurMz
    if (ionMass > precurMz - precursorMZExclude && ionMass < precurMz + precursorMZExclude) {
      continue;
    }
    //grass removal
    if (ionIntens < 0.05 * maxIonIntens) {
      continue;
    }
    ionMasses.push_back(ionMass);
  }

  //Determine which bin each amino acid mass is in
  vector<int> aaMassBin;
  for (int i = 0 ; i < nAA ; i++) {
    int binMass= (int)floor(MassConstants::mass2bin(aaMass[i]));
    aaMassBin.push_back(binMass);
  }

  //Need to add lines for matlab line 56?
  //Basically bounds the ion masses and ion mass bins so that
  // 1: ionMassBinB >=1
  // 2: ionMassBinB <= maxMassBin - max(aaMassBin)
  // Also need to figure out maxMassBin (currently hard coded to 8500)
  // need to do 4 times (+1/+2/B ion/ yion)

  int bIonMassBin;
  double bIonMass;
  double bIonIntens;

  //Populate ionMass and ionMassBin for b ions in 1+ charge state
  vector<double> ionMass;
  vector<int> ionMassBin;

  for(int ion=0 ; ion<ionMasses.size() ; ion++) {
    double tmpIonMass = ionMasses[ion];
    int binTmpIonMass = (int)floor(MassConstants::mass2bin(tmpIonMass));

    ionMass.push_back(tmpIonMass);
    ionMassBin.push_back(binTmpIonMass);
  }

  //find pairs of b ions in 1+ charge state
  for (int ion=0 ; ion<ionMasses.size() ; ion++) {
    bIonMass = ionMass[ion];
    bIonMassBin = ionMassBin[ion];

    for(int curAaMass = 0 ; curAaMass < nAA ; curAaMass++) {
      int newResMassBin = bIonMassBin + aaMassBin[curAaMass];

      //Find all ion mass bins that match newResMassBin
      int index = find(ionMassBin.begin(), ionMassBin.end(),newResMassBin) - ionMassBin.begin();
      for(int i = index ; i<ionMasses.size() ; i++) {
        if (newResMassBin != ionMassBin[i])
          continue;

        double ionMassDiff = ionMass[i] - bIonMass;
        double aaTolScore = residueToleranceMass - std::abs(ionMassDiff - aaMass[curAaMass]);

        if (aaTolScore >= 0) {
          //add evidence to matrix
          //&&-1 since all mass bins are index 1 istead of index 0
          residueEvidenceMatrix[curAaMass][newResMassBin-1] += aaTolScore;
        }
      }
    }
  }
  ionMass.clear();
  ionMassBin.clear();

  int yIonMassBin;
  double yIonMass;
  double yIonIntens;
  //find pairs of y ions in 1+ charge state
  for(int ion=0 ; ion<ionMasses.size() ; ion++) {
    //convert to equivalent b ion masses for ease of processing
    double tmpIonMass = precursorMass - ionMasses[ion] + (2.0 * massHMono);
    ionMass.push_back(tmpIonMass);

    //Determine which bin each ion mass is in
    int binTmpIonMass = (int)floor(MassConstants::mass2bin(tmpIonMass));
    ionMassBin.push_back(binTmpIonMass);
  }
  reverse(ionMass.begin(),ionMass.end());
  reverse(ionMassBin.begin(), ionMassBin.end());

  for(int ion=0 ; ion<ionMasses.size() ; ion++) {
    yIonMass = ionMass[ion];
    yIonMassBin = ionMassBin[ion];

    for(int curAaMass=0 ; curAaMass<nAA ; curAaMass++) {
      int newResMassBin = yIonMassBin + aaMassBin[curAaMass];
 
      //Find all ion mass bins that match newResMassBin
      int index = find(ionMassBin.begin(),ionMassBin.end(),newResMassBin) - ionMassBin.begin();
      for(int i=index ; i<ionMasses.size() ; i++) {
        if (newResMassBin != ionMassBin[i])
          continue;
        
        double ionMassDiff = ionMass[i] - yIonMass;
        double aaTolScore = residueToleranceMass - std::abs(ionMassDiff - aaMass[curAaMass]);

        if (aaTolScore >= 0) {
          //add evidence to matrix
          //&&-1 since all mass bins are index 1 istead of index 0
          residueEvidenceMatrix[curAaMass][newResMassBin-1] += aaTolScore;
        }
      }
    }
  }
  ionMass.clear();
  ionMassBin.clear();

  //find pairs of b ions in 2+ charge state
  for(int ion=0 ; ion<ionMasses.size() ; ion++) {
    double tmpIonMass = 2.0 * ionMasses[ion] - massHMono;
    int binTmpIonMass = (int)floor(MassConstants::mass2bin(tmpIonMass));

    ionMass.push_back(tmpIonMass);
    ionMassBin.push_back(binTmpIonMass);
  }

  for(int ion=0 ; ion<ionMasses.size() ; ion++) {
    bIonMass = ionMass[ion];
    bIonMassBin = ionMassBin[ion];

    for(int curAaMass=0 ; curAaMass<nAA ; curAaMass++) {
      int newResMassBin = bIonMassBin + aaMassBin[curAaMass];
      
      //Find all ion mass bins that match newResMassBin
      int index = find(ionMassBin.begin(),ionMassBin.end(),newResMassBin) - ionMassBin.begin();
      for(int i=index ; i<ionMasses.size() ; i++) {
        if (newResMassBin != ionMassBin[i])
          continue;

        double ionMassDiff = ionMass[i] - bIonMass;
        double aaTolScore = residueToleranceMass - std::abs(ionMassDiff - aaMass[curAaMass]);

        if(aaTolScore >= 0) {
          //add evidence to matrix
          //&&-1 since all mass bins are index 1 instead of index 0
          
          residueEvidenceMatrix[curAaMass][newResMassBin-1] += aaTolScore;
        }
      }
    }
  }
  ionMass.clear();
  ionMassBin.clear();

  //find pairs of y ions in 2+ charge state
  for(int ion=0 ; ion<ionMasses.size() ; ion++) {
    double tmpIonMass = precursorMass - (2.0 * ionMasses[ion] - massHMono) + (2.0 * massHMono);
    int binTmpIonMass = (int)floor(MassConstants::mass2bin(tmpIonMass));

    ionMass.push_back(tmpIonMass);
    ionMassBin.push_back(binTmpIonMass);
  }
  reverse(ionMass.begin(),ionMass.end());
  reverse(ionMassBin.begin(),ionMassBin.end());

  for(int ion=0 ; ion<ionMasses.size() ; ion++) {
    yIonMass = ionMass[ion];
    yIonMassBin = ionMassBin[ion];

    for(int curAaMass=0 ; curAaMass<nAA ; curAaMass++) {
      int newResMassBin = yIonMassBin + aaMassBin[curAaMass];
 
      //Find all ion mass bins that match newResMassBin
      int index = find(ionMassBin.begin(), ionMassBin.end(), newResMassBin) - ionMassBin.begin();
      for(int i=index ; i<ionMasses.size() ; i++) {
        if(newResMassBin != ionMassBin[i])
          continue;
   
        double ionMassDiff = ionMass[i] - yIonMass;
        double aaTolScore = residueToleranceMass - std::abs(ionMassDiff - aaMass[curAaMass]);

        if (aaTolScore >= 0) {
          //add evidence to matrix
          //&&-1 since all mass bins are index 1 instead of index 0
          residueEvidenceMatrix[curAaMass][newResMassBin-1] += aaTolScore;
        }
      }
    }
  }
  ionMass.clear();
  ionMassBin.clear();

  //Get maxEvidence value
  double maxEvidence = -1.0;
  for(int i =0 ; i<maxPrecurMassBin ; i++) {
    for(int curAaMass = 0; curAaMass<nAA;curAaMass++){
      if (residueEvidenceMatrix[curAaMass][i] > maxEvidence) {
        maxEvidence = residueEvidenceMatrix[curAaMass][i];
      }
    }
  }

  //discretize residue evidence so largest value is residueEvidenceIntScale
  double residueEvidenceIntScale = (double)granularityScale;
  for(int i=0 ; i<maxPrecurMassBin ; i++) {
    for(int curAaMass = 0; curAaMass<nAA;curAaMass++){
      if (residueEvidenceMatrix[curAaMass][i] > 0) {
        double residueEvidence = residueEvidenceMatrix[curAaMass][i];
        residueEvidenceMatrix[curAaMass][i] = round(residueEvidenceIntScale * residueEvidence / maxEvidence);
      }
    }
  }

//  std::cout << "maxEvidence: " << maxEvidence << std::endl;
/*
  int cnt=0;
  for(int i =0;i<maxPrecurMassBin;i++) {
    for(int curAaMass = 0; curAaMass<nAA;curAaMass++){
      if (residueEvidenceMatrix[curAaMass][i] != 0) {
//        std::cout << residueEvidenceMatrix[curAaMass][i] << std::endl;
        std::cout << "aa: " << curAaMass+1 << endl;
        std::cout << "massBin: " << i+1 << std::endl;
        std::cout << "evidence: " << residueEvidenceMatrix[curAaMass][i] << std::endl <<std::endl;
        cnt+=1;
      }
    }
  }
*/
}
