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

#include <math.h>
#include <stdio.h>
#include <gflags/gflags.h>
#include "spectrum_collection.h"
#include "spectrum_preprocess.h"

using namespace std;

#define NUM_SPECTRUM_REGIONS 10

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
  double partial_sums[end+1];
  for (int i = 0; i < end; ++i)
    partial_sums[i] = (total += observed[i]);
  partial_sums[end] = total;

  for (int i = 0; i < end; ++i) {
    int right_index = min(end, i + MAX_XCORR_OFFSET);
    int left_index = max(0, i - MAX_XCORR_OFFSET - 1);
    observed[i] -= multiplier * (partial_sums[right_index] - partial_sums[left_index]);
  }
}


// From Crux. Not optimized.
static void NormalizeEachRegion(double* observed, int max_mz,
                                double* max_intensity_per_region,
                                int region_selector) {
  int bin_idx = 0;
  int region_idx = 0;
  double max_intensity = max_intensity_per_region[region_idx];
  
  // normalize each region
  for(; bin_idx < max_mz; ++bin_idx){
    if(bin_idx >= region_selector*(region_idx+1) && region_idx < 9){
      ++region_idx;
      max_intensity = max_intensity_per_region[region_idx];;
    }

    // don't normalize if no peaks in region
    if(max_intensity != 0){
      // normalize intensity to max 50
      observed[bin_idx] = (observed[bin_idx] / max_intensity) * 50;
      
#ifdef VERBOSE
      // DEBUG
      printf("bin: %d, region idx: %d, obsered_mz: %.2lf\n", bin_idx,
             region_idx, observed[bin_idx]);
#endif
    }

    // no more peaks beyong the 10 regions mark, exit out
    if(bin_idx > 10*region_selector){
      return;
    }
  }
}

// Mostly from Crux but optimized in some places.
// Postprocessing at end is all new, including these calls:
//  SubtractBackground();
//  MakeInteger();
//  ComputeCache();
void ObservedPeakSet::PreprocessSpectrum(const Spectrum& spectrum,
					 int charge) {
  double peak_location = 0;
  int mz = 0;
  double intensity = 0;
  double bin_width = bin_width_mono;
  double precursor_mz = spectrum.PrecursorMZ();
  double experimental_mass_cut_off = precursor_mz*charge + 50;
  double max_peak_mz = spectrum.M_Z(spectrum.Size()-1);

  assert(MaxMZ::Global().MaxBin() > 0);
  if (experimental_mass_cut_off > FLAGS_max_mz)
    experimental_mass_cut_off = FLAGS_max_mz;

  // Judging from the code below, the highest peak we ever try to set
  // is never more than either experimental_mass_cut_off or
  // max_peak_mz. Things will work as long as max_mz_ is Inited to at
  // least this value.
  max_mz_.Init(min(experimental_mass_cut_off, max_peak_mz));
  // cache_end_ = max_mz_.CacheBinEnd() * NUM_PEAK_TYPES;
  cache_end_ = MaxMZ::Global().CacheBinEnd() * NUM_PEAK_TYPES;

  //memset(peaks_, 0, sizeof(double) * max_mz_.BackgroundBinEnd());
  memset(peaks_, 0, sizeof(double) * MaxMZ::Global().BackgroundBinEnd());
  
  // store the max intensity in each NUM_SPECTRUM_REGIONS regions to later normalize
  double max_intensity_per_region[NUM_SPECTRUM_REGIONS];
  memset(max_intensity_per_region, 0, sizeof(double) * NUM_SPECTRUM_REGIONS);
  int region_selector = int(max_peak_mz / NUM_SPECTRUM_REGIONS);

  int region = 0;
  for (int i = 0; i < spectrum.Size(); ++i) {
    double peak_location = spectrum.M_Z(i);
    
    // skip all peaks larger than experimental mass
    if(peak_location > experimental_mass_cut_off){
      continue;
    }
    
    // skip all peaks within precursor ion mz +/- 15
    if(peak_location < precursor_mz + 15 &&  peak_location > precursor_mz - 15){
      continue;
    }
    
    // map peak location to bin
    mz = (int)(peak_location / bin_width + 0.5);
    region = mz / region_selector;

    // don't let index beyond array
    if(region >= NUM_SPECTRUM_REGIONS){
      continue;
      // region = NUM_SPECTRUM_REGIONS;
    }

    // get intensity
    // sqrt the original intensity
    intensity = sqrt(spectrum.Intensity(i));

    // set intensity in array with correct mz, only if max peak in the bin
    if(peaks_[mz] < intensity) {
      peaks_[mz] = intensity;
            
      // check if this peak is max intensity in the region(one out of 10)
      if(max_intensity_per_region[region] < intensity){
        max_intensity_per_region[region] = intensity;
      }
    }    
  }
  
  // normalize each NUM_SPECTRUM_REGIONS regions to max intensity of 50
  NormalizeEachRegion(peaks_, max_mz_.MaxBin()+1, max_intensity_per_region, region_selector);
  
#ifdef DEBUG
  bool debug = (FLAGS_debug_spectrum_id == spectrum.SpectrumNumber()
                && (FLAGS_debug_charge == 0 || FLAGS_debug_charge == charge));
  if (debug) {
    cout << "GLOBAL MAX MZ: " << MaxMZ::Global().MaxBin() << ", " << MaxMZ::Global().BackgroundBinEnd()
         << ", " << MaxMZ::Global().CacheBinEnd() << endl;
    cout << "MAX MZ: " << max_mz_.MaxBin() << ", " << max_mz_.BackgroundBinEnd()
         << ", " << max_mz_.CacheBinEnd() << endl;
    ShowPeaks();
    cout << "====== SUBTRACTING BACKGROUND ======" << endl;
  }
#endif
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

  for (int i = max_mz_.BackgroundBinEnd() * NUM_PEAK_TYPES; i < cache_end_; ++i)
    cache_[i] = 0;

  for (int i = 0; i < max_mz_.CacheBinEnd(); ++i) {
    int flanks = Peak(PrimaryPeak, i);
    if (i > 0)
      flanks += Peak(FlankingPeak, i-1);
    if (i < max_mz_.CacheBinEnd() - 1)
      flanks += Peak(FlankingPeak, i+1);

    int Y1 = flanks;
    if (i > 16)
      Y1 += Peak(LossPeak, i-17);
    Peak(PeakCombinedY1, i) = Y1;

    int B1 = Y1;
    if (i > 17)
      B1 += Peak(LossPeak, i-18);
    if (i > 27)
      B1 += Peak(LossPeak, i-28);
    Peak(PeakCombinedB1, i) = B1;

    int Y2a = flanks;
    if (i > 8)
      Y2a += Peak(LossPeak, i-9);
    Peak(PeakCombinedY2a, i) = Y2a;

    int Y2b = flanks;
    if (i > 7)
      Y2b += Peak(LossPeak, i-8);
    Peak(PeakCombinedY2b, i) = Y2b;

    int B2a = Y2a;
    if (i > 13)
      B2a += Peak(LossPeak, i-14);
    Peak(PeakCombinedB2a, i) = B2a;

    int B2b = Y2b;
    if (i > 13)
      B2b += Peak(LossPeak, i-14);
    if (i > 8)
      B2b += Peak(LossPeak, i-9);
    Peak(PeakCombinedB2b, i) = B2b;
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

#ifdef DEBUG
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
